struct BranchSegment
    start::SVector{3, Float64}
    direction::SVector{3, Float64}
    length::Float64
    radius::Float64
    label::Symbol
end

struct GeometryGrid
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    dx::Float64
    dy::Float64
    dz::Float64
    mask::BitArray{3}
    wall_distance::Array{Float64, 3}
    centerline_distance::Array{Float64, 3}
    radius_field::Array{Float64, 3}
    axial_coordinate::Array{Float64, 3}
    segment_id::Array{Int, 3}
    segments::Vector{BranchSegment}
    bounds::NTuple{6, Float64}
end

function build_segments(params::GeometryParameters)
    θ = 0.5 * params.bifurcation_angle
    inlet_radius = 0.5 * params.inlet_diameter
    branch_radii = 0.5 .* params.branch_diameters
    inlet = BranchSegment(
        @SVector [0.0, 0.0, -params.inlet_length],
        @SVector [0.0, 0.0, 1.0],
        params.inlet_length,
        inlet_radius,
        :inlet,
    )
    upper = BranchSegment(
        @SVector [0.0, 0.0, 0.0],
        normalize(@SVector [sin(θ), 0.0, cos(θ)]),
        params.branch_length,
        branch_radii[1],
        :branch_a,
    )
    lower = BranchSegment(
        @SVector [0.0, 0.0, 0.0],
        normalize(@SVector [-sin(θ), 0.0, cos(θ)]),
        params.branch_length,
        branch_radii[2],
        :branch_b,
    )
    return [inlet, upper, lower]
end

function domain_bounds(params::GeometryParameters)
    θ = 0.5 * params.bifurcation_angle
    max_radius = 0.5 * max(params.inlet_diameter, maximum(params.branch_diameters))
    x_extent = params.branch_length * sin(θ) + max_radius + params.padding
    y_extent = max_radius + params.padding
    zmin = -params.inlet_length - max_radius - params.padding
    zmax = params.branch_length * cos(θ) + max_radius + params.padding
    return (-x_extent, x_extent, -y_extent, y_extent, zmin, zmax)
end

function capsule_coordinates(point::SVector{3, Float64}, segment::BranchSegment)
    relative = point - segment.start
    axial = clamp(dot(relative, segment.direction), 0.0, segment.length)
    closest = segment.start + axial * segment.direction
    radial = norm(point - closest)
    wall_distance = segment.radius - radial
    return wall_distance, radial, axial
end

function build_bifurcation_grid(params::GeometryParameters, numerics::NumericsParameters = NumericsParameters())
    bounds = domain_bounds(params)
    x = collect(range(bounds[1], bounds[2], length = numerics.nx))
    y = collect(range(bounds[3], bounds[4], length = numerics.ny))
    z = collect(range(bounds[5], bounds[6], length = numerics.nz))
    dx = x[2] - x[1]
    dy = y[2] - y[1]
    dz = z[2] - z[1]
    segments = build_segments(params)

    mask = falses(length(x), length(y), length(z))
    wall_distance = fill(-Inf, length(x), length(y), length(z))
    centerline_distance = fill(Inf, length(x), length(y), length(z))
    radius_field = zeros(length(x), length(y), length(z))
    axial_coordinate = zeros(length(x), length(y), length(z))
    segment_id = zeros(Int, length(x), length(y), length(z))

    @inbounds for k in eachindex(z), j in eachindex(y), i in eachindex(x)
        point = @SVector [x[i], y[j], z[k]]
        best_wall = -Inf
        best_radial = Inf
        best_axial = 0.0
        best_id = 0
        best_radius = 0.0
        for (id, segment) in enumerate(segments)
            wd, radial, axial = capsule_coordinates(point, segment)
            if wd > best_wall
                best_wall = wd
                best_radial = radial
                best_axial = axial
                best_id = id
                best_radius = segment.radius
            end
        end
        if best_wall > 0.0
            mask[i, j, k] = true
            wall_distance[i, j, k] = best_wall
            centerline_distance[i, j, k] = best_radial
            axial_coordinate[i, j, k] = best_axial
            segment_id[i, j, k] = best_id
            radius_field[i, j, k] = best_radius
        end
    end

    return GeometryGrid(
        x,
        y,
        z,
        dx,
        dy,
        dz,
        mask,
        wall_distance,
        centerline_distance,
        radius_field,
        axial_coordinate,
        segment_id,
        segments,
        bounds,
    )
end

function orthonormal_basis(direction::SVector{3, Float64})
    reference = abs(direction[3]) < 0.9 ? @SVector [0.0, 0.0, 1.0] : @SVector [0.0, 1.0, 0.0]
    n1 = normalize(cross(direction, reference))
    n2 = normalize(cross(direction, n1))
    return n1, n2
end

function tube_mesh(segment::BranchSegment; nθ::Int = 40, nz::Int = 28)
    n1, n2 = orthonormal_basis(segment.direction)
    vertices = Point3f[]
    faces = TriangleFace{Int}[]

    for iz in 0:nz
        s = segment.length * iz / nz
        center = segment.start + s * segment.direction
        for iθ in 0:(nθ - 1)
            θ = 2π * iθ / nθ
            point = center + segment.radius * (cos(θ) * n1 + sin(θ) * n2)
            push!(vertices, Point3f(point[1], point[2], point[3]))
        end
    end

    ring = nθ
    for iz in 0:(nz - 1)
        offset = iz * ring
        for iθ in 1:ring
            i1 = offset + iθ
            i2 = offset + (iθ % ring) + 1
            i3 = offset + ring + iθ
            i4 = offset + ring + (iθ % ring) + 1
            push!(faces, TriangleFace(i1, i2, i3))
            push!(faces, TriangleFace(i2, i4, i3))
        end
    end

    return GeometryBasics.Mesh(vertices, faces)
end

surface_meshes(params::GeometryParameters; nθ::Int = 40, nz::Int = 28) =
    [tube_mesh(segment; nθ = nθ, nz = nz) for segment in build_segments(params)]

cell_volume(grid::GeometryGrid) = grid.dx * grid.dy * grid.dz

function cell_centers(grid::GeometryGrid)
    points = Point3f[]
    @inbounds for k in eachindex(grid.z), j in eachindex(grid.y), i in eachindex(grid.x)
        if grid.mask[i, j, k]
            push!(points, Point3f(grid.x[i], grid.y[j], grid.z[k]))
        end
    end
    return points
end
