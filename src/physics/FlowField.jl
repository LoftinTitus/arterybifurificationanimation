struct FlowState
    velocity::Array{Float64, 4}
    speed::Array{Float64, 3}
    pressure::Array{Float64, 3}
    shear::Array{Float64, 3}
    fluid_fraction::Array{Float64, 3}
end

function segment_nominal_speeds(segments::Vector{BranchSegment}, flow::FlowParameters)
    inlet = segments[1]
    branches = segments[2:3]
    inlet_area = π * inlet.radius^2
    q_in = flow.inlet_mean_velocity * inlet_area
    conductances = [branch.radius^4 / branch.length for branch in branches]
    conductance_sum = sum(conductances)
    branch_flows = q_in .* conductances ./ conductance_sum
    branch_speeds = branch_flows ./ (π .* getfield.(branches, :radius) .^ 2)
    return [flow.inlet_mean_velocity, branch_speeds[1], branch_speeds[2]]
end

function segment_base_pressures(segments::Vector{BranchSegment}, flow::FlowParameters, speeds::Vector{Float64})
    q = [
        speeds[1] * π * segments[1].radius^2,
        speeds[2] * π * segments[2].radius^2,
        speeds[3] * π * segments[3].radius^2,
    ]
    drops = [8.0 * flow.viscosity * q[i] * segments[i].length / (π * segments[i].radius^4) for i in eachindex(segments)]
    junction_pressure = flow.outlet_pressure + max(drops[2], drops[3])
    return [junction_pressure + drops[1], flow.outlet_pressure + drops[2], flow.outlet_pressure + drops[3]], drops
end

function compute_flow_field(grid::GeometryGrid, phi::Array{Float64, 3}, params::SimulationParameters)
    nx, ny, nz = size(phi)
    velocity = zeros(Float64, nx, ny, nz, 3)
    speed = zeros(Float64, nx, ny, nz)
    pressure = fill(params.flow.outlet_pressure, nx, ny, nz)
    shear = zeros(Float64, nx, ny, nz)
    fluid_fraction = zeros(Float64, nx, ny, nz)

    nominal_speeds = segment_nominal_speeds(grid.segments, params.flow)
    base_pressures, nominal_drops = segment_base_pressures(grid.segments, params.flow, nominal_speeds)

    @inbounds for k in eachindex(grid.z), j in eachindex(grid.y), i in eachindex(grid.x)
        if !grid.mask[i, j, k]
            continue
        end

        seg_id = grid.segment_id[i, j, k]
        segment = grid.segments[seg_id]
        nominal_speed = nominal_speeds[seg_id]
        open_fraction = clamp(
            1.0 - params.clot.occlusion_sensitivity * phi[i, j, k],
            params.flow.minimum_open_fraction,
            1.0,
        )
        fluid_fraction[i, j, k] = open_fraction
        effective_radius = segment.radius * sqrt(open_fraction)
        radial = grid.centerline_distance[i, j, k]

        if radial >= effective_radius
            continue
        end

        speed_multiplier = min(params.flow.max_velocity_multiplier, open_fraction^(-0.5))
        mean_speed = nominal_speed * speed_multiplier
        local_speed = 2.0 * mean_speed * max(0.0, 1.0 - (radial / effective_radius)^2)
        speed[i, j, k] = local_speed
        velocity[i, j, k, 1] = local_speed * segment.direction[1]
        velocity[i, j, k, 2] = local_speed * segment.direction[2]
        velocity[i, j, k, 3] = local_speed * segment.direction[3]

        local_drop = nominal_drops[seg_id] / max(open_fraction^2, params.flow.minimum_open_fraction^2)
        pressure[i, j, k] = if seg_id == 1
            junction_pressure = base_pressures[1] - nominal_drops[1]
            junction_pressure + local_drop * (1.0 - grid.axial_coordinate[i, j, k] / segment.length)
        else
            params.flow.outlet_pressure + local_drop * (1.0 - grid.axial_coordinate[i, j, k] / segment.length)
        end

        shear[i, j, k] = 4.0 * params.flow.viscosity * mean_speed / max(effective_radius, eps())
    end

    return FlowState(velocity, speed, pressure, shear, fluid_fraction)
end
