struct SegmentCache
    nominal_speeds::Vector{Float64}
    nominal_drops::Vector{Float64}
    junction_pressure::Float64
end

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

function segment_pressure_drops(segments::Vector{BranchSegment}, flow::FlowParameters, speeds::Vector{Float64})
    q = [
        speeds[1] * π * segments[1].radius^2,
        speeds[2] * π * segments[2].radius^2,
        speeds[3] * π * segments[3].radius^2,
    ]
    drops = [8.0 * flow.viscosity * q[i] * segments[i].length / (π * segments[i].radius^4) for i in eachindex(segments)]
    junction_pressure = flow.outlet_pressure + max(drops[2], drops[3])
    return junction_pressure, drops
end

function build_segment_cache(grid::GeometryGrid, params::SimulationParameters)
    nominal_speeds = segment_nominal_speeds(grid.segments, params.flow)
    junction_pressure, nominal_drops = segment_pressure_drops(grid.segments, params.flow, nominal_speeds)
    return SegmentCache(nominal_speeds, nominal_drops, junction_pressure)
end

function allocate_flow_state(grid::GeometryGrid, params::SimulationParameters)
    nx, ny, nz = size(grid.mask)
    return FlowState(
        zeros(Float64, nx, ny, nz, 3),
        zeros(Float64, nx, ny, nz),
        fill(params.flow.outlet_pressure, nx, ny, nz),
        zeros(Float64, nx, ny, nz),
        zeros(Float64, nx, ny, nz),
    )
end

function compute_flow_field!(
    state::FlowState,
    grid::GeometryGrid,
    phi::Array{Float64, 3},
    params::SimulationParameters,
    segment_cache::SegmentCache,
    active_indices::Vector{CartesianIndex{3}},
)
    fill!(state.velocity, 0.0)
    fill!(state.speed, 0.0)
    fill!(state.pressure, params.flow.outlet_pressure)
    fill!(state.shear, 0.0)
    fill!(state.fluid_fraction, 0.0)

    @inbounds for idx in active_indices
        i, j, k = Tuple(idx)
        seg_id = grid.segment_id[idx]
        segment = grid.segments[seg_id]
        nominal_speed = segment_cache.nominal_speeds[seg_id]
        open_fraction = clamp(
            1.0 - params.clot.occlusion_sensitivity * phi[idx],
            params.flow.minimum_open_fraction,
            1.0,
        )
        state.fluid_fraction[idx] = open_fraction
        effective_radius = segment.radius * sqrt(open_fraction)
        radial = grid.centerline_distance[idx]

        local_drop = segment_cache.nominal_drops[seg_id] / max(open_fraction^2, params.flow.minimum_open_fraction^2)
        axial = grid.axial_coordinate[idx]
        state.pressure[idx] = if seg_id == 1
            segment_cache.junction_pressure + local_drop * (1.0 - axial / segment.length)
        else
            params.flow.outlet_pressure + local_drop * (1.0 - axial / segment.length)
        end

        if radial >= effective_radius
            continue
        end

        speed_multiplier = min(params.flow.max_velocity_multiplier, open_fraction^(-0.5))
        mean_speed = nominal_speed * speed_multiplier
        local_speed = 2.0 * mean_speed * max(0.0, 1.0 - (radial / effective_radius)^2)
        state.speed[idx] = local_speed
        state.velocity[i, j, k, 1] = local_speed * segment.direction[1]
        state.velocity[i, j, k, 2] = local_speed * segment.direction[2]
        state.velocity[i, j, k, 3] = local_speed * segment.direction[3]
        state.shear[idx] = 4.0 * params.flow.viscosity * mean_speed / max(effective_radius, eps())
    end

    return state
end

function compute_flow_field(grid::GeometryGrid, phi::Array{Float64, 3}, params::SimulationParameters)
    state = allocate_flow_state(grid, params)
    segment_cache = build_segment_cache(grid, params)
    active_indices = findall(grid.mask)
    return compute_flow_field!(state, grid, phi, params, segment_cache, active_indices)
end
