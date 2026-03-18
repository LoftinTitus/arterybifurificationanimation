struct SimulationOutputs
    params::SimulationParameters
    grid::GeometryGrid
    times::Vector{Float64}
    clot_volume::Vector{Float64}
    occlusion::Vector{Float64}
    mean_wall_shear::Vector{Float64}
    snapshot_times::Vector{Float64}
    clot_snapshots::Vector{Array{Float32, 3}}
    final_platelets::Array{Float32, 3}
    final_activated::Array{Float32, 3}
    final_agonist::Array{Float32, 3}
    final_bound::Array{Float32, 3}
    final_shear::Array{Float32, 3}
    final_velocity::Array{Float32, 4}
    final_pressure::Array{Float32, 3}
end

function compute_metrics(φ::Array{Float64, 3}, flow::FlowState, grid::GeometryGrid, cache::GridIndexCache)
    dV = cell_volume(grid)
    clot_sum = 0.0
    open_sum = 0.0
    @inbounds for idx in cache.active
        clot_sum += φ[idx]
        open_sum += flow.fluid_fraction[idx]
    end

    wall_shear_sum = 0.0
    @inbounds for idx in cache.wall
        wall_shear_sum += flow.shear[idx]
    end

    clot_volume = clot_sum * dV
    occlusion = 100.0 * (1.0 - open_sum / cache.active_count)
    mean_shear = isempty(cache.wall) ? 0.0 : wall_shear_sum / length(cache.wall)
    return clot_volume, occlusion, mean_shear
end

function snapshot!(store::Vector{Array{Float32, 3}}, field::Array{Float64, 3})
    push!(store, Float32.(field))
    return nothing
end

function total_platelets_snapshot(state::TransportState)
    total = Array{Float32, 3}(undef, size(state.resting))
    @inbounds for idx in eachindex(total)
        total[idx] = Float32(state.resting[idx] + state.activated[idx])
    end
    return total
end

function run_simulation(
    params::SimulationParameters = default_parameters();
    show_progress::Bool = false,
    progress_interval::Float64 = 0.05,
)
    grid = build_bifurcation_grid(params.geometry, params.numerics)
    cache = build_index_cache(grid, params)
    segment_cache = build_segment_cache(grid, params)
    transport = initialize_transport_state(grid, params)
    transport_next = allocate_transport_state(grid)
    clot = initialize_clot_state(grid)
    clot_next = allocate_clot_state(grid)
    flow = allocate_flow_state(grid, params)
    compute_flow_field!(flow, grid, clot.solid, params, segment_cache, cache.active)

    times = Float64[]
    clot_volume = Float64[]
    occlusion = Float64[]
    mean_wall_shear = Float64[]
    snapshot_times = Float64[]
    clot_snapshots = Array{Float32, 3}[]

    t = 0.0
    step_id = 0
    next_progress_report = progress_interval
    start_time = time()
    if show_progress
        println("Simulation progress: 0.0% (t = 0.00 / $(round(params.numerics.tfinal, digits = 2)) s)")
    end
    while t < params.numerics.tfinal
        dt = min(stable_timestep(flow, grid, params), params.numerics.tfinal - t)
        update_transport!(transport_next, transport, clot, flow, grid, params, cache, dt)
        update_clot!(clot_next, clot, transport_next, flow, grid, params, cache, dt)
        transport, transport_next = transport_next, transport
        clot, clot_next = clot_next, clot
        compute_flow_field!(flow, grid, clot.solid, params, segment_cache, cache.active)
        t += dt
        step_id += 1

        cv, occ, τm = compute_metrics(clot.solid, flow, grid, cache)
        push!(times, t)
        push!(clot_volume, cv)
        push!(occlusion, occ)
        push!(mean_wall_shear, τm)

        if step_id == 1 || step_id % params.numerics.save_every == 0 || t >= params.numerics.tfinal
            push!(snapshot_times, t)
            snapshot!(clot_snapshots, clot.solid)
        end

        if show_progress
            progress = t / params.numerics.tfinal
            if progress >= next_progress_report || t >= params.numerics.tfinal
                elapsed = time() - start_time
                println(
                    "Simulation progress: $(round(100 * progress, digits = 1))% " *
                    "(t = $(round(t, digits = 2)) / $(round(params.numerics.tfinal, digits = 2)) s, " *
                    "step = $(step_id), elapsed = $(round(elapsed, digits = 1)) s)",
                )
                next_progress_report += progress_interval
            end
        end
    end

    return SimulationOutputs(
        params,
        grid,
        times,
        clot_volume,
        occlusion,
        mean_wall_shear,
        snapshot_times,
        clot_snapshots,
        total_platelets_snapshot(transport),
        Float32.(transport.activated),
        Float32.(transport.agonist),
        Float32.(clot.bound),
        Float32.(flow.shear),
        Float32.(flow.velocity),
        Float32.(flow.pressure),
    )
end

function write_metrics_csv(results::SimulationOutputs, filepath::AbstractString)
    header = ["time_s", "clot_volume_m3", "occlusion_percent", "mean_wall_shear_pa"]
    data = hcat(results.times, results.clot_volume, results.occlusion, results.mean_wall_shear)
    open(filepath, "w") do io
        writedlm(io, permutedims(header), ',')
        writedlm(io, data, ',')
    end
    return filepath
end
