struct SimulationOutputs
    params::SimulationParameters
    grid::GeometryGrid
    times::Vector{Float64}
    clot_volume::Vector{Float64}
    occlusion::Vector{Float64}
    mean_wall_shear::Vector{Float64}
    snapshot_times::Vector{Float64}
    platelet_snapshots::Vector{Array{Float32, 3}}
    clot_snapshots::Vector{Array{Float32, 3}}
    shear_snapshots::Vector{Array{Float32, 3}}
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

function run_simulation(params::SimulationParameters = default_parameters())
    grid = build_bifurcation_grid(params.geometry, params.numerics)
    cache = build_index_cache(grid, params)
    segment_cache = build_segment_cache(grid, params)
    C = initialize_platelets(grid, params)
    Cnext = similar(C)
    φ = initialize_clot(grid)
    φnext = similar(φ)
    flow = allocate_flow_state(grid, params)
    compute_flow_field!(flow, grid, φ, params, segment_cache, cache.active)

    times = Float64[]
    clot_volume = Float64[]
    occlusion = Float64[]
    mean_wall_shear = Float64[]
    snapshot_times = Float64[]
    platelet_snapshots = Array{Float32, 3}[]
    clot_snapshots = Array{Float32, 3}[]
    shear_snapshots = Array{Float32, 3}[]

    t = 0.0
    step_id = 0
    while t < params.numerics.tfinal
        dt = min(stable_timestep(flow, grid, params), params.numerics.tfinal - t)
        update_platelets!(Cnext, C, φ, flow, grid, params, cache, dt)
        update_clot!(φnext, φ, Cnext, flow, grid, params, cache, dt)
        C, Cnext = Cnext, C
        φ, φnext = φnext, φ
        compute_flow_field!(flow, grid, φ, params, segment_cache, cache.active)
        t += dt
        step_id += 1

        cv, occ, τm = compute_metrics(φ, flow, grid, cache)
        push!(times, t)
        push!(clot_volume, cv)
        push!(occlusion, occ)
        push!(mean_wall_shear, τm)

        if step_id == 1 || step_id % params.numerics.save_every == 0 || t >= params.numerics.tfinal
            push!(snapshot_times, t)
            snapshot!(platelet_snapshots, C)
            snapshot!(clot_snapshots, φ)
            snapshot!(shear_snapshots, flow.shear)
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
        platelet_snapshots,
        clot_snapshots,
        shear_snapshots,
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
