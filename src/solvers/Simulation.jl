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

function wall_mask(grid::GeometryGrid, params::SimulationParameters)
    return grid.mask .& (grid.wall_distance .<= 2.0 * params.clot.wall_bandwidth)
end

function compute_metrics(φ::Array{Float64, 3}, flow::FlowState, grid::GeometryGrid, params::SimulationParameters)
    dV = cell_volume(grid)
    baseline_volume = sum(grid.mask) * dV
    clot_volume = sum(φ .* Float64.(grid.mask)) * dV
    open_volume = sum(flow.fluid_fraction .* Float64.(grid.mask)) * dV
    occlusion = 100.0 * (1.0 - open_volume / baseline_volume)
    wm = wall_mask(grid, params)
    mean_shear = any(wm) ? mean(flow.shear[wm]) : 0.0
    return clot_volume, occlusion, mean_shear
end

function snapshot!(store::Vector{Array{Float32, 3}}, field::Array{Float64, 3})
    push!(store, Float32.(field))
    return nothing
end

function run_simulation(params::SimulationParameters = default_parameters())
    grid = build_bifurcation_grid(params.geometry, params.numerics)
    C = initialize_platelets(grid, params)
    φ = initialize_clot(grid)
    flow = compute_flow_field(grid, φ, params)

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
        flow = compute_flow_field(grid, φ, params)
        dt = min(stable_timestep(flow, grid, params), params.numerics.tfinal - t)
        C = update_platelets(C, φ, flow, grid, params, dt)
        φ = update_clot(φ, C, flow, grid, params, dt)
        t += dt
        step_id += 1

        flow = compute_flow_field(grid, φ, params)
        cv, occ, τm = compute_metrics(φ, flow, grid, params)
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
