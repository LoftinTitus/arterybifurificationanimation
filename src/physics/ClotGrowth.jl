initialize_clot(grid::GeometryGrid) = zeros(Float64, size(grid.mask))

function deposition_response(τ::Float64, params::SimulationParameters)
    ratio = (τ / params.clot.tau_dep)^params.clot.dep_hill
    return 1.0 / (1.0 + ratio)
end

function erosion_response(τ::Float64, params::SimulationParameters)
    ratio = (τ / params.clot.tau_ero)^params.clot.ero_hill
    return ratio / (1.0 + ratio)
end

function wall_weight(distance::Float64, params::SimulationParameters)
    return exp(-(distance / params.clot.wall_bandwidth)^2)
end

function update_clot(
    φ::Array{Float64, 3},
    C::Array{Float64, 3},
    flow::FlowState,
    grid::GeometryGrid,
    params::SimulationParameters,
    dt::Float64,
)
    φnew = copy(φ)

    @inbounds for k in eachindex(grid.z), j in eachindex(grid.y), i in eachindex(grid.x)
        if !grid.mask[i, j, k]
            φnew[i, j, k] = 0.0
            continue
        end

        w = wall_weight(grid.wall_distance[i, j, k], params)
        fτ = deposition_response(flow.shear[i, j, k], params)
        gτ = erosion_response(flow.shear[i, j, k], params)
        growth = params.clot.k_dep * C[i, j, k] * w * fτ * (1.0 - φ[i, j, k])
        erosion = params.clot.k_ero * gτ * φ[i, j, k]
        φnew[i, j, k] = clamp(φ[i, j, k] + dt * (growth - erosion), 0.0, 1.0)
    end

    φnew[.!grid.mask] .= 0.0
    return φnew
end
