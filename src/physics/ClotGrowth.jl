initialize_clot(grid::GeometryGrid) = zeros(Float64, size(grid.mask))

function deposition_response(τ::Float64, params::SimulationParameters)
    ratio = (τ / params.clot.tau_dep)^params.clot.dep_hill
    return 1.0 / (1.0 + ratio)
end

function erosion_response(τ::Float64, params::SimulationParameters)
    ratio = (τ / params.clot.tau_ero)^params.clot.ero_hill
    return ratio / (1.0 + ratio)
end

function update_clot!(
    φnew::Array{Float64, 3},
    φ::Array{Float64, 3},
    C::Array{Float64, 3},
    flow::FlowState,
    grid::GeometryGrid,
    params::SimulationParameters,
    cache::GridIndexCache,
    dt::Float64,
)
    fill!(φnew, 0.0)

    @inbounds for idx in cache.active
        w = cache.wall_weight[idx]
        fτ = deposition_response(flow.shear[idx], params)
        gτ = erosion_response(flow.shear[idx], params)
        growth = params.clot.k_dep * C[idx] * w * fτ * (1.0 - φ[idx])
        erosion = params.clot.k_ero * gτ * φ[idx]
        φnew[idx] = clamp(φ[idx] + dt * (growth - erosion), 0.0, 1.0)
    end

    return φnew
end
