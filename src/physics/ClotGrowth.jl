initialize_clot(grid::GeometryGrid) = zeros(Float64, size(grid.mask))

function initialize_clot_state(grid::GeometryGrid)
    nx, ny, nz = size(grid.mask)
    return ClotState(
        zeros(Float64, nx, ny, nz),
        zeros(Float64, nx, ny, nz),
    )
end

function allocate_clot_state(grid::GeometryGrid)
    nx, ny, nz = size(grid.mask)
    return ClotState(
        zeros(Float64, nx, ny, nz),
        zeros(Float64, nx, ny, nz),
    )
end

deposition_response(τ::Float64, params::SimulationParameters) = low_shear_switch(τ, params)
erosion_response(τ::Float64, params::SimulationParameters) = erosion_switch(τ, params)

function update_clot!(
    next::ClotState,
    current::ClotState,
    transport::TransportState,
    flow::FlowState,
    grid::GeometryGrid,
    params::SimulationParameters,
    cache::GridIndexCache,
    dt::Float64,
)
    bound = current.bound
    φ = current.solid
    bound_next = next.bound
    φ_next = next.solid
    fill!(bound_next, 0.0)
    fill!(φ_next, 0.0)

    @inbounds for idx in cache.active
        τ = flow.shear[idx]
        w = cache.wall_weight[idx]
        pa = transport.activated[idx]
        local_bound = bound[idx]
        solid = φ[idx]
        bind = binding_rate(pa, local_bound, solid, w, τ, params)
        detach = bound_detachment_rate(local_bound, τ, params)
        compact = compaction_rate(local_bound, solid, w, params)
        erosion = params.clot.k_ero * erosion_response(τ, params) * solid
        bound_next[idx] = clamp(local_bound + dt * (bind - detach - compact), 0.0, params.clot.bound_capacity)
        φ_next[idx] = clamp(solid + dt * (compact - erosion), 0.0, 1.0)
    end

    return next
end
