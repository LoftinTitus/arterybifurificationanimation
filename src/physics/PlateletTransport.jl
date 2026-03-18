struct TransportState
    resting::Array{Float64, 3}
    activated::Array{Float64, 3}
    agonist::Array{Float64, 3}
end

struct ClotState
    bound::Array{Float64, 3}
    solid::Array{Float64, 3}
end

function initialize_transport_state(grid::GeometryGrid, params::SimulationParameters)
    nx, ny, nz = size(grid.mask)
    resting = zeros(Float64, nx, ny, nz)
    activated = zeros(Float64, nx, ny, nz)
    agonist = zeros(Float64, nx, ny, nz)

    inlet_activated = clamp(params.transport.inlet_activated_fraction, 0.0, 1.0) * params.transport.inlet_concentration
    inlet_resting = params.transport.inlet_concentration - inlet_activated

    @inbounds for idx in eachindex(grid.mask)
        if grid.mask[idx]
            resting[idx] = inlet_resting
            activated[idx] = inlet_activated
            agonist[idx] = params.transport.inlet_agonist
        end
    end

    return TransportState(resting, activated, agonist)
end

initialize_platelets(grid::GeometryGrid, params::SimulationParameters) =
    initialize_transport_state(grid, params).resting

function allocate_transport_state(grid::GeometryGrid)
    nx, ny, nz = size(grid.mask)
    return TransportState(
        zeros(Float64, nx, ny, nz),
        zeros(Float64, nx, ny, nz),
        zeros(Float64, nx, ny, nz),
    )
end

@inline function low_shear_switch(τ::Float64, params::SimulationParameters)
    ratio = (τ / params.clot.tau_dep)^params.clot.dep_hill
    return 1.0 / (1.0 + ratio)
end

@inline function erosion_switch(τ::Float64, params::SimulationParameters)
    ratio = (τ / params.clot.tau_ero)^params.clot.ero_hill
    return ratio / (1.0 + ratio)
end

@inline transport_porosity(φ::Float64, params::SimulationParameters) =
    clamp(1.0 - params.clot.occlusion_sensitivity * φ, params.flow.minimum_open_fraction, 1.0)

@inline function surface_site_weight(
    wall_weight::Float64,
    bound::Float64,
    φ::Float64,
    params::SimulationParameters,
)
    amplification = 1.0 + params.clot.aggregation_site_strength * (bound + φ)
    return clamp(wall_weight * amplification, 0.0, 1.0)
end

@inline function platelet_effective_diffusivity(
    τ::Float64,
    porosity::Float64,
    activated::Bool,
    params::SimulationParameters,
)
    base = activated ? params.transport.activated_diffusivity : params.transport.platelet_diffusivity
    return base * porosity^params.transport.platelet_porosity_exponent + params.transport.shear_diffusivity * τ
end

@inline agonist_effective_diffusivity(porosity::Float64, params::SimulationParameters) =
    params.transport.agonist_diffusivity * porosity^params.transport.agonist_porosity_exponent

@inline function activation_rate(
    resting::Float64,
    agonist::Float64,
    wall_weight::Float64,
    τ::Float64,
    params::SimulationParameters,
)
    shear_term = params.clot.k_activate_shear * τ / (params.clot.tau_activate + τ + eps())
    agonist_term = params.clot.k_activate_agonist * agonist / (params.clot.agonist_halfmax + agonist + eps())
    wall_term = params.clot.k_activate_wall * wall_weight
    return (wall_term + shear_term + agonist_term) * resting
end

@inline function binding_rate(
    activated::Float64,
    bound::Float64,
    φ::Float64,
    wall_weight::Float64,
    τ::Float64,
    params::SimulationParameters,
)
    occupancy = params.clot.bound_capacity > 0.0 ? clamp(bound / params.clot.bound_capacity, 0.0, 1.0) : 1.0
    return params.clot.k_dep *
           params.clot.binding_ratio *
           activated *
           surface_site_weight(wall_weight, bound, φ, params) *
           low_shear_switch(τ, params) *
           (1.0 - occupancy)
end

@inline function bound_detachment_rate(bound::Float64, τ::Float64, params::SimulationParameters)
    return params.clot.k_ero * params.clot.detachment_ratio * erosion_switch(τ, params) * bound
end

@inline function compaction_rate(
    bound::Float64,
    φ::Float64,
    wall_weight::Float64,
    params::SimulationParameters,
)
    return params.clot.k_dep *
           params.clot.compaction_ratio *
           bound *
           surface_site_weight(wall_weight, bound, φ, params) *
           (1.0 - φ)
end

@inline function agonist_source_rate(
    activated::Float64,
    bound::Float64,
    φ::Float64,
    wall_weight::Float64,
    params::SimulationParameters,
)
    site = surface_site_weight(wall_weight, bound, φ, params)
    return params.clot.k_dep * params.clot.agonist_release_ratio * site * (0.5 * activated + bound + 0.5 * φ)
end

function stable_timestep(flow::FlowState, grid::GeometryGrid, params::SimulationParameters)
    h = minimum((grid.dx, grid.dy, grid.dz))
    umax = maximum(flow.speed)
    adv_dt = umax > 0.0 ? params.numerics.cfl * h / umax : params.numerics.max_dt

    τmax = maximum(flow.shear)
    platelet_diff = platelet_effective_diffusivity(τmax, 1.0, true, params)
    diff_max = max(platelet_diff, params.transport.agonist_diffusivity)
    diff_dt = diff_max > 0.0 ? params.numerics.diffusion_cfl * h^2 / diff_max : params.numerics.max_dt

    react_rate = max(
        params.clot.k_dep * (params.clot.binding_ratio + params.clot.compaction_ratio + params.clot.agonist_release_ratio),
        params.clot.k_ero * (1.0 + params.clot.detachment_ratio),
        params.clot.k_activate_wall + params.clot.k_activate_shear + params.clot.k_activate_agonist,
        params.transport.agonist_decay,
    )
    react_dt = react_rate > 0.0 ? 0.20 / react_rate : params.numerics.max_dt
    return clamp(min(adv_dt, diff_dt, react_dt), params.numerics.min_dt, params.numerics.max_dt)
end

function update_transport!(
    next::TransportState,
    current::TransportState,
    clot::ClotState,
    flow::FlowState,
    grid::GeometryGrid,
    params::SimulationParameters,
    cache::GridIndexCache,
    dt::Float64,
)
    copyto!(next.resting, current.resting)
    copyto!(next.activated, current.activated)
    copyto!(next.agonist, current.agonist)

    resting = current.resting
    activated = current.activated
    agonist = current.agonist
    bound = clot.bound
    φ = clot.solid

    @inbounds Threads.@threads for n in eachindex(cache.interior)
        idx = cache.interior[n]
        i, j, k = Tuple(idx)
        pr = resting[idx]
        pa = activated[idx]
        ag = agonist[idx]
        solid = φ[idx]
        local_bound = bound[idx]
        porosity = transport_porosity(solid, params)

        pr_im1 = grid.mask[i - 1, j, k] ? resting[i - 1, j, k] : pr
        pr_ip1 = grid.mask[i + 1, j, k] ? resting[i + 1, j, k] : pr
        pr_jm1 = grid.mask[i, j - 1, k] ? resting[i, j - 1, k] : pr
        pr_jp1 = grid.mask[i, j + 1, k] ? resting[i, j + 1, k] : pr
        pr_km1 = grid.mask[i, j, k - 1] ? resting[i, j, k - 1] : pr
        pr_kp1 = grid.mask[i, j, k + 1] ? resting[i, j, k + 1] : pr

        pa_im1 = grid.mask[i - 1, j, k] ? activated[i - 1, j, k] : pa
        pa_ip1 = grid.mask[i + 1, j, k] ? activated[i + 1, j, k] : pa
        pa_jm1 = grid.mask[i, j - 1, k] ? activated[i, j - 1, k] : pa
        pa_jp1 = grid.mask[i, j + 1, k] ? activated[i, j + 1, k] : pa
        pa_km1 = grid.mask[i, j, k - 1] ? activated[i, j, k - 1] : pa
        pa_kp1 = grid.mask[i, j, k + 1] ? activated[i, j, k + 1] : pa

        ag_im1 = grid.mask[i - 1, j, k] ? agonist[i - 1, j, k] : ag
        ag_ip1 = grid.mask[i + 1, j, k] ? agonist[i + 1, j, k] : ag
        ag_jm1 = grid.mask[i, j - 1, k] ? agonist[i, j - 1, k] : ag
        ag_jp1 = grid.mask[i, j + 1, k] ? agonist[i, j + 1, k] : ag
        ag_km1 = grid.mask[i, j, k - 1] ? agonist[i, j, k - 1] : ag
        ag_kp1 = grid.mask[i, j, k + 1] ? agonist[i, j, k + 1] : ag

        ux = flow.velocity[i, j, k, 1] * porosity
        uy = flow.velocity[i, j, k, 2] * porosity
        uz = flow.velocity[i, j, k, 3] * porosity

        pr_dx = ux >= 0.0 ? (pr - pr_im1) / grid.dx : (pr_ip1 - pr) / grid.dx
        pr_dy = uy >= 0.0 ? (pr - pr_jm1) / grid.dy : (pr_jp1 - pr) / grid.dy
        pr_dz = uz >= 0.0 ? (pr - pr_km1) / grid.dz : (pr_kp1 - pr) / grid.dz
        pa_dx = ux >= 0.0 ? (pa - pa_im1) / grid.dx : (pa_ip1 - pa) / grid.dx
        pa_dy = uy >= 0.0 ? (pa - pa_jm1) / grid.dy : (pa_jp1 - pa) / grid.dy
        pa_dz = uz >= 0.0 ? (pa - pa_km1) / grid.dz : (pa_kp1 - pa) / grid.dz
        ag_dx = ux >= 0.0 ? (ag - ag_im1) / grid.dx : (ag_ip1 - ag) / grid.dx
        ag_dy = uy >= 0.0 ? (ag - ag_jm1) / grid.dy : (ag_jp1 - ag) / grid.dy
        ag_dz = uz >= 0.0 ? (ag - ag_km1) / grid.dz : (ag_kp1 - ag) / grid.dz

        advection_pr = ux * pr_dx + uy * pr_dy + uz * pr_dz
        advection_pa = ux * pa_dx + uy * pa_dy + uz * pa_dz
        advection_ag = ux * ag_dx + uy * ag_dy + uz * ag_dz

        lap_pr =
            (pr_ip1 - 2.0 * pr + pr_im1) / grid.dx^2 +
            (pr_jp1 - 2.0 * pr + pr_jm1) / grid.dy^2 +
            (pr_kp1 - 2.0 * pr + pr_km1) / grid.dz^2
        lap_pa =
            (pa_ip1 - 2.0 * pa + pa_im1) / grid.dx^2 +
            (pa_jp1 - 2.0 * pa + pa_jm1) / grid.dy^2 +
            (pa_kp1 - 2.0 * pa + pa_km1) / grid.dz^2
        lap_ag =
            (ag_ip1 - 2.0 * ag + ag_im1) / grid.dx^2 +
            (ag_jp1 - 2.0 * ag + ag_jm1) / grid.dy^2 +
            (ag_kp1 - 2.0 * ag + ag_km1) / grid.dz^2

        τ = flow.shear[idx]
        w = cache.wall_weight[idx]
        dpr = platelet_effective_diffusivity(τ, porosity, false, params)
        dpa = platelet_effective_diffusivity(τ, porosity, true, params)
        dag = agonist_effective_diffusivity(porosity, params)
        activation = activation_rate(pr, ag, w, τ, params)
        bind = binding_rate(pa, local_bound, solid, w, τ, params)
        detach = bound_detachment_rate(local_bound, τ, params)
        agonist_source = agonist_source_rate(pa, local_bound, solid, w, params)

        next.resting[idx] = clamp(
            pr + dt * (dpr * lap_pr - advection_pr - activation),
            0.0,
            params.transport.inlet_concentration,
        )
        next.activated[idx] = clamp(
            pa + dt * (dpa * lap_pa - advection_pa + activation - bind + detach),
            0.0,
            params.transport.inlet_concentration,
        )
        next.agonist[idx] = clamp(
            ag + dt * (dag * lap_ag - advection_ag + agonist_source - params.transport.agonist_decay * ag),
            0.0,
            1.0,
        )
    end

    inlet_activated = clamp(params.transport.inlet_activated_fraction, 0.0, 1.0) * params.transport.inlet_concentration
    inlet_resting = params.transport.inlet_concentration - inlet_activated
    @inbounds for idx in cache.inlet
        next.resting[idx] = inlet_resting
        next.activated[idx] = inlet_activated
        next.agonist[idx] = params.transport.inlet_agonist
    end

    return next
end

function total_mobile_platelets(state::TransportState)
    total = similar(state.resting)
    @inbounds for idx in eachindex(total)
        total[idx] = state.resting[idx] + state.activated[idx]
    end
    return total
end
