initialize_platelets(grid::GeometryGrid, params::SimulationParameters) =
    params.transport.inlet_concentration .* Float64.(grid.mask)

function low_shear_switch(τ::Float64, params::SimulationParameters)
    ratio = (τ / params.clot.tau_dep)^params.clot.dep_hill
    return 1.0 / (1.0 + ratio)
end

function platelet_sink(C::Float64, φ::Float64, wall_weight::Float64, τ::Float64, params::SimulationParameters)
    return params.transport.deposition_sink * C * wall_weight * low_shear_switch(τ, params) * (1.0 - φ)
end

function stable_timestep(flow::FlowState, grid::GeometryGrid, params::SimulationParameters)
    h = minimum((grid.dx, grid.dy, grid.dz))
    umax = maximum(flow.speed)
    adv_dt = umax > 0.0 ? params.numerics.cfl * h / umax : params.numerics.max_dt
    diff_dt = params.transport.platelet_diffusivity > 0.0 ?
        params.numerics.diffusion_cfl * h^2 / params.transport.platelet_diffusivity :
        params.numerics.max_dt
    react_rate = max(
        params.transport.deposition_sink + params.clot.k_dep * params.transport.inlet_concentration,
        params.clot.k_ero,
    )
    react_dt = react_rate > 0.0 ? 0.25 / react_rate : params.numerics.max_dt
    return clamp(min(adv_dt, diff_dt, react_dt), params.numerics.min_dt, params.numerics.max_dt)
end

function update_platelets!(
    Cnew::Array{Float64, 3},
    C::Array{Float64, 3},
    φ::Array{Float64, 3},
    flow::FlowState,
    grid::GeometryGrid,
    params::SimulationParameters,
    cache::GridIndexCache,
    dt::Float64,
)
    copyto!(Cnew, C)

    @inbounds for idx in cache.interior
        i, j, k = Tuple(idx)
        c = C[i, j, k]
        c_im1 = grid.mask[i - 1, j, k] ? C[i - 1, j, k] : c
        c_ip1 = grid.mask[i + 1, j, k] ? C[i + 1, j, k] : c
        c_jm1 = grid.mask[i, j - 1, k] ? C[i, j - 1, k] : c
        c_jp1 = grid.mask[i, j + 1, k] ? C[i, j + 1, k] : c
        c_km1 = grid.mask[i, j, k - 1] ? C[i, j, k - 1] : c
        c_kp1 = grid.mask[i, j, k + 1] ? C[i, j, k + 1] : c

        ux = flow.velocity[i, j, k, 1]
        uy = flow.velocity[i, j, k, 2]
        uz = flow.velocity[i, j, k, 3]

        dc_dx = ux >= 0.0 ? (c - c_im1) / grid.dx : (c_ip1 - c) / grid.dx
        dc_dy = uy >= 0.0 ? (c - c_jm1) / grid.dy : (c_jp1 - c) / grid.dy
        dc_dz = uz >= 0.0 ? (c - c_km1) / grid.dz : (c_kp1 - c) / grid.dz
        advection = ux * dc_dx + uy * dc_dy + uz * dc_dz

        laplacian =
            (c_ip1 - 2.0 * c + c_im1) / grid.dx^2 +
            (c_jp1 - 2.0 * c + c_jm1) / grid.dy^2 +
            (c_kp1 - 2.0 * c + c_km1) / grid.dz^2

        sink = platelet_sink(c, φ[idx], cache.wall_weight[idx], flow.shear[idx], params)
        Cnew[idx] = clamp(
            c + dt * (params.transport.platelet_diffusivity * laplacian - advection - sink),
            0.0,
            params.transport.inlet_concentration,
        )
    end

    @inbounds for idx in cache.inlet
        Cnew[idx] = params.transport.inlet_concentration
    end
    return Cnew
end
