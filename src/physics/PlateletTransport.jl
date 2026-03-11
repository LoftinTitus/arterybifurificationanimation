initialize_platelets(grid::GeometryGrid, params::SimulationParameters) =
    params.transport.inlet_concentration .* Float64.(grid.mask)

function inlet_region(grid::GeometryGrid, params::SimulationParameters)
    cutoff = params.transport.inlet_buffer_fraction * params.geometry.inlet_length
    return (grid.segment_id .== 1) .& (grid.axial_coordinate .<= cutoff)
end

function low_shear_switch(τ::Float64, params::SimulationParameters)
    ratio = (τ / params.clot.tau_dep)^params.clot.dep_hill
    return 1.0 / (1.0 + ratio)
end

function platelet_sink(C::Float64, φ::Float64, wall_distance::Float64, τ::Float64, params::SimulationParameters)
    wall_weight = exp(-(wall_distance / params.clot.wall_bandwidth)^2)
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

function update_platelets(
    C::Array{Float64, 3},
    φ::Array{Float64, 3},
    flow::FlowState,
    grid::GeometryGrid,
    params::SimulationParameters,
    dt::Float64,
)
    nx, ny, nz = size(C)
    Cnew = copy(C)

    @inbounds for k in 2:(nz - 1), j in 2:(ny - 1), i in 2:(nx - 1)
        if !grid.mask[i, j, k]
            Cnew[i, j, k] = 0.0
            continue
        end

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

        sink = platelet_sink(c, φ[i, j, k], grid.wall_distance[i, j, k], flow.shear[i, j, k], params)
        Cnew[i, j, k] = clamp(
            c + dt * (params.transport.platelet_diffusivity * laplacian - advection - sink),
            0.0,
            params.transport.inlet_concentration,
        )
    end

    Cnew[inlet_region(grid, params)] .= params.transport.inlet_concentration
    Cnew[.!grid.mask] .= 0.0
    return Cnew
end
