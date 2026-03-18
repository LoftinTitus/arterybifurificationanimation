const VOXEL_MARKER = Rect3f(Point3f(-0.5, -0.5, -0.5), Vec3f(1.0, 1.0, 1.0))

const TRANSPARENT_COLOR = RGBAf(0.82, 0.08, 0.10, 0.0)

function clot_voxel_cache(grid::GeometryGrid)
    active = findall(grid.mask)
    positions = Vector{Point3f}(undef, length(active))
    @inbounds for n in eachindex(active)
        i, j, k = Tuple(active[n])
        positions[n] = Point3f(grid.x[i], grid.y[j], grid.z[k])
    end
    return active, positions
end

function update_clot_colors!(
    colors::Vector{RGBAf},
    active::Vector{CartesianIndex{3}},
    φ::AbstractArray{<:Real, 3},
    threshold::Float64,
)
    @inbounds Threads.@threads for n in eachindex(active)
        value = Float32(φ[active[n]])
        colors[n] = value >= threshold ?
            RGBAf(0.82, 0.08, 0.10, clamp(0.25f0 + 0.75f0 * value, 0.25f0, 0.95f0)) :
            TRANSPARENT_COLOR
    end
    return colors
end

function plot_geometry(params::SimulationParameters = default_parameters())
    fig = Figure(size = (1100, 700))
    ax = Axis3(fig[1, 1], title = "Parameterized Arterial Bifurcation", aspect = :data)
    for mesh_data in surface_meshes(params.geometry)
        mesh!(ax, mesh_data, color = RGBAf(0.62, 0.78, 0.95, params.visualization.wall_alpha))
    end
    ax.azimuth[] = 0.45
    ax.elevation[] = 0.35
    return fig
end

function plot_simulation_summary(results::SimulationOutputs)
    fig = Figure(size = (1100, 700))
    ax1 = Axis(fig[1, 1], title = "Clot Volume", xlabel = "Time (s)", ylabel = "Volume (m^3)")
    ax2 = Axis(fig[1, 2], title = "Occlusion", xlabel = "Time (s)", ylabel = "Percent (%)")
    ax3 = Axis(fig[2, 1:2], title = "Mean Wall Shear Stress", xlabel = "Time (s)", ylabel = "Pa")
    lines!(ax1, results.times, results.clot_volume, color = :firebrick, linewidth = 3)
    lines!(ax2, results.times, results.occlusion, color = :darkorange, linewidth = 3)
    lines!(ax3, results.times, results.mean_wall_shear, color = :steelblue, linewidth = 3)
    return fig
end

function plot_fields(results::SimulationOutputs; snapshot_index::Int = length(results.clot_snapshots))
    idx = clamp(snapshot_index, 1, length(results.clot_snapshots))
    φ = results.clot_snapshots[idx]
    fig = Figure(size = (1200, 800))
    ax = Axis3(fig[1, 1], title = "Clot Growth in 3D", aspect = :data)
    for mesh_data in surface_meshes(results.params.geometry)
        mesh!(ax, mesh_data, color = RGBAf(0.65, 0.82, 0.98, results.params.visualization.wall_alpha))
    end
    active, positions = clot_voxel_cache(results.grid)
    colors = fill(TRANSPARENT_COLOR, length(active))
    update_clot_colors!(colors, active, φ, results.params.visualization.clot_threshold)
    meshscatter!(
        ax,
        positions;
        marker = VOXEL_MARKER,
        markersize = Vec3f(results.grid.dx, results.grid.dy, results.grid.dz),
        color = colors,
    )
    ax.azimuth[] = 0.45
    ax.elevation[] = 0.35
    Label(fig[0, 1], "Time = $(round(results.snapshot_times[idx], digits = 2)) s")
    return fig
end

function plot_slice(results::SimulationOutputs; field::Symbol = :platelets, snapshot_index::Int = length(results.clot_snapshots))
    idx = clamp(snapshot_index, 1, length(results.clot_snapshots))
    y_index = clamp(round(Int, results.params.visualization.slice_index * length(results.grid.y)), 1, length(results.grid.y))
    raw_field = if field == :platelets
        results.final_platelets
    elseif field == :activated
        results.final_activated
    elseif field == :agonist
        results.final_agonist
    elseif field == :bound
        results.final_bound
    elseif field == :shear
        results.final_shear
    elseif field == :clot
        results.clot_snapshots[idx]
    else
        results.final_bound
    end
    plane = Array(raw_field[:, y_index, :])
    mask = Array(results.grid.mask[:, y_index, :])
    plane[.!mask] .= NaN32
    fig = Figure(size = (1000, 500))
    ax = Axis(fig[1, 1], title = "Mid-plane Slice: $(field)", xlabel = "x (m)", ylabel = "z (m)")
    hm = heatmap!(
        ax,
        results.grid.x,
        results.grid.z,
        plane';
        colormap = field == :shear ? :plasma : (field == :clot || field == :bound ? :reds : :viridis),
        interpolate = false,
    )
    Colorbar(fig[1, 2], hm)
    return fig
end

function animate_clot_growth(results::SimulationOutputs, filepath::AbstractString)
    fig = Figure(size = (1200, 800))
    ax = Axis3(fig[1, 1], title = "Thrombus Formation Animation", aspect = :data)
    for mesh_data in surface_meshes(results.params.geometry)
        mesh!(ax, mesh_data, color = RGBAf(0.65, 0.82, 0.98, results.params.visualization.wall_alpha))
    end

    active, positions = clot_voxel_cache(results.grid)
    colors = fill(TRANSPARENT_COLOR, length(active))
    update_clot_colors!(colors, active, results.clot_snapshots[1], results.params.visualization.clot_threshold)
    color_obs = Observable(colors)
    meshscatter!(
        ax,
        positions;
        marker = VOXEL_MARKER,
        markersize = Vec3f(results.grid.dx, results.grid.dy, results.grid.dz),
        color = color_obs,
    )

    record(fig, filepath, eachindex(results.clot_snapshots)) do frame
        update_clot_colors!(
            colors,
            active,
            results.clot_snapshots[frame],
            results.params.visualization.clot_threshold,
        )
        notify(color_obs)
        ax.azimuth[] = 0.30 + 0.85 * (frame - 1) / max(length(results.clot_snapshots) - 1, 1)
        ax.elevation[] = 0.30
    end
    return filepath
end
