using Pkg
Pkg.activate(@__DIR__)

using GLMakie
using ThrombosisModel

params = default_parameters()
results = run_simulation(params)

mkpath(joinpath(@__DIR__, "output"))
write_metrics_csv(results, joinpath(@__DIR__, "output", "metrics.csv"))
GLMakie.save(joinpath(@__DIR__, "output", "summary.png"), plot_simulation_summary(results))
GLMakie.save(joinpath(@__DIR__, "output", "geometry.png"), plot_geometry(params))
GLMakie.save(joinpath(@__DIR__, "output", "fields.png"), plot_fields(results))
GLMakie.save(joinpath(@__DIR__, "output", "slice.png"), plot_slice(results, field = :shear))
animate_clot_growth(results, joinpath(@__DIR__, "output", "clot_growth.mp4"))
