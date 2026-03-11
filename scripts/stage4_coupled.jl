using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using GLMakie
using ThrombosisModel

params = default_parameters()
results = run_simulation(params)

mkpath(joinpath(@__DIR__, "..", "output"))
write_metrics_csv(results, joinpath(@__DIR__, "..", "output", "stage4_metrics.csv"))
GLMakie.save(joinpath(@__DIR__, "..", "output", "stage4_summary.png"), plot_simulation_summary(results))
animate_clot_growth(results, joinpath(@__DIR__, "..", "output", "stage4_clot_growth.mp4"))
