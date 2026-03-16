using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using GLMakie
using ThrombosisModel

params = default_parameters()
println("Starting full simulation...")
results = run_simulation(params; show_progress = true)

output_dir = joinpath(@__DIR__, "..", "output")
mkpath(output_dir)

println("Writing metrics.csv")
write_metrics_csv(results, joinpath(output_dir, "metrics.csv"))
println("Writing summary.png")
GLMakie.save(joinpath(output_dir, "summary.png"), plot_simulation_summary(results))
println("Writing slice_platelets.png")
GLMakie.save(joinpath(output_dir, "slice_platelets.png"), plot_slice(results, field = :platelets))
println("Writing slice_activated.png")
GLMakie.save(joinpath(output_dir, "slice_activated.png"), plot_slice(results, field = :activated))
println("Writing slice_agonist.png")
GLMakie.save(joinpath(output_dir, "slice_agonist.png"), plot_slice(results, field = :agonist))
println("Writing slice_bound.png")
GLMakie.save(joinpath(output_dir, "slice_bound.png"), plot_slice(results, field = :bound))
println("Writing slice_clot.png")
GLMakie.save(joinpath(output_dir, "slice_clot.png"), plot_slice(results, field = :clot))
println("Writing slice_shear.png")
GLMakie.save(joinpath(output_dir, "slice_shear.png"), plot_slice(results, field = :shear))
println("Writing clot_growth.mp4")
animate_clot_growth(results, joinpath(output_dir, "clot_growth.mp4"))

println("Simulation complete.")
println("Outputs written to: $(output_dir)")
