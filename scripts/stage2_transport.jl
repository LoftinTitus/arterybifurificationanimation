using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using GLMakie
using ThrombosisModel

params = default_parameters()
params = SimulationParameters(
    geometry = params.geometry,
    flow = params.flow,
    transport = params.transport,
    clot = ClotParameters(k_dep = 0.0, k_ero = 0.0, wall_bandwidth = params.clot.wall_bandwidth),
    numerics = NumericsParameters(tfinal = 6.0, save_every = 4),
    visualization = params.visualization,
)

results = run_simulation(params)
mkpath(joinpath(@__DIR__, "..", "output"))
GLMakie.save(joinpath(@__DIR__, "..", "output", "stage2_transport_slice.png"), plot_slice(results, field = :platelets))
