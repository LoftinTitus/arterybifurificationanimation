using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using GLMakie
using ThrombosisModel

params = default_parameters()
params = SimulationParameters(
    geometry = params.geometry,
    flow = params.flow,
    transport = params.transport,
    clot = ClotParameters(occlusion_sensitivity = 0.0),
    numerics = NumericsParameters(tfinal = 12.0, save_every = 6),
    visualization = params.visualization,
)

results = run_simulation(params)
mkpath(joinpath(@__DIR__, "..", "output"))
GLMakie.save(joinpath(@__DIR__, "..", "output", "stage3_clot_growth.png"), plot_fields(results))
