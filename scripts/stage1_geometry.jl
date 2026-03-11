using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using GLMakie
using ThrombosisModel

params = default_parameters()
fig = plot_geometry(params)
mkpath(joinpath(@__DIR__, "..", "output"))
GLMakie.save(joinpath(@__DIR__, "..", "output", "stage1_geometry.png"), fig)
