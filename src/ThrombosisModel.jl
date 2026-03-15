module ThrombosisModel

using DelimitedFiles
using GeometryBasics
using GLMakie
using LinearAlgebra
using StaticArrays
using Statistics

export GeometryParameters,
       FlowParameters,
       TransportParameters,
       ClotParameters,
       NumericsParameters,
       VisualizationParameters,
       SimulationParameters,
       GeometryGrid,
       FlowState,
       TransportState,
       ClotState,
       SimulationOutputs,
       default_parameters,
       build_bifurcation_grid,
       surface_meshes,
       compute_flow_field,
       initialize_platelets,
       initialize_transport_state,
       initialize_clot,
       initialize_clot_state,
       stable_timestep,
       run_simulation,
       write_metrics_csv,
       plot_geometry,
       plot_simulation_summary,
       plot_fields,
       plot_slice,
       animate_clot_growth

include("Parameters.jl")
include("geometry/BifurcationGeometry.jl")
include("physics/FlowField.jl")
include("physics/PlateletTransport.jl")
include("physics/ClotGrowth.jl")
include("solvers/Simulation.jl")
include("visualization/MakiePlots.jl")

end
