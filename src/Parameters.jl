Base.@kwdef struct GeometryParameters
    inlet_diameter::Float64 = 4.0e-3
    branch_diameters::NTuple{2, Float64} = (3.2e-3, 3.2e-3)
    bifurcation_angle::Float64 = deg2rad(55.0)
    inlet_length::Float64 = 18.0e-3
    branch_length::Float64 = 18.0e-3
    padding::Float64 = 3.0e-3
end

Base.@kwdef struct FlowParameters
    density::Float64 = 1060.0
    viscosity::Float64 = 3.5e-3
    inlet_mean_velocity::Float64 = 0.18
    outlet_pressure::Float64 = 0.0
    max_velocity_multiplier::Float64 = 2.5
    minimum_open_fraction::Float64 = 0.08
end

Base.@kwdef struct TransportParameters
    platelet_diffusivity::Float64 = 2.5e-10
    activated_diffusivity::Float64 = 2.8e-10
    agonist_diffusivity::Float64 = 1.1e-9
    shear_diffusivity::Float64 = 4.0e-11
    inlet_concentration::Float64 = 1.0
    inlet_activated_fraction::Float64 = 0.0
    inlet_agonist::Float64 = 0.0
    deposition_sink::Float64 = 0.18
    platelet_porosity_exponent::Float64 = 1.35
    agonist_porosity_exponent::Float64 = 0.85
    agonist_decay::Float64 = 0.45
    inlet_buffer_fraction::Float64 = 0.08
end

Base.@kwdef struct ClotParameters
    k_dep::Float64 = 0.16
    k_ero::Float64 = 0.05
    tau_dep::Float64 = 1.2
    tau_ero::Float64 = 4.0
    dep_hill::Float64 = 3.0
    ero_hill::Float64 = 2.0
    wall_bandwidth::Float64 = 4.5e-4
    occlusion_sensitivity::Float64 = 0.95
    solid_threshold::Float64 = 0.92
    k_activate_wall::Float64 = 0.22
    k_activate_shear::Float64 = 0.10
    tau_activate::Float64 = 2.0
    k_activate_agonist::Float64 = 0.35
    agonist_halfmax::Float64 = 0.12
    binding_ratio::Float64 = 1.0
    compaction_ratio::Float64 = 0.80
    detachment_ratio::Float64 = 0.50
    agonist_release_ratio::Float64 = 0.45
    aggregation_site_strength::Float64 = 0.65
    bound_capacity::Float64 = 1.0
    junction_length_scale::Float64 = 2.5e-3
end

Base.@kwdef struct NumericsParameters
    nx::Int = 44
    ny::Int = 44
    nz::Int = 96
    tfinal::Float64 = 18.0
    cfl::Float64 = 0.35
    diffusion_cfl::Float64 = 0.20
    min_dt::Float64 = 1.0e-3
    max_dt::Float64 = 5.0e-2
    save_every::Int = 8
end

Base.@kwdef struct VisualizationParameters
    clot_threshold::Float64 = 0.12
    wall_alpha::Float64 = 0.18
    slice_index::Float64 = 0.5
end

Base.@kwdef struct SimulationParameters
    geometry::GeometryParameters = GeometryParameters()
    flow::FlowParameters = FlowParameters()
    transport::TransportParameters = TransportParameters()
    clot::ClotParameters = ClotParameters()
    numerics::NumericsParameters = NumericsParameters()
    visualization::VisualizationParameters = VisualizationParameters()
end

default_parameters() = SimulationParameters()
