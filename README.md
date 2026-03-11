# ThrombosisModel.jl

This project implements a modular Julia foundation for modeling thrombus formation in a 3D arterial bifurcation. The current implementation is intentionally simplified so it can run on a personal workstation while still preserving the main physical couplings:

- parameterized Y-bifurcation geometry
- laminar Poiseuille-like flow in each branch
- 3D platelet convection-diffusion on a Cartesian grid
- wall-localized thrombus growth with shear-dependent deposition and erosion
- lumen narrowing through a clot-dependent open-fraction model
- Makie-based 3D plots and animation

## Project structure

```text
arterybifurificationanimation/
├── Project.toml
├── README.md
├── parameters.jl
├── simulation.jl
├── scripts/
│   ├── stage1_geometry.jl
│   ├── stage2_transport.jl
│   ├── stage3_clot_growth.jl
│   └── stage4_coupled.jl
├── src/
│   ├── ThrombosisModel.jl
│   ├── Parameters.jl
│   ├── geometry/
│   │   └── BifurcationGeometry.jl
│   ├── physics/
│   │   ├── FlowField.jl
│   │   ├── PlateletTransport.jl
│   │   └── ClotGrowth.jl
│   ├── solvers/
│   │   └── Simulation.jl
│   └── visualization/
│       └── MakiePlots.jl
└── output/
```

## Governing model

### 1. Geometry

The vessel lumen is represented as the union of three capsules:

- one inlet segment
- two daughter branches

Each segment is parameterized by length, radius, and axis direction. This gives a robust implicit geometry for grid-based transport and a separate tube surface mesh for visualization.

### 2. Flow model

The current solver does not solve full 3D Navier-Stokes. Instead, it assigns a steady laminar velocity field aligned with the local centerline:

```math
u(r) = 2 U_{mean}\left(1 - \left(\frac{r}{R_{eff}}\right)^2\right)
```

where:

- `r` is distance from the local centerline
- `R_eff = R sqrt(chi)` is the clot-reduced effective radius
- `chi = 1 - alpha phi` is the local open fraction

Branch flow split is based on a simple Poiseuille conductance law:

```math
Q_i \propto \frac{R_i^4}{L_i}
```

Wall shear stress is approximated as:

```math
tau_w = \frac{4 mu U_{mean}}{R_{eff}}
```

This field is reused in the transport and clot kinetics.

### 3. Platelet transport

Platelets satisfy:

```math
\frac{\partial C}{\partial t} + u \cdot \nabla C = D \nabla^2 C - R_p
```

with a deposition sink:

```math
R_p = k_p C w(d_w) f(tau) (1 - phi)
```

where `w(d_w)` localizes deposition near the wall and `f(tau)` suppresses deposition at high shear.

The spatial discretization is:

- first-order upwind advection
- second-order central diffusion
- explicit time stepping with CFL and reaction limits

### 4. Thrombus growth

The thrombus volume fraction `phi` evolves as:

```math
\frac{\partial phi}{\partial t}
= k_{dep} C w(d_w) f(tau) (1 - phi)
- k_{ero} g(tau) phi
```

with:

```math
f(tau) = \frac{1}{1 + (tau/tau_{dep})^{n_d}}
```

```math
g(tau) = \frac{(tau/tau_{ero})^{n_e}}{1 + (tau/tau_{ero})^{n_e}}
```

Interpretation:

- low shear favors platelet deposition and clot buildup
- high shear favors erosion
- growth is strongest within a wall-normal band of thickness `wall_bandwidth`

### 5. Lumen narrowing

Clot feedback enters the flow model through the open fraction:

```math
chi = clamp(1 - alpha phi, chi_{min}, 1)
```

This reduces the local effective lumen radius and therefore changes velocity and wall shear stress as the clot grows.

## Running the stages

Instantiate the environment first:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Then run:

```bash
julia --project=. scripts/stage1_geometry.jl
julia --project=. scripts/stage2_transport.jl
julia --project=. scripts/stage3_clot_growth.jl
julia --project=. scripts/stage4_coupled.jl
```

Or run the full pipeline:

```bash
julia --project=. simulation.jl
```

Generated outputs are written to `output/`.

## Modeling assumptions

This is a reduced-order model, not full CFD. The main assumptions are:

1. Blood is treated as a Newtonian fluid with steady laminar branch-wise flow.
2. The lumen is represented on a Cartesian grid using an implicit union-of-capsules geometry.
3. Platelet transport uses scalar concentration rather than discrete particle or bond kinetics.
4. Clot growth is represented as a continuum volume fraction field.
5. Wall shear stress is estimated from the local Poiseuille profile rather than computed from a full finite-element velocity gradient.
6. Lumen narrowing is local and algebraic, not a remeshed moving-boundary solve.

These assumptions keep the model computationally cheap while retaining meaningful couplings between transport, shear, and occlusion.

## Numerical stability guidance

Use these constraints when modifying parameters:

1. Keep the explicit advection CFL below about `0.3-0.4`.
2. Keep the diffusion limit near `D dt / h^2 < 0.2`.
3. If `k_dep` or `k_ero` are increased significantly, reduce `max_dt`.
4. If the geometry resolution is refined, expect the timestep to drop.
5. Avoid very low `minimum_open_fraction` values; they can create unrealistically large velocities.
6. If oscillations appear in `C`, switch to flux-limited advection or semi-Lagrangian transport.

## Extending to research-grade CFD

The clean extension path is:

1. Replace `compute_flow_field` with a finite-element Navier-Stokes solve using Gridap.jl or Ferrite.jl on an actual bifurcation mesh.
2. Compute wall shear stress directly from the velocity gradient at the wall.
3. Move from an implicit Cartesian grid to an unstructured mesh with boundary markers.
4. Treat thrombus as a porous Brinkman region:

   ```math
   -\nabla p + mu \nabla^2 u - beta(phi) u = 0
   ```

5. Add biochemical species such as thrombin, fibrinogen, and ADP using reaction-transport equations.
6. Introduce pulsatile inflow and outlet impedance models for more realistic hemodynamics.
7. Add adaptive remeshing or level-set geometry updates if the clot front needs sharper resolution.

## Suggested next upgrades

- fit deposition and erosion constants to literature data
- add asymmetric outlet diameters and flow split studies
- couple platelets to residence time metrics
- export VTK files for ParaView
- add parameter sweeps and sensitivity analysis

## Notes on package choices

Current code uses:

- `GeometryBasics.jl` for vessel surface meshes
- `GLMakie.jl` for 3D visualization and animation
- `StaticArrays.jl` for compact geometry calculations

For full CFD, add:

- `Gridap.jl` or `Ferrite.jl` for finite-element flow
- `DifferentialEquations.jl` for time integration and parameter estimation workflows
- `ModelingToolkit.jl` for symbolic model variants and reduced-order reformulations
