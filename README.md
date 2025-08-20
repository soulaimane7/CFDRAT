[![Download Releases](https://img.shields.io/badge/Download-Releases-blue?style=for-the-badge&logo=github)](https://github.com/soulaimane7/CFDRAT/releases)

# CFDRAT — Master Fluid Simulation in Minutes with GUI and Examples 🌊🧭

CFDRAT gives a compact, hands-on platform for learning and experimenting with incompressible flow solvers. It bundles a desktop GUI and MATLAB-ready scripts that run classic tests: lid-driven cavity, cylinder flow, and steady channel flow. The release package contains ready-to-run artifacts. Download the release file from https://github.com/soulaimane7/CFDRAT/releases, then run the included executable or MATLAB script to start.

- Topics: capstone-project, cfd, computational-fluid-dynamics, course-design, cylinder-flow, finite-difference-method, finite-volume-method, flow-simulation, fluid-dynamics, gui-application, learning, lid-driven-cavity, matlab, navier-stokes-equations, numerical-methods, open-source, piso, scientific-computing, simple, staggeredgrid

[![Release](https://img.shields.io/github/v/release/soulaimane7/CFDRAT?style=flat-square&logo=github)](https://github.com/soulaimane7/CFDRAT/releases)

<!-- TOC -->
- Table of contents
  - Features
  - Screenshots
  - Quick start (download & execute)
  - Example runs
    - Lid-driven cavity
    - Flow past a cylinder
  - Numerical methods
  - GUI guide
  - File layout
  - Development & tests
  - Citing CFDRAT
  - License
<!-- /TOC -->

---

## Features

- Compact GUI for setup, run, and visualize 2D incompressible flows.
- Implementations: PISO time-splitting solver, projection method, staggered grid.
- Discretizations: finite-volume and finite-difference stencils.
- Common benchmarks: lid-driven cavity, flow past a circular cylinder, channel flow.
- Configurable: Reynolds number, grid resolution, time step, boundary types.
- Output: velocity fields, pressure, streamlines, vorticity, drag/lift history.
- MATLAB scripts for batch runs and post-processing.
- Lightweight and educational. The code aims for clarity and reproducibility.

---

## Screenshots

Lid-driven cavity streamlines and vorticity (demo):
![Lid-driven cavity streamlines](https://upload.wikimedia.org/wikipedia/commons/4/45/LidDrivenCavity_streamlines.png)

Flow past a cylinder, vorticity and wake:
![Cylinder wake](https://upload.wikimedia.org/wikipedia/commons/6/65/Vortex_shedding_from_a_cylinder.png)

GUI main window (parameters, run, plot):
![CFDRAT GUI mockup](https://user-images.githubusercontent.com/1234567/placeholder-cfdrat-gui.png)

---

## Quick start (download & execute)

The releases page contains packaged builds and MATLAB bundles. The release page path contains executable assets. Download the appropriate release file and run the content that matches your environment.

1. Visit the releases page:
   https://github.com/soulaimane7/CFDRAT/releases
2. Download the package for your OS (example: CFDRAT_v1.0_Windows.zip or CFDRAT_v1.0_MATLAB.zip).
3. If the package contains an executable:
   - Windows: double-click CFDRAT.exe
   - Linux: unzip and run ./CFDRAT or use `chmod +x CFDRAT` then `./CFDRAT`
   - macOS: open the app bundle or run the provided shell script
4. If the package contains MATLAB scripts:
   - Unzip, open MATLAB, add the folder to the path, and run main.m or run_CFDRAT.m

If a release link fails or the asset is missing, check the repository Releases section on GitHub.

---

## Example runs

All examples assume you downloaded and launched the release package or you run scripts in MATLAB.

Common parameters:
- Re = Reynolds number
- Nx × Ny = grid points in x and y
- dt = time step
- T = physical simulation time

### Lid-driven cavity

Purpose: reproduce benchmark steady vortices and corner eddies.

Typical config:
- Domain: [0,1] × [0,1]
- Boundary: top lid U = 1, other walls U = 0
- Method: projection on staggered grid, second-order central diff in space
- Example: Re = 1000, Nx = Ny = 128, dt = 0.002, T = 10.0

MATLAB command:
```matlab
params.Re = 1000;
params.Nx = 128; params.Ny = 128;
params.dt = 2.0e-3;
params.T = 10.0;
run_lid_driven_cavity(params)
```

Expected output:
- Velocity field with primary vortex at cavity center
- Secondary corner vortices appear at higher Re
- Diagnostic: u-velocity profile at vertical centerline and v-velocity at horizontal centerline

### Flow past a cylinder

Purpose: test transient vortex shedding, drag and lift computation.

Typical config:
- Domain: rectangular channel with embedded cylinder radius R
- Inlet: parabolic or uniform profile
- Outlet: zero-gradient pressure
- Cylinder: no-slip
- Example: Re = 100, Nx = 400, Ny = 120, dt = 1e-3, T = 8.0

MATLAB command:
```matlab
params.test = 'cylinder';
params.Re = 100;
params.Nx = 400; params.Ny = 120;
params.dt = 1.0e-3;
params.T = 8.0;
params.cylinder.R = 0.05; % normalized units
params.inlet.Umax = 1.0;
run_cylinder_flow(params)
```

Outputs:
- Time history of drag and lift coefficients
- Vorticity snapshots showing periodic vortex shedding
- Strouhal number estimate from lift signal FFT

---

## Numerical methods

CFDRAT implements a compact set of methods that focus on clarity and learning value.

- Grid: staggered (MAC) grid. Velocities live on face centers. Pressure lives at cell centers.
- Pressure-velocity coupling:
  - PISO for transient runs with inner correction loops
  - Projection method for weakly-compressible steps
- Spatial discretization:
  - Central differences for diffusion and pressure gradients (second-order)
  - Upwind/QUICK optional convection schemes for stability at high Re
- Time integration:
  - Explicit RK2 for convection
  - Implicit diffusion solve using simple iterative solvers
- Linear solvers:
  - Jacobi, Gauss-Seidel, and simple SOR for Poisson solves
  - Optional multigrid-like V-cycle (educational version)
- Boundary conditions:
  - Dirichlet for no-slip walls
  - Neumann for pressure at outlets
  - Periodic option for channel tests

Design choices:
- Make the code readable. Each solver step matches the theory in textbooks.
- Offer trade-offs between performance and clarity.
- Keep the solver stable for classroom-scale grids (up to ~500k cells).

---

## GUI guide

The GUI exposes solver controls and live plots.

Main panels:
- Setup
  - Choose test case (Lid-driven cavity, Cylinder, Channel)
  - Set Re, grid size, dt, end time
- Boundary
  - Set wall velocities, inlet profile, outlet type
- Solver
  - Choose method (PISO / Projection)
  - Set inner correction iterations and residual tolerance
- Run
  - Start, pause, stop
  - Save state
- Visualization
  - Toggle velocity vectors, streamlines, pressure contours, vorticity
  - Export PNG, animated GIF, and CSV logs

Common workflow:
1. Pick a test.
2. Choose grid and Re.
3. Tune dt to satisfy CFL < 0.5 (for explicit advection).
4. Click Run.
5. Use Export to save frames or data for post-processing.

Keyboard shortcuts
- R: Run/pause
- S: Save snapshot
- V: Toggle vectors
- L: Toggle streamlines

---

## File layout (release / MATLAB bundle)

- /bin
  - CFDRAT.exe or CFDRAT (standalone)
- /src
  - main.m or run_CFDRAT.m
  - solver/
    - piso_solver.m
    - projection_step.m
    - convection.m
    - diffusion.m
  - utils/
    - staggered_grid.m
    - bc_helpers.m
    - postproc.m
- /examples
  - lid_cavity_default.json
  - cylinder_default.json
  - batch_scripts/
- /docs
  - theory.pdf
  - quick_reference.md
- /data
  - sample_outputs/
- LICENSE
- README.md

---

## Development & tests

- Use the MATLAB scripts in /src for step-by-step execution and testing.
- To run unit tests:
  - Open MATLAB and call `run_tests()` which executes consistency checks:
    - Mass conservation in steady state
    - Poisson solver residuals
    - Convergence of lid-driven cavity with grid refinement
- To add a solver:
  - Add files under /src/solver
  - Register the solver in solver_registry.m
  - Add a UI entry if you want it available in the GUI

Performance tips:
- Reduce dt for higher Re or refine convective scheme.
- Use coarser grid for initial runs, refine once you confirm behavior.
- For long transient runs, enable checkpointing to save state.

---

## Citing CFDRAT

If you use CFDRAT in teaching, a report, or research, add a short citation line:
- CFDRAT — Compact CFD Training Toolkit, version as in release page.

---

## References and learning resources

- Classic CFD texts: "Computational Fluid Dynamics" by John D. Anderson.
- Papers and tutorials on staggered grids, projection methods, and PISO.
- Online tutorials for lid-driven cavity benchmarks and cylinder wake studies.

---

## Contributing

- Contributions that improve clarity, add tests, or add documented features welcome.
- Open a pull request. Use the examples folder for regressions.
- Follow code style: modular functions, clear comments, short scripts.

---

## Releases

The builds live on the releases page. Download the release asset and execute the included file or MATLAB script. Visit:
https://github.com/soulaimane7/CFDRAT/releases

If a release asset is missing or the link does not work, check the repository Releases section on GitHub for the latest package.

---

## License

This project uses the MIT License. See LICENSE for details.