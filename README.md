# 1D Heat Rod Simulator

This project is an interactive simulation of 1D heat conduction in a rod, written entirely in Julia using `GLMakie` for visualization and `DifferentialEquations.jl` for numerical integration.

Users can:
- Choose between several initial temperature profiles (constant, sine, Gaussian, linear gradient)
- Set Dirichlet or Neumann boundary conditions
- Select from a list of ODE solvers (RK4, Implicit Euler, Trapezoid, CNAB2, SBDF2, KenCarp4)
- Run the simulation and visualize the result as a time-animated temperature profile

This simulator is designed to demonstrate the physical behavior of heat transfer using PDE discretization and solver choice.

## Features
- Finite difference spatial discretization of the 1D heat equation
- Interactive parameter selection (material presets, time settings, BCs)
- Fully animated output using GLMakie sliders
- Optional data export
- Compatible with latest Julia version (v1.10+)

## Dependencies
Listed in `Project.toml`. Key packages:
- GLMakie
- DifferentialEquations
- ModelingToolkit
- NativeFileDialog
- DelimitedFiles

## Usage
1. Install Julia (1.10 or later).
2. Clone the repository.
3. Open `main.jl` in your Julia environment.
4. Run the script and configure the simulation via the GUI.

## Example
A user can model a copper rod with Dirichlet ends and a sine wave initial condition, then choose the Crank–Nicolson method for time integration. The resulting animation shows diffusion across time.

## License
MIT License

## Version
`v0.8.0` — Minor UI fixes and better time stepping. Initial public release.

