# ComputationalScienceHackathon

Repo for the Deltares Computational Science Hackathon (2025), hosted by [Computational Science NL](https://www.linkedin.com/company/computationalsciencenl) and featuring a 1D shallow water modeling challenge from [Deltares](https://www.deltares.nl/).

⚠️ _Note: This code was developed under time constraints typical of a hackathon setting. While functional, it is not fully refactored or optimized. See the `src/` directory for main scripts and solver logic._

## Simulations
### Periodic Boundaries
| Peak Water, Wavy Bottom without Discharge | Peak Water, Wavy Bottom with Horizontal Discharge |
|------------------------|----------------------------------|
| ![peak_water_wavy_bottom](https://github.com/VictorianHues/DeltaresHackathon2025/blob/main/gifs/swe_peak_water_wavy_bottom.gif?raw=true) | ![peak_water_wavy_bottom_discharge](https://github.com/VictorianHues/DeltaresHackathon2025/blob/main/gifs/swe_peak_water_wavy_bottom_discharge.gif?raw=true) |


| Peak Water, Peak Bottom without Discharge | Peak Water, Peak Bottom with Horizontal Discharge |
|------------------------|----------------------------------|
| ![swe_peak_water_peak_bottom](https://github.com/VictorianHues/DeltaresHackathon2025/blob/main/gifs/swe_peak_water_peak_bottom.gif?raw=true) | ![peak_water_wavy_bottom_discharge](https://github.com/VictorianHues/DeltaresHackathon2025/blob/main/gifs/swe_peak_water_peak_bottom_discharge.gif?raw=true) |

| Flat Water, Inclined Bottom with Horizontal Discharge | Flat Water, Wavy Bottom with Horizontal Discharge |
|------------------------|----------------------------------|
| ![swe_flat_water_incline_bottom_discharge](https://github.com/VictorianHues/DeltaresHackathon2025/blob/main/gifs/swe_flat_water_incline_bottom_discharge.gif?raw=true) | ![swe_flat_water_wavy_bottom_discharge](https://github.com/VictorianHues/DeltaresHackathon2025/blob/main/gifs/swe_flat_water_wavy_bottom_discharge.gif?raw=true) |



## Table of Contents

1. [Overview](#overview)
2. [Usage and Installation](#usage-and-installation)
3. [Implementation](#implementation)
4. [Contributing](#contributing)
5. [License](#license)

---

## Overview

This repository contains the code developed during the Deltares Computational Science Hackathon 2025. The challenge focused on simulating the **1D shallow water equations** over variable topography using **implicit time-stepping**.

The solver is implemented in **Julia** using the **SciML ecosystem**, particularly [`DAEProblem`](https://docs.sciml.ai/DiffEqDocs/stable/types/dae_types/) and the `IDA()` solver from Sundials. The core model includes:

- Central finite difference spatial discretization
- DAE formulation of the shallow water equations
- Implicit time integration via `IDA()`
- Support for periodic boundary conditions
- Animated visualization of wave evolution

---

## Usage and Installation

To run the simulations and generate plots, follow these steps:

1. Clone the repository:

    ```sh
    gh repo clone VictorianHues/DeltaresHackathon2025
    ```

2. Activate the Julia environment in the root directory:

    ```sh
    julia --project=.
    ```

3. **Install Dependencies**:

    From the Julia REPL:

    ```julia
    using Pkg
    Pkg.instantiate()
    Pkg.activate(".")
    ```

4. **Run the main script**: This will run the solver and generate an animated `.gif` of the wave surface evolution for a wide range of intial conditions and parameters. Remove and add portions of the `main()` function to select specific models.

    ```julia
    include("src/main.jl")
    ```

---

## Implementation

### Source Code

All implementation is located in the `src/` folder, including:

- `main.jl` - Main script to set parameters and initial conditions
- `solver.jl` - DAE residual function (`swe_dae_residual!`).
- `model.jl` - Timeloop function to model the solution over a specific timespan (`timeloop`)
- `plotting.jl` - Plotting functions used for animation and data analysis

### Features Implemented

- ✅ DAE formulation using `DAEProblem` with full residual specification
- ✅ Periodic and Dirichlet boundary conditions
- ✅ Initial conditions with Gaussian wave pulse and optional uniform flow
- ✅ Animated visualization of water surface

---

## Contributing

If you would like to contribute to this project, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Make your changes and commit them with clear and concise messages.
4. Push your changes to your forked repository.
5. Create a pull request to the main repository.

---

## License

This project is licensed under the MIT License. See the LICENSE file for more details.
