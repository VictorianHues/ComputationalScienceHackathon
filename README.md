# ComputationalScienceHackathon

Repo for the Deltares Computational Science Hackathon

## Table of Contents

1. [Overview](#overview)
2. [Usage and Installation](#usage-and-installation)
3. [Implementation](#implementation)
4. [Contributing](#contributing)
5. [License](#license)

## Overview


## Usage and Installation

To run the simulations and generate plots, follow these steps:

1. Clone the repository:

    ```sh
    gh repo clone VictorianHues/ComputationalScienceHackathon
    ```

2. Activate the Julia Environment in the root directory:

    ```sh
    julia --project=.
    ```

3. **Install Dependencies**: Ensure you have all the required dependencies installed. You can install them from the Julia REPL:

    ```sh
    using Pkg
    Pkg.instantiate()
    Pkg.activate(".")
    Pkg.develop(path=".")
    ```

4. **Load the main script in the REPL**:

    ```sh
    include("src/main.jl")
    ```

5. **Run Main**:

    ```sh
    main()
    ```





8. **Run Unit Tests**:

    ```sh
    using Pkg
    Pkg.test()
    ```


## Implementation



### Source Code



### Scripts


## Contributing

If you would like to contribute to this project, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Make your changes and commit them with clear and concise messages.
4. Push your changes to your forked repository.
5. Create a pull request to the main repository.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.