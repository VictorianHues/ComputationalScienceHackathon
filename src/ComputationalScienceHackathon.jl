"""
Module Entry Point for Computational Science Hackathon Project
This module includes various components such as utilities, solvers, models, and plotting functionalities.
"""
module ComputationalScienceHackathon

# Include additional source files
include("utils.jl")
include("solver.jl")
include("model.jl")
include("plotting.jl")

# Optionally, export functions or types you want users to access directly
# export my_function, MyType

export solver, model, plotting, solve_jacobi, solve_groundwater_1d

end