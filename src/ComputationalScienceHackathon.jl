module ComputationalScienceHackathon

# Include additional source files
include("utils.jl")
include("solver.jl")
include("model.jl")
include("plotting.jl")

# Optionally, export functions or types you want users to access directly
# export my_function, MyType

export solver, model, plotting

end