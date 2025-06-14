# Julia script to solve the 1D shallow water equations as a DAE problem
using NonlinearSolve, LinearAlgebra, Parameters, Plots, Sundials


# --- 1. Parameter setup ---
function make_parameters()
    g = 9.81                     # gravitational acceleration
    N = 200                      # number of grid points
    x = range(0.0, 5; length=N) |> collect  # domain points
    xmax = maximum(x)            # maximum x in domain
    D = 10.0                     # domain depth

    zb = -D .+ 0.4 .* sin.(2*π*5 .* x ./ xmax .* ((N - 1) / N))  # wavy bottom 
    # zb = fill(-D, N) # flat bottom

    tstart = 0.0
    tstop = 1.0

    cf = 0.00232  # friction coefficient

    return (; g, N, x, D, zb, tstart, tstop, cf)
end

# --- 2. Initial condition ---
function initial_conditions(params)
    @unpack N, x, zb = params
    xmax = maximum(x)
    h = 0.1 .* exp.(-100 .* ((x ./ xmax .- 0.5) .* xmax).^2) .- zb  # from Table I
    q = zeros(N) 

    return h, q
end

function main()
    params = make_parameters()
    
    solution = timeloop(params) 
    plot_solution(solution, params)
end

main()
