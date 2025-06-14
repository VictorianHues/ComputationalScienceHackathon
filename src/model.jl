using Plots
using DifferentialEquations
using NonlinearSolve
using LinearAlgebra
using UnPack

using ComputationalScienceHackathon


function timeloop(params)
    # Unpack parameters 

    # set up initial conditions
    h0, q0 = initial_conditions(params)
    u0 = vcat(h0, q0)
    du0 = zeros(2N)  # Initial guess for du/dt

    tspan = (tstart, tstop) # defines the start and end times for the simulation

    # Specify differentiable variables as (true) -> all variables
    differential_vars = trues(2N)

    dae_prob = DAEProblem(
        swe_dae_residual!, du0, u0, tspan, params;
        differential_vars=differential_vars
    )
    sol = solve(dae_prob, IDA(), reltol=1e-8, abstol=1e-8) # solves the DAE problem using default settings

    # --- 5. a Live Plots ---

    return sol # return solution object
end