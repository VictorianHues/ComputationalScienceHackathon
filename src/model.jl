using NonlinearSolve, LinearAlgebra, Parameters, Plots, Sundials

using ComputationalScienceHackathon

function timeloop(params, boundary_conditions="dirichlet_neumann")
    @unpack g, N, x, D, zb, tstart, tstop, h0, q0  = params

    u0 = vcat(h0, q0)
    du0 = zeros(2N)

    tspan = (tstart, tstop)

    differential_vars = trues(2N)

    if boundary_conditions == "dirichlet_neumann"
        dae_prob = DAEProblem(
        swe_dae_residual_dirichlet_neumann!, du0, u0, tspan, params;
        differential_vars=differential_vars
        )
    elseif boundary_conditions == "periodic"
        dae_prob = DAEProblem(
            swe_dae_residual!, du0, u0, tspan, params;
            differential_vars=differential_vars
        )
    else
        error("Unknown boundary condition type: $boundary_conditions")
    end

    
    sol = solve(dae_prob, IDA())

    return sol
end