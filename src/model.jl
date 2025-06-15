using NonlinearSolve, LinearAlgebra, Parameters, Plots, Sundials

using ComputationalScienceHackathon

function timeloop(params)
    @unpack g, N, x, D, zb, tstart, tstop, h0, q0  = params

    u0 = vcat(h0, q0)
    du0 = zeros(2N)

    tspan = (tstart, tstop)

    differential_vars = trues(2N)

    dae_prob = DAEProblem(
        swe_dae_residual!, du0, u0, tspan, params;
        differential_vars=differential_vars
    )
    sol = solve(dae_prob, IDA())

    return sol
end