using NonlinearSolve, LinearAlgebra, Parameters, Plots, Sundials

using ComputationalScienceHackathon


# --- 4. Time integration ---
function timeloop(params)
    # Unpack parameters
    @unpack g, N, x, D, zb, tstart, tstop = params

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
    #sol = solve(dae_prob, IDA(), reltol=1e-8, abstol=1e-8) # solves the DAE problem using default settings
    sol = solve(dae_prob, FBDF(), reltol=1e-8, abstol=1e-8)

    # --- 5. a Live Plots ---
    anim = @animate for i in 1:10:length(sol.t)
        h = sol.u[i][1:N]
        plot(x, h, ylim=(9.5, 10.75), xlabel="x", ylabel="Water height h", title="Time = $(round(sol.t[i], digits=2)) s")
        plot!(x, zb .+ 20, label="Bed floor zb", linestyle=:dash, color=:black)
    end
    gif(anim, plot_name, fps=10)

    return sol # return solution object
end