using Plots
using DifferentialEquations
using NonlinearSolve
using LinearAlgebra
using UnPack
using Sundials

function make_parameters()
    # Gravitational constant
    g = 9.81

    # Number of grid points
    nx = 200

    # Spatial domain
    x_start = 0.0
    x_stop = 5.0
    x = range(x_start, x_stop; length=nx)

    # Domain depth (flat bottom reference level)
    D = 10.0

    # Bed level
    zb = [-D + 0.4 * sin(2π * xi / x_stop * (nx - 1) / (nx * 5)) for xi in x]

    # Initial and final time
    t_initial = 0.0
    t_final = 1.0

    cf = 0.00232         # bed friction coefficient

    return (; g, N=nx, x, D, zb, t_initial, t_final, cf)
end

function initial_conditions(params)
    x = params.x
    zb = params.zb
    xmax = maximum(x)
    h = [0.1 * exp(-100 * ((xi / xmax - 0.5) * xmax)^2) - zbi for (xi, zbi) in zip(x, zb)]
    q = zeros(length(x))
    return h, q
end

function swe_dae_residual!(residual, du, u, p, t)
    @unpack g, N, x, zb, cf = p

    dx = x[2] - x[1]  # uniform grid spacing

    # Extract variables and their time derivatives
    h = @view u[1:N]
    q = @view u[N+1:2N]

    dhdt = @view du[1:N]
    dqdt = @view du[N+1:2N]

    # Fill residuals
    @views begin
        for i in 1:N
            # Periodic indexing
            iL = mod1(i - 1, N)
            iR = mod1(i + 1, N)

            # Mass conservation (continuity equation)
            dqdx = (q[iR] - q[iL]) / (2dx)
            residual[i] = dhdt[i] + dqdx

            # Momentum conservation
            q2h_R = q[iR]^2 / h[iR]
            q2h_L = q[iL]^2 / h[iL]
            dq2h_dx = (q2h_R - q2h_L) / (2dx)

            zeta_R = h[iR] + zb[iR]
            zeta_L = h[iL] + zb[iL]
            dzeta_dx = (zeta_R - zeta_L) / (2dx)

            friction = -cf * q[i] * abs(q[i]) / h[i]^2

            residual[N + i] = dqdt[i] + dq2h_dx + (g * h[i] * dzeta_dx) + friction
        end
    end

    return nothing
end

function timeloop(params)
    @unpack N, t_initial, t_final = params

    # Set up initial state vectors
    h0, q0 = initial_conditions(params)
    u0 = vcat(h0, q0)
    du0 = zeros(2N)  # Initial guess for du/dt

    # Time span for simulation
    tspan = (t_initial, t_final)

    # Specify that all variables are differential (time-dependent)
    differential_vars = trues(2N)

    # Define the DAE problem
    dae_prob = DAEProblem(
        swe_dae_residual!, du0, u0, tspan, params;
        differential_vars = differential_vars
    )

    cb = DiscreteCallback(log_tolerance_condition, log_tolerance_affect!)

    # Solve the problem using IDA solver (good for stiff DAE systems)
    #sol = solve(dae_prob, IDA(), reltol=1e-8, abstol=1e-8, callback=cb)
    sol = solve(dae_prob, IDA(), callback=cb)

    return sol  # return the solution object
end

function animate_solution(sol, params; filename="swe_animation.gif")
    N = params.N
    x = params.x

    anim = @gif for (i, t) in enumerate(sol.t)
        h = sol.u[i][1:N]
        plot(x, h, ylim=(0, 0.15),
             xlabel="x", ylabel="Water Depth h",
             title="t = $(round(t, digits=3)) s", lw=2)
    end every 2  # Save every 2nd time step to control GIF size

    # Save the animation silently
    gif(anim, filename, fps = 20)
end

# function main()
#     # Step 1: Load simulation parameters
#     parameters = make_parameters()
#     N = parameters.N

#     # Step 2: Generate initial conditions
#     h0, q0 = initial_conditions(parameters)
#     u0 = vcat(h0, q0)                  # State vector: [h; q]
#     du0 = zeros(length(u0))            # Initial guess for time derivative
#     tspan = (parameters.t_initial, parameters.t_final)

#     # Step 3: Define the DAE problem
#     differential_vars = trues(length(u0))  # all variables are differential
#     prob = DAEProblem(swe_dae_residual!, du0, u0, tspan, parameters; differential_vars=differential_vars)

#     # Step 4: Solve the problem using an implicit DAE solver
#     sol = solve(prob, IDA())

#     # Step 5: Plot h(x) at final time
#     h_final = sol[end][1:N]  # Extract water depth at final time
#     x = parameters.x
#     plot(x, h_final, xlabel="x", ylabel="Water Depth h", title="Final Water Depth h(x, t_final)", lw=2)
# end

function main()
    params = make_parameters()
    sol = timeloop(params)
    animate_solution(sol, params)  # Saves swe_animation.gif in current directory
end

main()
