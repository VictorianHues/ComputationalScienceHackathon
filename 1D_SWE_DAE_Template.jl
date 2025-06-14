# Julia script to solve the 1D shallow water equations as a DAE problem
using NonlinearSolve, LinearAlgebra, Parameters, Plots, Sundials

# This code is a template for solving the 1D shallow water equations (SWE) as a DAE problem.
# The setup is as follows:
# 1. Define parameters for the simulation, including gravity, number of grid points,
#   spatial domain, and bottom topography.
# 2. Set up initial conditions for water height and momentum.
# 3. Define the DAE residual function that describes the SWE.
# 4. Implement a time loop to solve the DAE problem using Sundials' IDA solver.
# 5. Plot the results.
# 6. Function calls to start the simulation.

#Boundry option choice 
const CONFIG = Dict{Symbol, Any}()

function set_boundary!(type::Symbol; q_left=0.1, h_right=10.0, robin_alpha=1.0, robin_target=10.0)
    CONFIG[:boundary_type] = type
    if type == :dirichlet_neumann
        CONFIG[:q_left] = q_left
        CONFIG[:robin_alpha] = robin_alpha
        CONFIG[:robin_target] = robin_target
    elseif type == :periodic
        # Optional: clean up any previous values
        pop!(CONFIG, :q_left, nothing)
        pop!(CONFIG, :h_right, nothing)
    else
        error("Unsupported boundary type: $type")
    end
end

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

    return (; g, N, x, D, zb, tstart, tstop)
end

# --- 2. Initial condition ---
function initial_conditions(params)
    @unpack N, x, zb = params
    xmax = maximum(x)
    h = 0.1 .* exp.(-100 .* ((x ./ xmax .- 0.5) .* xmax).^2) .- zb  # from Table I
    q = zeros(N)  # initially at rest

    return h, q
end

# --- 3. DAE residual function ---
# Note: the "!" at the end of the function name indicates that the function modifies 
# its arguments (convention in Julia)
function swe_dae_residual!(residual, du, u, p, t)
    @unpack g, N, x, zb = p
    dx = x[2] - x[1]  # uniform grid assumed
    cf = 0.00232  # friction coefficient

    # Extract h and q from state vector u
    h = @view u[1:N]
    q = @view u[N+1:2N]

    dhdt = @view du[1:N]
    dqdt = @view du[N+1:2N]

    # Allocate views into residual
    rh = @view residual[1:N]
    rq = @view residual[N+1:2N]

    # Apply interior stencil (2nd-order central difference)
    for i in 2:N-1
        # Equation 5a: ∂h/∂t + ∂q/∂x = 0
        ∂q∂x = (q[i+1] - q[i-1]) / (2*dx)
        rh[i] = dhdt[i] + ∂q∂x

        # Equation 5b
        dq2_over_h_dx = ((q[i+1]^2 / h[i+1]) - (q[i-1]^2 / h[i-1])) / (2*dx)
        ∂ζ∂x = ((h[i+1] + zb[i+1]) - (h[i-1] + zb[i-1])) / (2*dx)
        friction = cf * q[i] * abs(q[i]) / h[i]^2

        rq[i] = dqdt[i] + dq2_over_h_dx + g * h[i] * ∂ζ∂x + friction
    end

    # --- Boundary Conditions: Periodic or Neumann/Dirichlet---
    
    boundary_type = CONFIG[:boundary_type]

    if boundary_type == :periodic
        # i = 1 (left boundary)
        ∂q∂x = (q[2] - q[N]) / (2*dx)
        rh[1] = dhdt[1] + ∂q∂x

        dq2_over_h_dx = ((q[2]^2 / h[2]) - (q[N]^2 / h[N])) / (2*dx)
        ∂ζ∂x = ((h[2] + zb[2]) - (h[N] + zb[N])) / (2*dx)
        friction = cf * q[1] * abs(q[1]) / h[1]^2
        rq[1] = dqdt[1] + dq2_over_h_dx + g * h[1] * ∂ζ∂x + friction

        # i = N (right boundary)
        ∂q∂x = (q[1] - q[N-1]) / (2*dx)
        rh[N] = dhdt[N] + ∂q∂x

        dq2_over_h_dx = ((q[1]^2 / h[1]) - (q[N-1]^2 / h[N-1])) / (2*dx)
        ∂ζ∂x = ((h[1] + zb[1]) - (h[N-1] + zb[N-1])) / (2*dx)
        friction = cf * q[N] * abs(q[N]) / h[N]^2
        rq[N] = dqdt[N] + dq2_over_h_dx + g * h[N] * ∂ζ∂x + friction

    elseif boundary_type == :dirichlet_neumann
        q_left = CONFIG[:q_left]
        α = CONFIG[:robin_alpha]
        b = CONFIG[:robin_target]

        # --- Left boundary ---
        ∂q∂x = (q[2] - q_left) / (2dx)
        rh[1] = dhdt[1] + ∂q∂x
        dq2_over_h_dx = ((q[2]^2 / h[2]) - (q_left^2 / h[1])) / (2dx)
        ∂ζ∂x = ((h[2] + zb[2]) - (h[1] + zb[1])) / dx
        friction = cf * q[1] * abs(q[1]) / h[1]^2
        rq[1] = dqdt[1] + dq2_over_h_dx + g * h[1] * ∂ζ∂x + friction

        # --- Right boundary ---
        rh[N] = h[N] + α * (h[N] - h[N-1]) / dx - b

        # Use one-sided spatial diff for q at right boundary
        ∂q∂x = (q[N] - q[N-1]) / dx
        dq2_over_h_dx = ((q[N]^2 / h[N]) - (q[N-1]^2 / h[N-1])) / dx
        ∂ζ∂x = ((h[N] + zb[N]) - (h[N-1] + zb[N-1])) / dx
        friction = cf * q[N] * abs(q[N]) / h[N]^2
        rq[N] = dqdt[N] + dq2_over_h_dx + g * h[N] * ∂ζ∂x + friction
    else
        error("Unsupported boundary type: $boundary_type")
    end
end

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
    sol = solve(dae_prob, IDA(), reltol=1e-8, abstol=1e-8) # solves the DAE problem using default settings

    # --- 5. a Live Plots ---
    anim = @animate for i in 1:10:length(sol.t)
        h = sol.u[i][1:N]
        plot(x, h, ylim=(-12, 12), xlabel="x", ylabel="Water height h", title="Time = $(round(sol.t[i], digits=2)) s")
        plot!(x, zb, label="Bed floor zb", linestyle=:dash, color=:black)
    end
    gif(anim, "h_evolution_D_N.gif", fps=10)

    return sol # return solution object
end

# --- 5. b Plotting results ---
function plot_solution(sol, params)
    @unpack N, x = params
    h_series = [sol.u[i][1:N] for i in 1:length(sol.t)]

    anim = @animate for i in 1:10:length(sol.t)
        plot(x, h_series[i], ylim=(0, 12), xlabel="x", ylabel="Water height h",
             title="Time = $(round(sol.t[i], digits=2)) s", legend=false)
    end
    gif(anim, "shallow_water_D_N.gif", fps=10)
end

# --- 6. Main script ---
#choose the boundry
set_boundary!(:dirichlet_neumann, q_left=0.1, h_right=10.0, robin_alpha=1.0, robin_target=10.0)
# Set up parameters
params = make_parameters()
# Call the time loop function
solution = timeloop(params) 
plot_solution(solution, params)
