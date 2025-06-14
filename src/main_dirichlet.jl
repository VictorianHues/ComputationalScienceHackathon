# Julia script to solve the 1D shallow water equations as a DAE problem
using NonlinearSolve, LinearAlgebra, Parameters, Plots, Sundials

# using ComputationalScienceHackathon

# --- 1. Parameter setup ---
function make_parameters()
    g = 9.81                     # gravitational acceleration
    N = 200                      # number of grid points
    x = range(0.0, 5; length=N) |> collect  # domain points
    xmax = maximum(x)            # maximum x in domain
    D = 10.0                     # domain depth

    #zb = -D .+ 0.4 .* sin.(2*π*5 .* x ./ xmax .* ((N - 1) / N))  # wavy bottom 
    zb = fill(-D, N) # flat bottom

    tstart = 0.0
    tstop = 1.0

    cf = 0.00232  # friction coefficient

    h = 0.1 .* exp.(-100 .* ((x ./ xmax .- 0.5) .* xmax).^2) .- zb  # from Table I
    q = zeros(N) 

    h_bc = [h[1], h[end]]   # h at left and right boundaries
    q_bc = [q[1], q[end]]   # q at left and right boundaries



    return (; g, N, x, D, zb, tstart, tstop, cf, h_bc, q_bc)
end

# --- 2. Initial condition ---
function initial_conditions(params)
    @unpack N, x, zb = params
    xmax = maximum(x)
    h = 0.1 .* exp.(-100 .* ((x ./ xmax .- 0.5) .* xmax).^2) .- zb  # from Table I
    q = zeros(N) 

    return h, q
end

function swe_dae_residual!(residual, du, u, p, t)
    @unpack g, N, x, zb, cf, h_bc, q_bc = p  # Add h_bc, q_bc for Dirichlet BCs
    dx = x[2] - x[1]  # uniform grid assumed

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
        dqdx = (q[i+1] - q[i-1]) / (2*dx)
        rh[i] = dhdt[i] + dqdx

        # Equation 5b
        dq2_over_h_dx = ((q[i+1]^2 / h[i+1]) - (q[i-1]^2 / h[i-1])) / (2*dx)
        dzdx = ((h[i+1] + zb[i+1]) - (h[i-1] + zb[i-1])) / (2*dx)
        friction = cf * q[i] * abs(q[i]) / h[i]^2

        rq[i] = dqdt[i] + dq2_over_h_dx + g * h[i] * dzdx + friction
    end

    # --- Boundary Conditions: Dirichlet (fixed h and q at boundaries) ---
    # i = 1 (left boundary)
    rh[1] = h[1] - h_bc[1]      # Dirichlet for h at left
    rq[1] = q[1] - q_bc[1]      # Dirichlet for q at left

    # i = N (right boundary)
    rh[N] = h[N] - h_bc[2]      # Dirichlet for h at right
    rq[N] = q[N] - q_bc[2]      # Dirichlet for q at right

    return nothing
end

function timeloop(params)
    @unpack g, N, x, D, zb, tstart, tstop = params

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
    sol = solve(dae_prob, IDA()) # solves the DAE problem using default settings
    #sol = solve(dae_prob, Trapezoid(), reltol=1e-8, abstol=1e-8)

    # --- 5. a Live Plots ---
    anim = @animate for i in 1:10:length(sol.t)
        h = sol.u[i][1:N]
        plot(x, h, ylim=(9.5, 10.75), xlabel="x", ylabel="Water height h", title="Time = $(round(sol.t[i], digits=2)) s")
        plot!(x, zb .+ 20, label="Bed floor zb", linestyle=:dash, color=:black)
    end
    gif(anim, "dirichlet_flat.gif", fps=10)

    return sol # return solution object
end

function main()
    params = make_parameters()
    
    solution = timeloop(params) 
    #plot_solution(solution, params)
end

main()
