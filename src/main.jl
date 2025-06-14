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

function swe_dae_residual!(residual, du, u, p, t)
    @unpack g, N, x, zb, cf = p
    dx = x[2] - x[1]  # uniform grid assumed
    #cf = 0.00232  # friction coefficient

    h = @view u[1:N]
    q = @view u[N+1:2N]

    dhdt = @view du[1:N]
    dqdt = @view du[N+1:2N]

    rh = @view residual[1:N]
    rq = @view residual[N+1:2N]

    for i in 2:N-1
        dqdx = (q[i+1] - q[i-1]) / (2*dx)
        rh[i] = dhdt[i] + dqdx

        dq2_over_h_dx = ((q[i+1]^2 / h[i+1]) - (q[i-1]^2 / h[i-1])) / (2*dx)
        dzdx = ((h[i+1] + zb[i+1]) - (h[i-1] + zb[i-1])) / (2*dx)
        friction = cf * q[i] * abs(q[i]) / h[i]^2

        rq[i] = dqdt[i] + dq2_over_h_dx + g * h[i] * dzdx + friction
    end

    dqdx = (q[2] - q[N]) / (2*dx)
    rh[1] = dhdt[1] + dqdx

    dq2_over_h_dx = ((q[2]^2 / h[2]) - (q[N]^2 / h[N])) / (2*dx)
    dzdx = ((h[2] + zb[2]) - (h[N] + zb[N])) / (2*dx)
    friction = cf * q[1] * abs(q[1]) / h[1]^2
    rq[1] = dqdt[1] + dq2_over_h_dx + g * h[1] * dzdx + friction

    dqdx = (q[1] - q[N-1]) / (2*dx)
    rh[N] = dhdt[N] + dqdx

    dq2_over_h_dx = ((q[1]^2 / h[1]) - (q[N-1]^2 / h[N-1])) / (2*dx)
    dzdx = ((h[1] + zb[1]) - (h[N-1] + zb[N-1])) / (2*dx)
    friction = cf * q[N] * abs(q[N]) / h[N]^2
    rq[N] = dqdt[N] + dq2_over_h_dx + g * h[N] * dzdx + friction

    return nothing
end

function timeloop(params)
    @unpack g, N, x, D, zb, tstart, tstop = params

    h0, q0 = initial_conditions(params)
    u0 = vcat(h0, q0)
    du0 = zeros(2N)

    tspan = (tstart, tstop)

    # Specify differentiable variables as (true) -> all variables
    differential_vars = trues(2N)

    dae_prob = DAEProblem(
        swe_dae_residual!, du0, u0, tspan, params;
        differential_vars=differential_vars
    )
    sol = solve(dae_prob, IDA()) # solves the DAE problem using default settings

    # --- 5. a Live Plots ---
    anim = @animate for i in 1:10:length(sol.t)
        h = sol.u[i][1:N]
        plot(x, h, ylim=(9.5, 10.75), xlabel="x", ylabel="Water height h", title="Time = $(round(sol.t[i], digits=2)) s")
        plot!(x, zb .+ 20, label="Bed floor zb", linestyle=:dash, color=:black)
    end
    gif(anim, "periodic.gif", fps=10)

    return sol # return solution object
end

function compute_total_mass(sol, params)
    @unpack N, x = params
    dx = x[2] - x[1]
    mass = [sum(sol.u[i][1:N]) * dx for i in eachindex(sol.t)]
    return sol.t, mass
end

function compute_total_energy(sol, params)
    @unpack N, x, g = params
    dx = x[2] - x[1]

    energy = Float64[]
    for u in sol.u
        h = u[1:N]
        q = u[N+1:2N]
        E_kin = q.^2 ./ (2 .* h)
        E_pot = 0.5 .* g .* h.^2
        push!(energy, sum(E_kin .+ E_pot) * dx)
    end
    return sol.t, energy
end

function plot_final_state(sol, params)
    @unpack x, zb, N = params
    h_final = sol.u[end][1:N]
    q_final = sol.u[end][N+1:2N]

    plt = plot(x, h_final, label="Water height h", lw=2)
    plot!(plt, x, zb .+ 20, label="Bed elevation", linestyle=:dash)
    plot!(plt, xlabel="x", ylabel="h", title="Final State at t = $(round(sol.t[end], digits=2)) s")
    savefig(plt, "final_state.png")
end

function plot_conservation(sol, params)
    t1, mass = compute_total_mass(sol, params)
    t2, energy = compute_total_energy(sol, params)

    p1 = plot(t1, mass, lw=2, xlabel="Time", ylabel="Total Mass", title="Mass Conservation")
    p2 = plot(t2, energy, lw=2, xlabel="Time", ylabel="Total Energy", title="Energy Conservation")
    plt = plot(p1, p2, layout=(1,2), size=(900,400))
    savefig(plt, "conservation.png")
end

function main()
    params = make_parameters()
    
    solution = timeloop(params) 

    plot_final_state(solution, params)
    plot_conservation(solution, params)
    #plot_solution(solution, params)
end

main()
