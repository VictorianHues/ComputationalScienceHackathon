using Plots

using ComputationalScienceHackathon

function main()
    println("Starting Computational Science Hackathon project...")

    solver()

    max_iter = 1000
    n_steps = 100
    tolerance = 1e-6
    c_new, iteration, delta = solve_jacobi(max_iter, n_steps, tolerance)

    println("Jacobi solver completed:")
    println("Final concentration array: ", c_new)
    println("Iterations: ", iteration)
    println("Final delta: ", delta)
    println("Project completed successfully!")





    # Parameters
    K = 10.0       # Hydraulic conductivity [m/day]
    h0 = 100.0     # Head at x = 0
    hL = 90.0      # Head at x = L
    L = 100.0      # Length of domain [m]
    N = 50         # Number of internal grid points

    # Solve and plot
    x, h = solve_groundwater_1d(K, h0, hL, L, N)
    plot(x, h, xlabel="x (m)", ylabel="Head h (m)", lw=2, title="1D Steady Groundwater Flow")


end

main()