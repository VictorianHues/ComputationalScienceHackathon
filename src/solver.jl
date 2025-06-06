using Plots

using ComputationalScienceHackathon

function solver()
    println("Solver Source File")
    value = 1
    println("Value: ", value)
    return value
end

function solve_jacobi(max_iter, n_steps, tolerance)
    c_old = zeros(n_steps, 10)  # Initialize concentration array
    c_new = zeros(n_steps, 10)  # New concentration array

    display(heatmap(c_old, title="Initial Concentration", xlabel="Column", ylabel="Row"))

    # Perform Jacobi iterations
    c_new, iteration, delta = jacobi(c_old, c_new, max_iter, n_steps, tolerance)

    display(heatmap(c_new, title="New Concentration", xlabel="Column", ylabel="Row"))

    return c_new, iteration, delta
end

function jacobi(c_old, c_new, max_iter, n_steps, tolerance)
    n_cols = size(c_old, 2)
    delta = 0.0
    iteration = 0
    for iter in 1:max_iter
        # Set boundary conditions
        c_old[:, 1] .= 0.0  # Left boundary
        c_old[:, end] .= 1.0  # Right boundary

        for i in 2:n_steps-1
            for j in 2:n_cols-1
                c_new[i, j] = 0.25 * (c_old[i+1, j] + c_old[i-1, j] + c_old[i, j+1] + c_old[i, j-1])
            end
        end

        # Enforce boundary conditions on c_new as well
        c_new[:, 1] .= 0.0
        c_new[:, end] .= 1.0

        delta = maximum(abs.(c_new - c_old))
        if delta < tolerance
            println("Converged after $iter iterations")
            iteration = iter
            break
        end
        c_old .= c_new
        iteration = iter
    end
    return c_new, iteration, delta
end