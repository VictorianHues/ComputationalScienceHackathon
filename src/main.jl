# Julia script to solve the 1D shallow water equations as a DAE problem
using NonlinearSolve, LinearAlgebra, Parameters, Plots, Sundials

using ComputationalScienceHackathon

function make_parameters(bed_level="flat")
    g = 9.81                     # gravitational acceleration
    N = 200                      # number of grid points
    x = range(0.0, 5; length=N) |> collect  # domain points
    xmax = maximum(x)            # maximum x in domain
    D = 10.0                     # domain depth

    tstart = 0.0
    tstop = 1.0

    cf = 0.00232  # friction coefficient

    # Bed level definition
    if bed_level == "flat"
        zb = fill(-D, N)  # flat bottom
    elseif bed_level == "wavy"
        zb = -D .+ 5.0 .* sin.(2*π*3 .* x ./ xmax .* ((N - 1) / N))
    elseif bed_level == "peak"
        zb = -D .+ 5.0 .* exp.(-0.5 .* ((x .- xmax/2) ./ (0.2 * xmax)).^2)  # normal distribution peak
    elseif bed_level == "incline"
        zb = -D .+ 5.0 .* (x ./ xmax)  # linear incline
    elseif bed_level == "decline"
        zb = -D .+ 5.0 .* (1 .- x ./ xmax)  # linear decline
    else
        error("Unknown bed level type: $bed_level")
    end

    return (; g, N, x, D, zb, tstart, tstop, cf)
end

function initial_conditions(params, u, h_initial="peak")
    @unpack N, x, zb = params
    xmax = maximum(x)

    # Initialize water height h
    if h_initial == "peak"
        h = 1.0 .* exp.(-100 .* ((x ./ xmax .- 0.5) .* xmax).^2) .- zb 
    elseif h_initial == "flat"
        h = fill(0.1, N) .- zb
    elseif h_initial == "wavy"
        h = 0.1 .+ 1.0 .* sin.(2*π*5 .* x ./ xmax .* ((N - 1) / N)) .- zb
    else
        error("Unknown initial condition type: $h_initial")
    end

    q = h .* u

    return h, q
end

function main()
    # Define experiment configurations
    bed_levels = ["wavy", "flat", "peak", "incline", "decline"]
    water_inits = ["peak", "flat"]
    boundary_conditions = ["dirichlet_neumann", "periodic"]
    u_initials = [0.0, 0.5]
    discharge_labels = ["", "_discharge"]

    for (i, u_initial) in enumerate(u_initials)
        for water_init in water_inits
            for bed_level in bed_levels
                for bc in boundary_conditions
                    # Only run "peak" water_init for u_initial == 0.0 and u_initial == 0.5
                    # Only run "flat" water_init for u_initial == 0.5
                    if (water_init == "flat" && u_initial == 0.0)
                        continue
                    end

                    discharge_label = discharge_labels[i]
                    file_name = "swe_$(water_init)_water_$(bed_level)_bottom_$(bc)$(discharge_label).gif"
                    params = make_parameters(bed_level)
                    h0, q0 = initial_conditions(params, u_initial, water_init)
                    params = merge(params, (; h0, q0))
                    solution = timeloop(params, bc)
                    plot_solution_animation(solution, params, file_name)
                end
            end
        end
    end


end

main()
