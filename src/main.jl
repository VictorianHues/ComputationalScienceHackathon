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

    if bed_level == "flat"
        zb = fill(-D, N)  # flat bottom
    elseif bed_level == "wavy"
        zb = -D .+ 1.0 .* sin.(2*π*5 .* x ./ xmax .* ((N - 1) / N))
    elseif bed_level == "peak"
        zb = -D .+ 1.0 .* exp.(-0.5 .* ((x .- xmax/2) ./ (0.2 * xmax)).^2)  # normal distribution peak
    elseif bed_level == "incline"
        zb = -D .+ 1.0 .* (x ./ xmax)  # linear incline
    elseif bed_level == "decline"
        zb = -D .+ 1.0 .* (1 .- x ./ xmax)  # linear decline
    else
        error("Unknown bed level type: $bed_level")
    end

    return (; g, N, x, D, zb, tstart, tstop, cf)
end

function initial_conditions(params, u, h_initial="peak")
    @unpack N, x, zb = params
    xmax = maximum(x)

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
    u_initial = 0.0
    file_name = "swe_single_peak_wavy_bottom.gif"
    params = make_parameters("wavy")
    h0, q0 = initial_conditions(params, u_initial, "peak")
    params = merge(params, (; h0, q0))
    solution = timeloop(params) 
    plot_solution_animation(solution, params, file_name)

    file_name = "swe_single_peak_flat_bottom.gif"
    params = make_parameters("flat")
    params = merge(params, (; h0, q0))
    solution = timeloop(params) 
    plot_solution_animation(solution, params, file_name)

    file_name = "swe_single_peak_peak_bottom.gif"
    params = make_parameters("peak")
    params = merge(params, (; h0, q0))
    solution = timeloop(params) 
    plot_solution_animation(solution, params, file_name)

    file_name = "swe_single_peak_incline_bottom.gif"
    params = make_parameters("incline")
    params = merge(params, (; h0, q0))
    solution = timeloop(params) 
    plot_solution_animation(solution, params, file_name)

    file_name = "swe_single_peak_decline_bottom.gif"
    params = make_parameters("decline")
    params = merge(params, (; h0, q0))
    solution = timeloop(params) 
    plot_solution_animation(solution, params, file_name)


    u_initial = 0.5
    file_name = "swe_single_flat_wavy_bottom_discharge.gif"
    params = make_parameters("wavy")
    h0, q0 = initial_conditions(params, u_initial, "flat")
    params = merge(params, (; h0, q0))
    solution = timeloop(params) 
    plot_solution_animation(solution, params, file_name)

    file_name = "swe_single_flat_flat_bottom_discharge.gif"
    params = make_parameters("flat")
    params = merge(params, (; h0, q0))
    solution = timeloop(params) 
    plot_solution_animation(solution, params, file_name)

    file_name = "swe_single_flat_peak_bottom_discharge.gif"
    params = make_parameters("peak")
    params = merge(params, (; h0, q0))
    solution = timeloop(params) 
    plot_solution_animation(solution, params, file_name)

    file_name = "swe_single_flat_incline_bottom_discharge.gif"
    params = make_parameters("incline")
    params = merge(params, (; h0, q0))
    solution = timeloop(params)
    plot_solution_animation(solution, params, file_name)
    file_name = "swe_single_flat_decline_bottom_discharge.gif"
    params = make_parameters("decline")
    params = merge(params, (; h0, q0))
    solution = timeloop(params)
    plot_solution_animation(solution, params, file_name)


    u_initial = 0.5
    file_name = "swe_single_peak_wavy_bottom_discharge.gif"
    params = make_parameters("wavy")
    h0, q0 = initial_conditions(params, u_initial, "peak")
    params = merge(params, (; h0, q0))
    solution = timeloop(params) 
    plot_solution_animation(solution, params, file_name)

    file_name = "swe_single_peak_flat_bottom_discharge.gif"
    params = make_parameters("flat")
    params = merge(params, (; h0, q0))
    solution = timeloop(params) 
    plot_solution_animation(solution, params, file_name)

    file_name = "swe_single_peak_peak_bottom_discharge.gif"
    params = make_parameters("peak")
    params = merge(params, (; h0, q0))
    solution = timeloop(params) 
    plot_solution_animation(solution, params, file_name)

    file_name = "swe_single_peak_incline_bottom_discharge.gif"
    params = make_parameters("incline")
    params = merge(params, (; h0, q0))
    solution = timeloop(params)
    plot_solution_animation(solution, params, file_name)
    file_name = "swe_single_peak_decline_bottom_discharge.gif"
    params = make_parameters("decline")
    params = merge(params, (; h0, q0))
    solution = timeloop(params)
    plot_solution_animation(solution, params, file_name)


end

main()
