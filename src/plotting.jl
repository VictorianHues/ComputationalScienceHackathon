using ComputationalScienceHackathon

function plot_solution_animation(sol, params, file_name)
    @unpack N, x, zb = params

    anim = @animate for i in 1:2:length(sol.t)
        h = sol.u[i][1:N]
        zeta = h .+ zb 
        plot(x, zeta, ylim=(minimum(zb), 2.0),
            xlabel="x", ylabel="Water surface elevation ζ = h + zb",
            title="Time = $(round(sol.t[i], digits=2)) s", lw=2, label="Water Surface", legend=:bottomright)

        plot!(x, zb, label="Bed floor zb", linestyle=:dash, color=:black, legend=:bottomright)
    end

    file_path = joinpath("gifs", file_name)

    gif(anim, file_path, fps=15)
end