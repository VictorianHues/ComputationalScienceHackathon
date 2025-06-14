using ComputationalScienceHackathon

function plot_solution(sol, params)
    @unpack N, x, zb = params
    h_series = [sol.u[i][1:N] for i in 1:length(sol.t)]

    # anim = @animate for i in 1:10:length(sol.t)
    anim = @animate for i in 1:10:length(sol.t)
        plot(x, h_series[i], ylim=(9.5, 10.5), xlabel="x", ylabel="Water height h",
             title="Time = $(round(sol.t[i], digits=2)) s", legend=false)
        plot!(x, zb .+ 20, label="Bed floor zb", linestyle=:dash, color=:black)
    end
    gif(anim, "shallow_water.gif", fps=10)
end