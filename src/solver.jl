using Plots
using DifferentialEquations
using NonlinearSolve
using LinearAlgebra
using UnPack

using ComputationalScienceHackathon


function swe_dae_residual!(residual, du, u, p, t)
    @unpack g, N, x, D, zb, t_initial, t_final, cf = p
    dx = x[2] - x[1]  # uniform grid assumed

    h = @view u[1:N]
    q = @view u[N+1:2N]

    dhdt = @view du[1:N]
    dqdt = @view du[N+1:2N]

    @views begin
        for i in 1:N
            # Periodic boundary indices
            iL = mod1(i - 1, N)
            iR = mod1(i + 1, N)

            # Mass conservation
            dqdx = (q[iR] - q[iL]) / (2 * dx)
            residual[i] = dhdt[i] + dqdx

            # Momentum conservation
            dq2h_dx = ((q[iR]^2 / h[iR]) - (q[iL]^2 / h[iL])) / (2 * dx)
            dzetadx = ((h[iR] + zb[iR]) - (h[iL] + zb[iL])) / (2 * dx)
            shear = -cf * q[i] * abs(q[i]) / h[i]^2

            residual[N + i] = dqdt[i] + dq2h_dx + g * h[i] * dzetadx + shear
        end
    end

    return nothing
end