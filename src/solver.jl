using NonlinearSolve, LinearAlgebra, Parameters, Plots, Sundials

using ComputationalScienceHackathon

function swe_dae_residual!(residual, du, u, p, t)
    @unpack g, N, x, zb, cf = p
    dx = x[2] - x[1]  # uniform grid assumed

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


function swe_dae_residual_dirichlet_neumann!(residual, du, u, p, t)
    @unpack g, N, x, zb, cf = p
    dx = x[2] - x[1]  # uniform grid assumed

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

    # Dirichlet boundary condition at the left end
    dqdx = (q[2] - q[1]) / dx
    rh[1] = dhdt[1] + dqdx
    dq2_over_h_dx = (q[2]^2 / h[2] - q[1]^2 / h[1]) / dx
    dzdx = (h[2] + zb[2] - h[1] - zb[1]) / dx
    friction = cf * q[1] * abs(q[1]) / h[1]^2
    rq[1] = dqdt[1] + dq2_over_h_dx + g * h[1] * dzdx + friction
    # Neumann boundary condition at the right end
    dqdx = (q[N] - q[N-1]) / dx
    rh[N] = dhdt[N] + dqdx
    dq2_over_h_dx = (q[N]^2 / h[N] - q[N-1]^2 / h[N-1]) / dx
    dzdx = (h[N] + zb[N] - h[N-1] - zb[N-1]) / dx
    friction = cf * q[N] * abs(q[N]) / h[N]^2
    rq[N] = dqdt[N] + dq2_over_h_dx + g * h[N] * dzdx + friction

    return nothing
end