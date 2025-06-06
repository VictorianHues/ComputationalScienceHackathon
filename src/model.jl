using LinearAlgebra

using ComputationalScienceHackathon

function model()
    println("Model Source File")

end

# Function to solve 1D steady groundwater flow
function solve_groundwater_1d(K::Float64, h0::Float64, hL::Float64, L::Float64, N::Int)
    Δx = L / (N + 1)
    x = range(0, L, length=N+2)

    # Construct system: A * h = b
    A = zeros(N, N)
    b = zeros(N)

    for i in 1:N
        A[i, i] = -2
        if i > 1
            A[i, i-1] = 1
        end
        if i < N
            A[i, i+1] = 1
        end
    end

    # Incorporate Dirichlet BCs into b
    b[1] -= h0
    b[end] -= hL

    # Solve system
    h_internal = A \ b

    # Combine full solution
    h = vcat(h0, h_internal, hL)
    return x, h
end