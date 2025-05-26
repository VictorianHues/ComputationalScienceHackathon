using Test
using ComputationalScienceHackathon

@testset "Solver Tests" begin
    expected_value = 1
    result = solver()
    println("Solver result: ", result)
    @test result == expected_value
end