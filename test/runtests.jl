using CellMetabolism
using Test
using Aqua
using JET

@testset "CellMetabolism.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(CellMetabolism)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(CellMetabolism; target_defined_modules = true)
    end
    # Write your tests here.
end
