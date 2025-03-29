@testitem "Code quality (Aqua.jl)" begin
    using Aqua
    Aqua.test_all(CellMetabolism)
end
@testitem "Code linting (JET.jl)" begin
    using JET
    JET.test_package(CellMetabolism; target_defined_modules=true)
end
