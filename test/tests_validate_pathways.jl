@testitem "MetabolicPathways Validation Tests" begin
    using LabelledArrays, BenchmarkTools, OrdinaryDiffEq, CellMetabolismBase

@test_nowarn CellMetabolismBase.validate_MetabolicPathway(
    CellMetabolism.glycolysis_pathway,
    CellMetabolism.glycolysis_init_conc,
    CellMetabolism.glycolysis_params,
    )
end
