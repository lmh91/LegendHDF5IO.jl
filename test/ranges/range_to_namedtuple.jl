using LegendHDF5IO: _nt_to_range, _range_to_nt

@testset "Range <-> NamedTuple" begin
    @testset "UnitRange" begin
        @test _nt_to_range(_range_to_nt(1:10)) == 1:10
    end
    @testset "StepRange" begin
        @test _nt_to_range(_range_to_nt(1:2:10)) == 1:2:10
    end
    @testset "LinRange" begin
        @test _nt_to_range(_range_to_nt(LinRange(1, 10, 3))) == LinRange(1, 10, 3)
    end
    @testset "StepRangeLen" begin
        @test _nt_to_range(_range_to_nt(StepRangeLen(1.0:0.3:1.7))) == StepRangeLen(1.0:0.3:1.7)
    end
end