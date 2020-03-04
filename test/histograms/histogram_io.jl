using StatsBase

using LegendHDF5IO: _nt_to_histogram, _histogram_to_nt

@testset "Histogram <-> NamedTuple" begin
    h = fit(Histogram, rand(10))
    @test _nt_to_histogram(_histogram_to_nt(h)) == h
    h = fit(Histogram, (rand(10), rand(10)))
    @test _nt_to_histogram(_histogram_to_nt(h)) == h
    h = fit(Histogram, (rand(10), rand(10)), (0:0.2:1, Float64[0, 0.5, 1]))
    @test _nt_to_histogram(_histogram_to_nt(h)) == h
end