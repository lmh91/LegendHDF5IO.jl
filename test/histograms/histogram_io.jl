using StatsBase

using LegendHDF5IO: _nt_to_histogram, _histogram_to_nt

@testset "Histogram <-> NamedTuple" begin
    h = fit(Histogram, rand(10))
    @test _nt_to_histogram(_histogram_to_nt(h)) == h
    h = fit(Histogram, (rand(10), rand(10)))
    @test _nt_to_histogram(_histogram_to_nt(h)) == h
end