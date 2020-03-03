_range_to_nt(r::AbstractRange) = (
    start = first(r),
    stop = last(r),
    step = step(r),
    length = length(r)
)
_nt_to_range(nt::NamedTuple) =
    range(nt.start, nt.stop, length = nt.length)
    
_range_to_nt(r::UnitRange) = (
    start = first(r),
    stop = last(r)
)
_nt_to_range(nt::NamedTuple{(:start, :stop)}) =
    UnitRange(nt.start, nt.stop)

_range_to_nt(r::StepRange) = (
    start = first(r),
    step = step(r),
    stop = last(r)
)
_nt_to_range(nt::NamedTuple{(:start, :step, :stop)}) =
    StepRange(nt.start, nt.step, nt.stop)

_range_to_nt(r::LinRange) = (
    start = first(r),
    stop = last(r),
    length = length(r)
)
_nt_to_range(nt::NamedTuple{(:start, :stop, :length)}) =
    LinRange(nt.start, nt.stop, nt.length)


function _nt_to_histogram(nt::NamedTuple)
    return StatsBase.Histogram(
        tuple(map(b ->_nt_to_range(b.binedges), nt.binning)...), 
        nt.weights, 
        nt.binning[1].closedleft ? :left : :right, 
        nt.isdensity
    )  
end
    
function _histogram_to_nt(h::StatsBase.Histogram)
    n::Int = ndims(h.weights)
    axs_sym = Symbol.(["axis_$(i)" for i in Base.OneTo(n)])
    axs = [(
        binedges = _range_to_nt(h.edges[i]),
        closedleft = h.closed == :left 
    ) for i in Base.OneTo(n)]
    return (
        ndims = n,
        binning = NamedTuple{tuple(axs_sym...)}(axs),
        weights = h.weights,
        isdensity = h.isdensity,
        closed = h.closed    
    )
end
        

function LegendDataTypes.writedata(
    output::HDF5.DataFile, name::AbstractString,
    x::Histogram,
    fulldatatype::DataType = typeof(x)
) where {T}
    @assert fulldatatype == typeof(x)
    writedata(output, name, _histogram_to_nt(x))
end


function LegendDataTypes.readdata(
    input::HDF5.DataFile, name::AbstractString,
    AT::Type{<:Histogram}
)
    _nt_to_histogram(readdata(input, name))
end

