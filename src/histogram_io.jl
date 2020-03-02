_nt_to_histogram(nt::NamedTuple) =
    StatsBase.Histogram(values(nt.edges), nt.weights, nt.closed, nt.isdensity)  
    
function _histogram_to_nt(h::StatsBase.Histogram) 
    n = ndims(h.weights)
    edges_vec = collect.(h.edges)
    edges_sym = Symbol.(["dim$(i)" for i in 1:n])
    nt_edges = NamedTuple{tuple(edges_sym...)}(edges_vec)
    return (
        n = ndims(h.weights),
        ndims = 3,
        weights = h.weights,
        edges = nt_edges,
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

