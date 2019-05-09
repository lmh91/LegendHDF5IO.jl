# This file is a part of LegendHDF5IO.jl, licensed under the MIT License (MIT


function to_table(x::AbstractVector{<:RDWaveform})
    TypedTables.Table(
        t0 = first.(x.t),
        dt = step.(x.t),
        value = x.v
    )
end


function _dtt02range(dt::RealQuantity, t0::RealQuantity, value::AbstractVector)
    # TODO: Handle different units for dt and t0
    t0 .+ (0:(size(value, 1) - 1)) .* dt
end

function from_table(tbl, ::Type{<:AbstractVector{<:RDWaveform}})
    StructArray{RDWaveform}((
        tbl.value,
        _dtt02range.(tbl.dt, tbl.t0, tbl.value)
    ))
end


function LegendDataTypes.writedata(
    output::HDF5.DataFile, name::AbstractString,
    x::AbstractVector{<:RDWaveform},
    fulldatatype::DataType = typeof(x)
) where {T}
    @assert fulldatatype == typeof(x)
    writedata(output, name, to_table(x))
end


function LegendDataTypes.readdata(
    input::HDF5.DataFile, name::AbstractString,
    AT::Type{<:AbstractVector{<:RDWaveform}}
)
    tbl = readdata(input, name, TypedTables.Table{<:NamedTuple{(:t0, :dt, :value)}})
    from_table(tbl, AbstractVector{<:RDWaveform})
end

