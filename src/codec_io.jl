# This file is a part of LegendHDF5IO.jl, licensed under the MIT License (MIT


# ToDo: Generalize to any Codec, currently specific to RadwareSigcompress
function to_table(x::AbstractVector{<:EncodedData})
    TypedTables.Table(
        bytes = x.bytes,
        codec_shift = (c -> c.shift).(x.codec)
    )
end


# ToDo: Generalize to any Codec, currently specific to RadwareSigcompress
function from_table(tbl, ::Type{<:AbstractVector{<:EncodedData}})
    codec_vec = RadwareSigcompress.(tbl.codec_shift)
    bytes_vec = tbl.bytes
    StructArray{EncodedData{Vector{Int32,},eltype(codec_vec),Vector{UInt8}}}((codec_vec, bytes_vec))
end


function LegendDataTypes.writedata(
    output::HDF5.DataFile, name::AbstractString,
    x::AbstractVector{<:EncodedData},
    fulldatatype::DataType = typeof(x)
) where {T}
    @assert fulldatatype == typeof(x)
    writedata(output, name, to_table(x))
end


function LegendDataTypes.readdata(
    input::HDF5.DataFile, name::AbstractString,
    AT::Type{<:AbstractVector{<:EncodedData}}
)
    tbl = readdata(input, name, TypedTables.Table{<:NamedTuple{(:bytes, :codec_shift)}})
    from_table(tbl, AbstractVector{<:EncodedData})
end
