# This file is a part of LegendHDF5IO.jl, licensed under the MIT License (MIT


#!!! TODO: LegendHDF5File, LegendHDF5Input, LegendHDF5Output
# with Base.read/write, open/close, etc., and atomic file names.

# use :semantictype attribute for high-level meaning of datatypes


const default_compression = ("shuffle", (), "deflate", 3)


const datatype_regexp = r"""^(([A-Za-z_]*)(<([0-9,]*)>)?)(\{(.*)\})?$"""
const arraydims_regexp = r"""^<([0-9,]*)>$"""

function _eldatatype_from_string(s::Union{Nothing,AbstractString}, element_dset_size::NTuple{L_dset,Integer}) where L_dset
    N = length(element_dset_size)
    if s == nothing || s == ""
        N == 0 || throw(ErrorException("Array element specified as scalar, but $N unused trailing dimensions in dataset"))
        RealQuantity
    else
        datatype_from_string(s, element_dset_size)
    end
end



Base.@propagate_inbounds _tuple_droplast(x::NTuple{N,Any}, ::Val{M}) where {N,M} =
    Base.ntuple(i -> x[i], Val{N-M}())

_namedtuple_type(members::AbstractVector{<:AbstractString}) = NamedTuple{(Symbol.(members)...,)}


function datatype_from_string(s::AbstractString, dset_size::NTuple{L_dset,Integer}) where L_dset
    N_dset = length(dset_size)
    #@info "datatype_from_string(\"$s\", $dset_size)"

    if s == "real"
        #!!! N_dset == 0 || throw(ErrorException("Array element specified as scalar, but $N_dset unused trailing dimensions in dataset"))
        RealQuantity
    elseif s == "string"
        # dset_size == () || ...
        String
    elseif s == "symbol"
        # dset_size == () || ...
        Symbol
    elseif haskey(_datatype_dict, s)
        _datatype_dict[s]
    else
        m = match(datatype_regexp, s)
        m isa Nothing && throw(ErrorException("Invalid datatype string \"$s\""))
        tp = m[2]
        content = m[6]
        if tp == "struct"
            #!!! @assert dset_size == ()
            _namedtuple_type(split(content, ","))
        elseif tp == "table"
            @assert dset_size == ()
            TypedTables.Table{<:_namedtuple_type(split(content, ","))}
        else
            dims = parse.(Int, split(m[4], ","))
            eltp = content
            if tp == "array_of_equalsized_arrays"
                length(dims) == 2 || throw(ErrorException("Invalid dims $dims for datatype \"$tp\""))
                M = dims[1]; N = dims[2]
                L = M + N
                #!!! L <= L_dset || throw(ErrorException("Expected $M + $N dimensions, but dataset has $L_dset dimensions"))
                # T isa Array || throw(ErrorException("Datatype \"$tp\" with \"array\" nested in \"array_of_equalsized_arrays\" currently not supported"))
                #!!! T = _eldatatype_from_string(eltp, _tuple_droplast(dset_size, Val(L)))
                T = _eldatatype_from_string(eltp, dset_size)
                AbstractArrayOfSimilarArrays{<:T,M,N}
            else
                length(dims) == 1 || throw(ErrorException("Invalid dims $dims for datatype \"$tp\""))
                N = dims[1]
                #!!! N <= N_dset || throw(ErrorException("Expected $dims dimensions, but dataset has $N_dset dimensions"))
                #!!! T = _eldatatype_from_string(eltp, _tuple_droplast(dset_size, Val(N)))
                T = _eldatatype_from_string(eltp, dset_size)
                if tp == "fixedsize_array"
                    N == 1  || throw(ErrorException("Datatype fixedsize_array with $dims dims currently not supported\"$tp\""))
                    # T <: RealQuantity || throw(ErrorException("Element type \"$eltp\" in datatype \"$tp\" currently not supported"))
                    #@info dset_size
                    StaticArray{Tuple{dset_size[1]},<:T,1} where {TPL}
                elseif tp == "array"
                    AbstractArray{<:T,N}
                else
                    throw(ErrorException("Unknown datatype \"$tp\""))
                end
            end
        end
    end
end


function _inner_datatype_to_string(::Type{T}) where T
    s = datatype_to_string(T)
    isempty(s) ? "" : "{$s}"
end


datatype_to_string(::Type{<:RealQuantity}) = "real"

datatype_to_string(T::Type{<:Enum{U}}) where {U} = "enum{"*join(broadcast(x -> "$(string(x))=$(U(x))", instances.(T)), ",")*"}"

datatype_to_string(::Type{<:AbstractString}) = "string"

datatype_to_string(::Type{<:Symbol}) = "symbol"

datatype_to_string(::Type{<:AbstractArray{T,N}}) where {T,N} =
    "array<$N>$(_inner_datatype_to_string(T))"

datatype_to_string(::Type{<:StaticArray{TPL,T,N}}) where {TPL,T<:RealQuantity,N} =
    "fixedsize_array<$N>$(_inner_datatype_to_string(T))"

datatype_to_string(::Type{<:ArrayOfSimilarArrays{T,M,N}}) where {T,M,N} =
    "array_of_equalsized_arrays<$M,$N>$(_inner_datatype_to_string(T))"

datatype_to_string(::Type{<:NamedTuple{K}}) where K = "struct{$(join(K,","))}"

datatype_to_string(::Type{<:TypedTables.Table{<:NamedTuple{K}}}) where K = "table{$(join(K,","))}"



function _eltype(dset::HDF5.HDF5Dataset)
    dtype = HDF5.datatype(dset)
    try
        HDF5.hdf5_to_julia_eltype(dtype)
    finally
        close(dtype)
    end
end


function _cumulative_length(A::VectorOfArrays)
    elem_ptr = ArraysOfArrays.internal_element_ptr(A)
    elem_ptr[(firstindex(elem_ptr) + 1):end] .- first(elem_ptr)
end

function _element_ptrs(clen::Vector{<:Integer})
    vcat([1], clen .+ 1)
end


# _without_units(A::Array, unit::Unitful.Unitlike) = ustrip.(uconvert.(unit, A))
# _with_units(A::Array, unit::Unitful.Unitlike) = A * unit


function hasattribute(
    obj::Union{HDF5.HDF5Dataset, HDF5.HDF5.DataFile}, key::Symbol
) where {T<:Union{AbstractString,Real}}
    key_str = String(key)
    attributes = HDF5.attrs(obj)
    HDF5.exists(attributes, key_str)
end


function getattribute(
    obj::Union{HDF5.HDF5Dataset, HDF5.HDF5.DataFile}, key::Symbol, ::Type{T}
) where {T<:Union{AbstractString,Real,HDF5.HDF5ReferenceObj}}
    key_str = String(key)
    attributes = HDF5.attrs(obj)
    x = read(attributes[key_str])
    x isa T ? x : convert(T, x)
end

function getattribute(
    obj::Union{HDF5.HDF5Dataset, HDF5.HDF5.DataFile}, key::Symbol, default_value::T
) where {T<:Union{AbstractString,Real}}
    if hasattribute(obj, key)
        getattribute(obj, key, T)
    else
        default_value
    end
end


function setattribute!(
    obj::Union{HDF5.HDF5Dataset, HDF5.HDF5Group}, key::Symbol,
    value::Union{AbstractString,Real,HDF5.HDF5ReferenceObj}
)
    HDF5.attrs(obj)[String(key)] = value
end



h5ref(ds::HDF5.HDF5Dataset) = HDF5.HDF5ReferenceObj(ds.file, HDF5.name(ds))

h5deref(ref::HDF5.HDF5ReferenceObj, context::Union{HDF5.HDF5Dataset, HDF5.HDF5Group}) = context.file[ref]


LegendDataTypes.getunits(dset::HDF5.HDF5Dataset) = units_from_string(getattribute(dset, :units, ""))

LegendDataTypes.setunits!(dset::HDF5.HDF5Dataset, units::Unitful.Unitlike) = setattribute!(dset, :units, units_to_string(units))


default_datatype(dset::HDF5.HDF5Dataset) = AbstractArray{<:RealQuantity,length(size(dset))}
default_datatype(df::HDF5.DataFile) = NamedTuple{(Symbol.(names(df))...,)}

_size(dset::HDF5.HDF5Dataset) = size(dset)
_size(df::HDF5.DataFile) = ()

function getdatatype(input::Union{HDF5.HDF5Dataset, HDF5.DataFile})
    dtstr = getattribute(input, :datatype, "")
    sz = _size(input)
    isempty(dtstr) ? default_datatype(input) : datatype_from_string(dtstr, sz)
end

function setdatatype!(output::Union{HDF5.HDF5Dataset, HDF5.DataFile}, datatype::Type)
    setattribute!(output, :datatype, datatype_to_string(datatype))
end



function getinfo(dset::HDF5.HDF5Dataset)
    (
        datatype = getattribute(dset, :datatype, ""),
        units = getunits(dset),
        cumsum_length_ds = getattribute(dset, :cumsum_length_ds, ""),
    )
end



function LegendDataTypes.readdata(input::Union{HDF5.HDF5Dataset, HDF5.DataFile}, name::AbstractString)
    datatype = getdatatype(input[name])
    readdata(input, name, datatype)
end


function LegendDataTypes.writedata(
    output::HDF5.DataFile, name::AbstractString,
    x::Union{T,AbstractArray{T}},
    fulldatatype::DataType = typeof(x)
) where {T<:RealQuantity}
    units = unit(eltype(x))
    # @debug "name" name
    if units == NoUnits
        # @info "without units"
        #output[name, "shuffle", (), "deflate", 3] = x
        output[name] = x
    else
        # @debug "with units"
        output[name] = ustrip.(x)
        setunits!(output[name], units)
    end
    if fulldatatype != Nothing
        setdatatype!(output[name], fulldatatype)
    end
    nothing
end

function LegendDataTypes.readdata(
    input::HDF5.DataFile, name::AbstractString,
    ::Type{<:Union{RealQuantity,AbstractArray}}
)
    dset = input[name]
    data = read(dset)#::AT
    units = getunits(dset)
    if units == NoUnits
        data
    else
        data * units
    end
end


function LegendDataTypes.writedata(
    output::HDF5.DataFile, name::AbstractString,
    x::AbstractArray{<:Enum{T}},
    fulldatatype::DataType = typeof(x)
) where {T}
    writedata(output, name, reinterpret(T, x), fulldatatype)
end

function LegendDataTypes.readdata(
    input::HDF5.DataFile, name::AbstractString,
    AT::Type{<:AbstractArray{<:Enum}}
)
    ET = AT.body.parameters[1].ub
    data = readdata(input, name, AbstractArray{RealQuantity})
    ET.(data)
end


function LegendDataTypes.writedata(
    output::HDF5.DataFile, name::AbstractString,
    x::AbstractArray{T,N},
    fulldatatype::DataType = typeof(x)
) where {TPL,T<:StaticArray{TPL,<:RealQuantity},N}
    writedata(output, name, flatview(x), fulldatatype)
end


function LegendDataTypes.readdata(
    input::HDF5.DataFile, name::AbstractString,
    AT::Type{<:AbstractArray{<:StaticArray{Tuple{N}}}}
) where {N}
    data = readdata(input, name, AbstractArray{RealQuantity,2})
    nestedview(data, SVector{N})
end


function LegendDataTypes.writedata(
    output::HDF5.DataFile, name::AbstractString,
    x::AbstractString,
    fulldatatype::DataType = typeof(x)
)
    output[name] = x
    if fulldatatype != Nothing
        setdatatype!(output[name], fulldatatype)
    end
    nothing
end

function LegendDataTypes.readdata(
    input::HDF5.DataFile, name::AbstractString,
    ::Type{<:String}
)
    dset = input[name]
    read(dset)#::AT
end


function LegendDataTypes.writedata(
    output::HDF5.DataFile, name::AbstractString,
    x::Symbol,
    fulldatatype::DataType = typeof(x)
)
    writedata(output, name, String(x), fulldatatype)
end

function LegendDataTypes.readdata(
    input::HDF5.DataFile, name::AbstractString,
    ::Type{<:Symbol}
)
    Symbol(readdata(input, name, String))
end


function LegendDataTypes.readdata(
    input::HDF5.DataFile, name::AbstractString,
    AT::Type{<:AbstractArray{<:AbstractArray}}
)
    data = readdata(input, name, AbstractArray)
    dset = input[name]
    clen = read(h5deref(getattribute(dset, :cumsum_length, HDF5.HDF5ReferenceObj), dset))
    data_vec = VectorOfVectors(data, _element_ptrs(clen))

    if hasattribute(dset, :codec)
        codec_name = Symbol(getattribute(dset, :codec, String))
        C = LegendDataTypes.array_codecs[codec_name]
        codec = read_from_properties(getattribute, dset, C)
        n = length(data_vec)
        codec_vec = StructArray(fill(codec, n)) # ToDo: Improve
        decoded_length = read(h5deref(getattribute(dset, :decoded_length, HDF5.HDF5ReferenceObj), dset))
        size_vec = broadcast(x -> (Int(x),), decoded_length)

        StructArray{EncodedArray{Int32,1,C,Vector{UInt8}}}((
            codec_vec,
            size_vec,
            data_vec
        ))
    else
        data_vec
    end
end

# Hack:
function LegendDataTypes.readdata(
    input::HDF5.DataFile, name::AbstractString,
    AT::Type{<:AbstractArray{<:AbstractArray{<:StaticVector{3}}}}
)
    data = readdata(input, name, AbstractArray)
    nested_data = copy(nestedview(data, SVector{3}))
    dset = input[name]
    clen = read(h5deref(getattribute(dset, :cumsum_length, HDF5.HDF5ReferenceObj), dset))
    VectorOfVectors(nested_data, _element_ptrs(clen))
end

function LegendDataTypes.writedata(
    output::HDF5.DataFile, name::AbstractString,
    x::VectorOfArrays{T,1},
    fulldatatype::DataType = typeof(x)
) where {T,N}
    writedata(output, name, flatview(x), fulldatatype)
    dset = output[name]
    output["$(name)_clen"] = _cumulative_length(x)
    setattribute!(dset, :cumsum_length, h5ref(output["$(name)_clen"]))
    nothing
end

function LegendDataTypes.writedata(
    output::HDF5.DataFile, name::AbstractString,
    x::VectorOfEncodedArrays{T,1,C},
    fulldatatype::DataType = typeof(x)
) where {T,N,C}
    codec = first(x.codec)
    @inbounds for c in x.codec
        c != codec && throw("Can't write VectorOfEncodedArrays that has non-uniform codec parameters")
    end

    writedata(output, name, x.encoded, fulldatatype)
    dset = output[name]
    decoded_length = Int32.(prod.(x.size))
    output["$(name)_declen"] = decoded_length
    decsize_vec_ds = output["$(name)_declen"]
    setattribute!(dset, :decoded_length, h5ref(decsize_vec_ds))
 
    codec_name = LegendDataTypes.array_codecs[C]
    setattribute!(dset, :codec, String(codec_name))
    write_to_properties!(setattribute!, dset, codec)

    nothing
end


function LegendDataTypes.writedata(
    output::HDF5.DataFile, name::AbstractString,
    x::AbstractArrayOfSimilarArrays{T,M,N},
    fulldatatype::DataType = typeof(x)
) where {T,M,N}
    writedata(output, name, flatview(x), fulldatatype)
    nothing
end

# Hack:
function LegendDataTypes.readdata(
    input::HDF5.DataFile, name::AbstractString,
    AT::Type{<:AbstractVectorOfSimilarVectors}
)
    nestedview(readdata(input, name, AbstractArray{<:RealQuantity}))
end


function LegendDataTypes.writedata(
    output::HDF5.DataFile, name::AbstractString,
    x::NamedTuple,
    fulldatatype::DataType = typeof(x)
)
    for k in keys(x)
        y = x[k]
        writedata(output, "$name/$k", y, typeof(y))
    end
    setdatatype!(output[name], fulldatatype)
    nothing
end

function LegendDataTypes.readdata(
    input::HDF5.DataFile, name::AbstractString,
    AT::Type{<:NamedTuple}
)
    syms = AT.body.parameters[1]
    NamedTuple{syms}(map(k -> readdata(input, "$name/$k"), syms))
end


function LegendDataTypes.writedata(
    output::HDF5.DataFile, name::AbstractString,
    x::Any,
    fulldatatype::DataType = typeof(x)
)
    Tables.istable(x) || throw(ArgumentError("Value to write, of type $(typeof(x)),is not a table"))
    cols = Tables.columns(x)
    writedata(output, name, cols, fulldatatype)
end

function LegendDataTypes.readdata(
    input::HDF5.DataFile, name::AbstractString,
    AT::Type{<:TypedTables.Table}
)
    # Hack:
    TypedTables.Table(readdata(input, name, NamedTuple{AT.var.ub.body.parameters[1]}))
end
