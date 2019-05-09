# This file is a part of LegendHDF5IO.jl, licensed under the MIT License (MIT).

__precompile__(true)

module LegendHDF5IO

using ArraysOfArrays
using ElasticArrays
using LegendDataTypes
using RadiationDetectorSignals
using StaticArrays
using StructArrays
using Tables
using Unitful

import HDF5
import TypedTables

using LegendDataTypes: readdata, writedata, getunits, setunits!, units_from_string, units_to_string

using RadiationDetectorSignals: RealQuantity, ArrayOfDims, AosAOfDims, SArrayOfDims,
    recursive_ndims


include("generic_io.jl")
include("radsig_io.jl")


function __init__()
    _datatype_dict[LegendHDF5IO.datatype_to_string(EventType)] = EventType
end

end # module
