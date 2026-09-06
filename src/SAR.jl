"""
    SAR

Read SAR acquisitions into one type, whatever the sensor.

[`open_sar`](@ref) returns a [`Radar`](@ref) carrying an [`Identification`](@ref) and a
[`RadarGeometry`](@ref) — the slant-range/azimuth geometry in the units a geometry kernel wants — plus
a [`StateVectors`](@ref) record, read on first access.

NISAR-format HDF5 products are read natively. The metadata of one is a few tens of kilobytes near the
start of a file tens of gigabytes long, so a remote product can be opened without transferring it.
"""
module SAR

import Dates
using Dates: DateTime
using HDF5: h5open, ishdf5, read_attribute
using StaticArrays: SVector

export open_sar, orbit, nlines, nsamples, start_datetime, stop_datetime
export radar_coordinate, image_pair, repeat_interval
export Radar, Identification, RadarGeometry, StateVectors
export LocalFile, RemoteHTTP, RemoteS3

# `LookSide`, `LookLeft` and `LookRight` are deliberately not exported: a geometry package consuming
# this one defines its own, and exporting both makes the name ambiguous at every call site.

include("types.jl")
include("nisar.jl")
# `remote.jl` before `source.jl`: the remote sources are what `open_sar` dispatches on.
include("remote.jl")
include("source.jl")
include("geometry_api.jl")

end
