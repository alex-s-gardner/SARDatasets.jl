"""
    SAR

Read SAR acquisitions into one type, whatever the sensor.

[`open_sar`](@ref) returns a [`Radar`](@ref) carrying an [`Identification`](@ref) and a
[`RadarGeometry`](@ref) — the slant-range/azimuth geometry in the units a geometry kernel wants — plus
an [`Orbit`](@ref), read on first access.

NISAR-format HDF5 products are read natively. The metadata of one is a few tens of kilobytes near the
start of a file tens of gigabytes long, so a remote product can be opened without transferring it.
"""
module SAR

using Dates: DateTime
using HDF5: h5open, ishdf5, read_attribute
using StaticArrays: SVector

export open_sar, orbit, nlines, nsamples, start_datetime, stop_datetime
export Radar, Identification, RadarGeometry, Orbit
export LookSide, LookLeft, LookRight
export LocalFile

include("types.jl")
include("nisar.jl")
include("source.jl")

end
