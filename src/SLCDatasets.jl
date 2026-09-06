"""
    SLCDatasets

Read single-look complex SAR products into one type, whatever the sensor.

[`open_slc`](@ref) returns an [`SLC`](@ref) carrying an [`Identification`](@ref) and a
[`RadarGeometry`](@ref) — the slant-range/azimuth geometry in the units a geometry kernel wants — plus
a [`StateVectors`](@ref) record, read on first access.

The scope is SLCs in radar geometry, not SAR products generally. An interferogram or a covariance
product carries a multilooked grid, and a geocoded SLC carries map coordinates rather than a slant-range
axis — neither is described by [`RadarGeometry`](@ref), so neither is read here.

NISAR-format HDF5 range-Doppler products (RSLC) are read natively. The metadata of one is a few tens of
kilobytes near the start of a file tens of gigabytes long, so a remote product can be opened without
transferring it.
"""
module SLCDatasets

import Dates
using Dates: DateTime
using HDF5: h5open, ishdf5, read_attribute
using StaticArrays: SVector

export open_slc, orbit, nlines, nsamples, start_datetime, stop_datetime
export repeat_interval, epoch_offset
export SLC, Identification, RadarGeometry, StateVectors
export LocalFile, RemoteHTTP, RemoteS3

# `LookSide`, `LookLeft` and `LookRight` are deliberately not exported: a geometry package consuming
# this one defines its own, and exporting both makes the name ambiguous at every call site.

include("types.jl")
include("nisar.jl")
# `remote.jl` before `source.jl`: the remote sources are what `open_slc` dispatches on.
include("remote.jl")
include("source.jl")
include("pairing.jl")

end
