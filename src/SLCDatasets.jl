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

Sentinel-1 IW SLCs are read from a `.SAFE` directory or the zip of one, taking geometry from the
annotation XML and state vectors from a POEORB or RESORB `.EOF` file, which a Sentinel-1 product does
not carry. Only the annotation is read, so the measurement TIFFs cost nothing. A TOPS product is three
subswaths of bursts rather than one image; [`open_slc`](@ref) describes either the mosaic across them or
one individual burst.
"""
module SLCDatasets

import Dates
using Dates: DateTime
using EzXML: parsexml, readxml, root, findfirst, findall, nodecontent, eachelement
using HDF5: h5open, ishdf5, read_attribute
using Mmap: mmap
using StaticArrays: SVector
using ZipArchives: ZipReader, zip_name, zip_nentries, zip_readentry

export open_slc, orbit, nlines, nsamples, start_datetime, stop_datetime
export repeat_interval, epoch_offset, nbursts
export SLC, Identification, RadarGeometry, StateVectors
export LocalFile, RemoteHTTP, RemoteS3

# `LookSide`, `LookLeft` and `LookRight` are deliberately not exported: a geometry package consuming
# this one defines its own, and exporting both makes the name ambiguous at every call site.

include("types.jl")
include("nisar.jl")
include("sentinel1.jl")
# `remote.jl` before `source.jl`: the remote sources are what `open_slc` dispatches on.
include("remote.jl")
include("source.jl")
include("pairing.jl")

end
