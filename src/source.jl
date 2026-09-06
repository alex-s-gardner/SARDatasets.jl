# Where a product's bytes come from, and the single entry point that opens one.
#
# A source is separate from a backend: the backend knows a sensor's group and dataset names, the
# source knows how to get a readable file. Splitting them means a new sensor costs a backend and a new
# access route costs a source, rather than one costing both.

"""
    AbstractSARSource

Where a product's bytes come from. [`LocalFile`](@ref) is the only route that needs no network.
"""
abstract type AbstractSARSource end

"""
    LocalFile(path)

A product already on disk.
"""
struct LocalFile <: AbstractSARSource
    path::String
end

LocalFile(path::AbstractString) = LocalFile(String(path))

"""
    localpath(src::AbstractSARSource) -> String

A path that can be opened, materializing whatever the source needs to make that true.
"""
localpath(src::LocalFile) = src.path

"""
    open_sar(src; frequency = nothing) -> Radar

Open a SAR acquisition.

`src` is an [`AbstractSARSource`](@ref), or a path or URL, which is resolved to one. The sensor
backend is chosen by inspecting the product rather than by the caller naming it, so the same call
serves every supported format.

`frequency` selects the sub-band for a sensor that has them, defaulting to the first the product
lists.

The identification and geometry records are read here; the orbit and the image bands are read on
first access.

# Examples

```julia
s = open_sar("NISAR_L1_PR_RSLC_....h5")
s.geometry.prf
```
"""
function open_sar(src::AbstractSARSource; frequency = nothing)
    path = localpath(src)
    isfile(path) || throw(ArgumentError("`$path` is not a readable file"))
    ishdf5(path) || throw(ArgumentError(
        "`$path` is not an HDF5 file; the only formats currently read are NISAR-format HDF5 " *
        "products"))
    band = nisar_band(path)
    product_type = nisar_product_type(path, band)
    freq = frequency === nothing ? default_frequency(path, band) : String(frequency)
    backend = NisarBackend(path, band, product_type, freq)
    return Radar(backend, read_identification(backend), read_geometry(backend))
end

open_sar(path::AbstractString; kwargs...) = open_sar(source_for(path); kwargs...)

"""
    source_for(spec) -> AbstractSARSource

The source a path or URL denotes.
"""
function source_for(spec::AbstractString)
    return LocalFile(spec)
end

"""
    default_frequency(path, band) -> String

The first sub-band a NISAR product lists, as the one to read when the caller names none.
"""
function default_frequency(path::AbstractString, band::AbstractString)
    return h5open(path, "r") do h
        name = string(SCIENCE_ROOT, "/", band, "/identification/listOfFrequencies")
        haskey(h, name) || return "A"
        listed = read(h[name])
        isempty(listed) && throw(ArgumentError("`$path` lists no frequencies"))
        return _string(first(listed))
    end
end
