# The single entry point that opens a product, and the dispatch that picks its reader.
#
# A source is separate from a backend: the backend knows a sensor's group and dataset names, the
# source knows how to get a readable file. Splitting them means a new sensor costs a backend and a new
# access route costs a source, rather than one costing both.

"""
    open_slc(src; frequency = nothing) -> SLC

Open a SAR acquisition.

`src` is an [`AbstractSLCSource`](@ref), or a path or URL, which is resolved to one. The sensor
backend is chosen by inspecting the product rather than by the caller naming it, so the same call
serves every supported format.

`frequency` selects the sub-band for a sensor that has them, defaulting to the first the product
lists.

The identification and geometry records are read here; the state vectors are read on first access.

# Examples

```julia
s = open_slc("NISAR_L1_PR_RSLC_....h5")
s.geometry.prf
```
"""
function open_slc(src::AbstractSLCSource; frequency = nothing)
    path = localpath(src)
    isfile(path) || throw(ArgumentError("`$path` is not a readable file"))
    ishdf5(path) || throw(ArgumentError(
        "`$path` is not an HDF5 file; the only formats currently read are NISAR-format HDF5 " *
        "products"))
    band = nisar_band(path)
    product_type = nisar_product_type(path, band)
    freq = frequency === nothing ? default_frequency(path, band) : String(frequency)
    backend = NisarBackend(path, band, product_type, freq)
    return SLC(backend, read_identification(backend), read_geometry(backend))
end

# A partly fetched product fails inside HDF5, reading zeros out of the hole past the prefetch window.
# Reporting that as-is would leave a caller staring at an HDF5 stack trace, so the two remote sources
# translate it into the knob that fixes it.
function open_slc(src::Union{RemoteHTTP,RemoteS3}; frequency = nothing)
    path = localpath(src)
    return try
        open_slc(LocalFile(path); frequency)
    catch e
        e isa ArgumentError && rethrow()
        throw(_prefetch_error(src, path, e))
    end
end

open_slc(path::AbstractString; kwargs...) = open_slc(source_for(path); kwargs...)

"""
    source_for(spec) -> AbstractSLCSource

The source a path, URL or S3 URI denotes.
"""
function source_for(spec::AbstractString)
    startswith(spec, "s3://") && return RemoteS3(spec)
    (startswith(spec, "http://") || startswith(spec, "https://")) && return RemoteHTTP(spec)
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
