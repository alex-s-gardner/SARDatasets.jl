# The single entry point that opens a product, and the dispatch that picks its reader.
#
# A source is separate from a backend: the backend knows a sensor's group and dataset names, the
# source knows how to get a readable file. Splitting them means a new sensor costs a backend and a new
# access route costs a source, rather than one costing both.

"""
    open_slc(src; frequency = nothing, orbit = nothing, swath = nothing, swaths = nothing,
             burst = nothing, polarization = nothing) -> SLC

Open a SAR acquisition.

`src` is an [`AbstractSLCSource`](@ref), or a path or URL, which is resolved to one. The sensor
backend is chosen by inspecting the product rather than by the caller naming it, so the same call
serves every supported format.

The identification and geometry records are read here; the state vectors are read on first access.

# NISAR

`frequency` selects the sub-band, defaulting to the first the product lists.

# Sentinel-1

`orbit` is the path to the POEORB or RESORB `.EOF` file holding the state vectors, and is required:
a Sentinel-1 product does not contain its own orbit.

A TOPS product is three subswaths of bursts rather than one image, so there are two geometries to
ask for. By default the result spans `swaths` — all three — the way hyp3-autorift's
`loadMetadataSlc` does. Passing `burst` instead describes that single burst of the single subswath
`swath`, matching `loadMetadata`. Use [`nbursts`](@ref) for the burst count.

`polarization` defaults to the first of `vv`, `vh`, `hh`, `hv` the product carries.

# Examples

```julia
s = open_slc("NISAR_L1_PR_RSLC_....h5")
s.geometry.prf

mosaic = open_slc("S1A_IW_SLC_....zip"; orbit = "S1A_OPER_AUX_POEORB_....EOF")
one    = open_slc("S1A_IW_SLC_....zip"; orbit = "...EOF", swath = 2, burst = 3)
```
"""
function open_slc(src::AbstractSLCSource; frequency = nothing, orbit = nothing,
                  swath = nothing, swaths = nothing, burst = nothing, polarization = nothing)
    path = localpath(src)
    ispath(path) || throw(ArgumentError("`$path` is not a readable file"))
    if is_safe_product(path)
        return _open_sentinel1(path; orbit, swath, swaths, burst, polarization)
    end
    isfile(path) || throw(ArgumentError("`$path` is not a readable file"))
    ishdf5(path) || throw(ArgumentError(
        "`$path` is neither an HDF5 file nor a Sentinel-1 SAFE product; those are the only " *
        "formats currently read"))
    band = nisar_band(path)
    product_type = nisar_product_type(path, band)
    freq = frequency === nothing ? default_frequency(path, band) : String(frequency)
    backend = NisarBackend(path, band, product_type, freq)
    return SLC(backend, read_identification(backend), read_geometry(backend))
end

function _open_sentinel1(path::AbstractString; orbit, swath, swaths, burst, polarization)
    orbit === nothing && throw(ArgumentError(
        "`$path` is a Sentinel-1 product, whose state vectors live in a separate POEORB/RESORB " *
        "`.EOF` file rather than in the product. Pass `orbit = \"<file>.EOF\"`"))
    orbit_path = String(orbit)
    isfile(orbit_path) || throw(ArgumentError("`$orbit_path` is not a readable orbit file"))

    swath !== nothing && swaths !== nothing && throw(ArgumentError(
        "pass either `swath` for one subswath or `swaths` for several, not both"))
    # A single burst belongs to one subswath, so naming a burst without a subswath would silently
    # pick one; requiring `swath` alongside `burst` keeps the choice the caller's.
    burst !== nothing && swath === nothing && swaths === nothing && throw(ArgumentError(
        "`burst` selects a burst of one subswath, so it needs `swath = 1`, `2` or `3` as well"))

    selected = swath !== nothing ? [Int(swath)] :
               swaths !== nothing ? sort(collect(Int, swaths)) : [1, 2, 3]
    isempty(selected) && throw(ArgumentError("`swaths` is empty; name at least one subswath"))
    all(s -> 1 <= s <= 3, selected) || throw(ArgumentError(
        "Sentinel-1 IW subswaths are 1, 2 and 3; got $(join(selected, ", "))"))
    allunique(selected) || throw(ArgumentError(
        "`swaths` repeats a subswath: $(join(selected, ", "))"))
    burst !== nothing && length(selected) == 1 || burst === nothing || throw(ArgumentError(
        "a burst belongs to one subswath, but $(length(selected)) were named"))

    pol = polarization === nothing ? default_polarization(path) : lowercase(String(polarization))
    available = safe_polarizations(path)
    pol in available || throw(ArgumentError(
        "`$path` carries polarization$(length(available) == 1 ? " " : "s ")" *
        "$(join(uppercase.(available), ", ")), not $(uppercase(pol))"))

    backend = Sentinel1Backend(path, orbit_path, pol, selected,
                               burst === nothing ? nothing : Int(burst))
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
