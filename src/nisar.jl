# The NISAR backend: HDF5 products under `/science/{LSAR,SSAR}/<productType>`.
#
# Group paths are built the way isce3's `nisar.products.readers.Base` builds them — a root path found
# by probing for the sensor band, then the product type, then a fixed subgroup name. The product type
# is read from the identification group and selects the reader, which is also how isce3's
# `open_product` dispatches. A product whose type is `SLC` is an early-mission RSLC.

const SCIENCE_ROOT = "science"
const SENSOR_BANDS = ("LSAR", "SSAR")

# Speed of light, for converting the processed center frequency to a wavelength. This is the value
# isce3 uses (`isce3.core.speed_of_light`), and the conversion must match it to the bit.
const SPEED_OF_LIGHT = 299792458.0

"""
    NisarBackend <: AbstractSARBackend

A NISAR-format HDF5 product: a file, the sensor band group inside it, and the product type.

`frequency` selects the sub-band (`"A"` or `"B"`); a product lists the ones it carries at
`identification/listOfFrequencies`.
"""
struct NisarBackend <: AbstractSARBackend
    path::String
    band::String
    product_type::String
    frequency::String
end

product_path(b::NisarBackend) = string(SCIENCE_ROOT, "/", b.band, "/", b.product_type)
metadata_path(b::NisarBackend) = string(product_path(b), "/metadata")
swath_path(b::NisarBackend) = string(product_path(b), "/swaths")
frequency_path(b::NisarBackend) = string(swath_path(b), "/frequency", b.frequency)
identification_path(b::NisarBackend) = string(SCIENCE_ROOT, "/", b.band, "/identification")
orbit_path(b::NisarBackend) = string(metadata_path(b), "/orbit")

# A geocoded product stores its samples on a map grid under `grids` rather than `swaths`. Only the
# radar-geometry products carry the slant-range/azimuth axes this package reports, so the geocoded
# ones are recognized in order to be rejected with a clear message rather than a missing-group error.
const GEOCODED_TYPES = ("GSLC", "GCOV", "GUNW", "GOFF")

"""
    nisar_band(path) -> String

The sensor band group of a NISAR product: `"LSAR"` or `"SSAR"`.
"""
function nisar_band(path::AbstractString)
    return h5open(path, "r") do h
        haskey(h, SCIENCE_ROOT) || throw(ArgumentError(
            "`$path` has no `/$SCIENCE_ROOT` group, so it is not a NISAR-format product"))
        for band in SENSOR_BANDS
            haskey(h[SCIENCE_ROOT], band) && return band
        end
        throw(ArgumentError("`$path` has none of $(join(SENSOR_BANDS, ", ")) under " *
                            "`/$SCIENCE_ROOT`, so it is not a NISAR-format product"))
    end
end

"""
    nisar_product_type(path, band) -> String

The product type recorded in a NISAR product's identification group, as the group name that holds
it. Early-mission products name the group `SLC` where later ones name it `RSLC`; both report
`RSLC`.
"""
function nisar_product_type(path::AbstractString, band::AbstractString)
    return h5open(path, "r") do h
        g = h[string(SCIENCE_ROOT, "/", band)]
        haskey(g, "identification") || throw(ArgumentError(
            "`$path` has no `identification` group under `/$SCIENCE_ROOT/$band`"))
        declared = _string(read(g["identification/productType"]))
        # The declared type names the group except for the early-mission spelling.
        declared == "SLC" && return haskey(g, "RSLC") ? "RSLC" : "SLC"
        return declared
    end
end

# Every scalar string in a NISAR product is stored as fixed-length bytes, and HDF5.jl surfaces those
# variously as a `String` or as a byte vector depending on the dataset. Trailing NULs are padding.
_string(x::AbstractString) = String(rstrip(x, '\0'))
_string(x::AbstractVector{UInt8}) = _string(String(x))
_string(x::AbstractArray) = _string(only(x))

"""
    parse_cf_epoch(units) -> DateTime

The reference instant of a CF-convention units string such as
`"seconds since 2025-10-28T00:00:00"`.

Only seconds are accepted: the products this reads record times in seconds, and silently rescaling a
different unit would corrupt every time by a constant factor.
"""
function parse_cf_epoch(units::AbstractString)
    s = _string(units)
    m = match(r"^\s*(\w+)\s+since\s+(.+?)\s*$", s)
    m === nothing && throw(ArgumentError(
        "cannot read a reference epoch from the units string \"$s\"; expected " *
        "\"<unit> since <timestamp>\""))
    unit = lowercase(m[1])
    unit in ("second", "seconds") || throw(ArgumentError(
        "time units are \"$unit\" in \"$s\", but only seconds are supported"))
    stamp = replace(m[2], " " => "T")
    # A trailing zone designator is dropped: these products are UTC, and `DateTime` carries no zone.
    stamp = replace(stamp, r"(Z|[+-]\d{2}:?\d{2})$" => "")
    return DateTime(stamp)
end

# Products record `zeroDopplerStartTime` to nanoseconds, which is finer than `DateTime`'s
# milliseconds. `Identification` keeps the string so nothing is lost; this converts only where
# millisecond resolution is enough, dropping the extra digits rather than rounding into them.
function _truncated_datetime(s::AbstractString)
    t = _string(s)
    t = replace(t, r"(Z|[+-]\d{2}:?\d{2})$" => "")
    m = match(r"^(.*\.\d{1,3})\d*$", t)
    return DateTime(m === nothing ? t : m[1])
end

function read_identification(b::NisarBackend)
    return h5open(b.path, "r") do h
        g = h[identification_path(b)]
        get_str(name, default = "") = haskey(g, name) ? _string(read(g[name])) : default
        return Identification(
            get_str("missionId"),
            get_str("productType"),
            Int(only(read(g["absoluteOrbitNumber"]))),
            get_str("orbitPassDirection"),
            get_str("lookDirection"),
            get_str("zeroDopplerStartTime"),
            get_str("zeroDopplerEndTime"),
            get_str("boundingPolygon"),
        )
    end
end

function read_geometry(b::NisarBackend)
    b.product_type in GEOCODED_TYPES && throw(ArgumentError(
        "`$(b.path)` is a $(b.product_type) product, which stores its samples on a map grid and " *
        "carries no slant-range/azimuth geometry; open a radar-geometry product (RSLC) instead"))
    return h5open(b.path, "r") do h
        freq = h[frequency_path(b)]
        swaths = h[swath_path(b)]
        ident = h[identification_path(b)]

        slant_range = read(freq["slantRange"])
        zd_time = read(swaths["zeroDopplerTime"])
        zd_spacing = only(read(swaths["zeroDopplerTimeSpacing"]))
        center_frequency = only(read(freq["processedCenterFrequency"]))

        look = _string(read(ident["lookDirection"]))
        side = lowercase(look) == "left" ? LookLeft :
               lowercase(look) == "right" ? LookRight :
               throw(ArgumentError(
                   "`$(b.path)` records lookDirection \"$look\"; expected \"Left\" or \"Right\""))

        # The azimuth times and the orbit's state vector times share one epoch in this format, so the
        # epoch is read once here and reported for both.
        epoch = parse_cf_epoch(read_attribute(swaths["zeroDopplerTime"], "units"))

        return RadarGeometry(
            first(slant_range),
            last(slant_range),
            only(read(freq["slantRangeSpacing"])),
            SPEED_OF_LIGHT / center_frequency,
            1 / zd_spacing,
            first(zd_time),
            last(zd_time),
            length(zd_time),
            length(slant_range),
            side,
            epoch,
        )
    end
end

function read_orbit(b::NisarBackend)
    return h5open(b.path, "r") do h
        g = h[orbit_path(b)]
        time = Vector{Float64}(read(g["time"]))
        # Stored (N, 3) row-major, which HDF5.jl presents as (3, N).
        pos = read(g["position"])
        vel = read(g["velocity"])
        size(pos, 1) == 3 || throw(ArgumentError(
            "orbit position in `$(b.path)` has leading dimension $(size(pos, 1)), expected 3"))
        axes(pos) == axes(vel) || throw(DimensionMismatch(
            "orbit position and velocity in `$(b.path)` have axes $(axes(pos)) and $(axes(vel))"))
        size(pos, 2) == length(time) || throw(DimensionMismatch(
            "orbit in `$(b.path)` has $(length(time)) times but $(size(pos, 2)) state vectors"))
        epoch = parse_cf_epoch(read_attribute(g["time"], "units"))
        get_str(name, default) = haskey(g, name) ? _string(read(g[name])) : default
        return StateVectors(
            time,
            [SVector{3,Float64}(pos[1, i], pos[2, i], pos[3, i]) for i in axes(pos, 2)],
            [SVector{3,Float64}(vel[1, i], vel[2, i], vel[3, i]) for i in axes(vel, 2)],
            epoch,
            get_str("interpMethod", "Hermite"),
            get_str("orbitType", "Custom"),
        )
    end
end
