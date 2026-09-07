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
    NisarBackend <: AbstractSLCBackend

A NISAR-format HDF5 product: a file, the sensor band group inside it, and the product type.

`frequency` selects the sub-band (`"A"` or `"B"`); a product lists the ones it carries at
`identification/listOfFrequencies`.
"""
struct NisarBackend <: AbstractSLCBackend
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

# A geocoded product stores its samples on a map grid under `grids` rather than `swaths`, so it has no
# slant-range axis for `RadarGeometry` to describe and is out of this package's scope. Recognized by name
# in order to say that, rather than failing on a missing group.
const GEOCODED_TYPES = ("GSLC", "GCOV", "GUNW", "GOFF")

"""
    nisar_band(h) -> String

The sensor band group of an open NISAR product: `"LSAR"` or `"SSAR"`.
"""
function nisar_band(h)::String
    haskey(h, SCIENCE_ROOT) || throw(ArgumentError(
        "`$(HDF5.filename(h))` has no `/$SCIENCE_ROOT` group, so it is not a NISAR-format product"))
    science = _group(h, SCIENCE_ROOT)
    for band in SENSOR_BANDS
        haskey(science, band) && return band
    end
    throw(ArgumentError("`$(HDF5.filename(h))` has none of $(join(SENSOR_BANDS, ", ")) under " *
                        "`/$SCIENCE_ROOT`, so it is not a NISAR-format product"))
end

nisar_band(path::AbstractString) = h5open(nisar_band, path, "r")

"""
    nisar_product_type(h, band) -> String

The product type recorded in an open NISAR product's identification group, as the group name that
holds it. Early-mission products name the group `SLC` where later ones name it `RSLC`; both report
`RSLC`.
"""
function nisar_product_type(h, band::AbstractString)::String
    g = _group(h, string(SCIENCE_ROOT, "/", band))
    haskey(g, "identification") || throw(ArgumentError(
        "`$(HDF5.filename(h))` has no `identification` group under `/$SCIENCE_ROOT/$band`"))
    declared = _text(g, "identification/productType")
    # The declared type names the group except for the early-mission spelling.
    declared == "SLC" && return haskey(g, "RSLC") ? "RSLC" : "SLC"
    return declared
end

nisar_product_type(path::AbstractString, band::AbstractString) =
    h5open(h -> nisar_product_type(h, band), path, "r")

# Indexing an HDF5 file or group gives back `Union{Attribute,Dataset,Datatype,Group}`, so `read` on the
# result infers as `Any` and every value derived from it dispatches at runtime. These accessors assert
# what the layout says a name is, which both keeps inference concrete and turns a product whose
# structure differs from the NISAR layout into an error naming the path rather than a `MethodError`
# deeper in.

@noinline _not_a(kind, parent, name) = throw(ArgumentError(
    "`$name` under `$(HDF5.name(parent))` is a $(nameof(typeof(parent))) member that is not a " *
    "$kind, so this is not a NISAR-format product"))

function _dataset(parent, name::AbstractString)
    d = parent[name]
    d isa HDF5.Dataset || _not_a("dataset", parent, name)
    return d
end

function _group(parent, name::AbstractString)
    g = parent[name]
    g isa HDF5.Group || _not_a("group", parent, name)
    return g
end

# `read` on a dataset infers as `Any`, because a dataset's element type is only known once the file is
# open. So each accessor below converts to the type the layout calls for and annotates that, which keeps
# inference concrete through the arithmetic that follows. The conversion is explicit rather than
# `read(d, T)`, which reinterprets the stored bytes as `T` and errors on a width mismatch instead of
# converting — a product storing a count as `Int32` or times as `Float32` is read correctly here.

# A dataset's value where only a scalar is wanted.
_scalar(parent, name::AbstractString)::Float64 = _only_number(_dataset(parent, name))

# A whole-number field, stored variously as an integer or a float across product generations.
function _scalar_int(parent, name::AbstractString)
    d = _dataset(parent, name)
    v = _only_number(d)
    isinteger(v) || throw(ArgumentError(
        "`$(HDF5.name(d))` holds $v where a whole number is expected, so this is not a " *
        "NISAR-format product"))
    return Int(v)
end

# A scalar dataset holds a number, a zero-dimensional array, or a one-element one depending on how the
# product was written; `only` covers the arrays and the bare number falls through.
_only_number(d)::Float64 = _as_number(read(d))
_as_number(x::Number) = Float64(x)
_as_number(x::AbstractArray) = Float64(only(x))

# A fixed-length string dataset, surfaced as a `String` or as bytes; `_string` narrows either and strips
# the NUL padding.
_text(parent, name::AbstractString)::String = _string(read(_dataset(parent, name)))

# The first and last element of a one-dimensional dataset, and its length. Reading the dataset whole
# would transfer the entire axis — 27,000 azimuth times for a real RSLC — where only the ends are used.
function _axis_bounds(parent, name::AbstractString)
    d = _dataset(parent, name)
    n = length(d)
    n > 0 || throw(ArgumentError("`$(HDF5.name(d))` is empty, so it has no bounds"))
    # Only the endpoints are read: the axis itself is tens of thousands of samples.
    return _as_number(d[1]), _as_number(d[n]), n
end

# `read_attribute` has no type-asserting form, so its `Any` is narrowed here rather than at each call.
function _units_epoch(parent, name::AbstractString)
    d = _dataset(parent, name)
    haskey(HDF5.attributes(d), "units") || throw(ArgumentError(
        "`$(HDF5.name(d))` has no `units` attribute, so the epoch its times are measured from is " *
        "unknown"))
    units = read_attribute(d, "units")
    return parse_cf_epoch(_string(units))
end

# Every read in the eager path goes through one open file handle. Opening a NISAR product costs about
# as much as the metadata reads themselves, so a reader that reopened per query would spend most of its
# time in `h5open`.
function read_identification(b::NisarBackend)
    return h5open(h -> read_identification(b, h), b.path, "r")
end

function read_identification(b::NisarBackend, h)
    g = _group(h, identification_path(b))
    # An absent optional field reads as empty rather than failing: a product missing one is still
    # usable, where one missing `absoluteOrbitNumber` is not.
    get_str(name)::String = haskey(g, name) ? _text(g, name) : ""
    return Identification(
        get_str("missionId"),
        get_str("productType"),
        _scalar_int(g, "absoluteOrbitNumber"),
        get_str("orbitPassDirection"),
        get_str("lookDirection"),
        get_str("zeroDopplerStartTime"),
        get_str("zeroDopplerEndTime"),
        get_str("boundingPolygon"),
    )
end

read_geometry(b::NisarBackend) = h5open(h -> read_geometry(b, h), b.path, "r")

function read_geometry(b::NisarBackend, h)
    b.product_type in GEOCODED_TYPES && throw(ArgumentError(
        "`$(b.path)` is a $(b.product_type) product, which stores its samples on a map grid and so " *
        "carries no slant-range/azimuth geometry. This package reads SLCs in radar geometry; open an " *
        "RSLC, or use a raster library for a geocoded product"))
    freq = _group(h, frequency_path(b))
    swaths = _group(h, swath_path(b))
    ident = _group(h, identification_path(b))

    near_range, far_range, nsamples = _axis_bounds(freq, "slantRange")
    sensing_start, sensing_stop, nlines = _axis_bounds(swaths, "zeroDopplerTime")

    look = _text(ident, "lookDirection")
    side = lowercase(look) == "left" ? LookLeft :
           lowercase(look) == "right" ? LookRight :
           throw(ArgumentError(
               "`$(b.path)` records lookDirection \"$look\"; expected \"Left\" or \"Right\""))

    return RadarGeometry(
        near_range,
        far_range,
        _scalar(freq, "slantRangeSpacing"),
        SPEED_OF_LIGHT / _scalar(freq, "processedCenterFrequency"),
        1 / _scalar(swaths, "zeroDopplerTimeSpacing"),
        sensing_start,
        sensing_stop,
        nlines,
        nsamples,
        side,
        # The azimuth times and the orbit's state vector times share one epoch in this format, so the
        # epoch is read once here and reported for both.
        _units_epoch(swaths, "zeroDopplerTime"),
    )
end

read_orbit(b::NisarBackend) = h5open(h -> read_orbit(b, h), b.path, "r")

function read_orbit(b::NisarBackend, h)
    g = _group(h, orbit_path(b))
    time = _read_vector(g, "time")
    # Stored (N, 3) row-major, which HDF5.jl presents as (3, N).
    pos = _read_matrix(g, "position")
    vel = _read_matrix(g, "velocity")
    size(pos, 1) == 3 || throw(ArgumentError(
        "orbit position in `$(b.path)` has leading dimension $(size(pos, 1)), expected 3"))
    axes(pos) == axes(vel) || throw(DimensionMismatch(
        "orbit position and velocity in `$(b.path)` have axes $(axes(pos)) and $(axes(vel))"))
    size(pos, 2) == length(time) || throw(DimensionMismatch(
        "orbit in `$(b.path)` has $(length(time)) times but $(size(pos, 2)) state vectors"))
    get_str(name, default)::String = haskey(g, name) ? _text(g, name) : default
    return StateVectors(
        time,
        _columns_as_svectors(pos),
        _columns_as_svectors(vel),
        _units_epoch(g, "time"),
        get_str("interpMethod", "Hermite"),
        get_str("orbitType", "Custom"),
    )
end

# A product may store these as `Float32`, so the element type is converted rather than propagated into
# `StateVectors`, whose fields are `Float64`. `convert` is a no-op when the file already holds `Float64`.
_read_vector(parent, name::AbstractString)::Vector{Float64} =
    convert(Vector{Float64}, read(_dataset(parent, name))::AbstractVector)

_read_matrix(parent, name::AbstractString)::Matrix{Float64} =
    convert(Matrix{Float64}, read(_dataset(parent, name))::AbstractMatrix)

# The (3, N) array HDF5.jl hands back, as N 3-vectors.
function _columns_as_svectors(a::AbstractMatrix{Float64})
    Base.require_one_based_indexing(a)
    return [SVector{3,Float64}(a[1, i], a[2, i], a[3, i]) for i in axes(a, 2)]
end
