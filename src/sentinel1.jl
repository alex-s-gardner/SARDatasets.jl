# The Sentinel-1 backend: TOPS products as a `.SAFE` directory or the zip of one.
#
# Two things make this format unlike the NISAR one. The geometry lives in per-subswath annotation XML
# rather than in one dataset tree, and the state vectors are not in the product at all — they come from
# a separately distributed POEORB/RESORB `.EOF` file, which is why `orbit` is a required keyword here and
# absent from the NISAR path.
#
# A TOPS product is also not one image. Each subswath is a strip of bursts at its own slant range, so
# there are two defensible geometries and this reads both: the mosaic spanning the subswaths, and one
# individual burst. They are `loadMetadataSlc` and `loadMetadata` in hyp3-autorift's `testGeogrid.py`,
# which is the reference this reproduces field for field.

# Sentinel-1 looks right, on every acquisition of the mission. The annotation carries no look-direction
# field, so unlike the NISAR path this is a constant rather than something read from the product.
const S1_LOOK_SIDE = LookRight

# State vectors are kept from one minute before a burst starts to one minute after it ends. Widening
# this changes which vectors an interpolator sees near the ends of the window, so it matches
# `s1reader`'s `PADDING_SHORT` exactly.
const S1_ORBIT_PADDING = 60.0

# The epoch every Sentinel-1 time is reported against is the burst's own sensing start less two days.
# The offset is arbitrary but not free to change: `s1reader` picks it in `get_burst_orbit` and
# `doppler_poly1d_to_lut2d` hardcodes a matching one, so a consumer pairing this package's times with
# an ISCE3-derived Doppler LUT needs the same value.
const S1_EPOCH_OFFSET_DAYS = 2

"""
    Sentinel1Backend <: AbstractSLCBackend

A Sentinel-1 TOPS product: the container, the orbit file, and which part of the product to describe.

`burst === nothing` selects the mosaic across `swaths`; an integer selects that 1-based burst of the
single subswath in `swaths`.
"""
struct Sentinel1Backend <: AbstractSLCBackend
    path::String
    orbit_path::String
    polarization::String
    swaths::Vector{Int}
    burst::Union{Nothing,Int}
end

"""
    is_safe_product(path) -> Bool

Whether `path` looks like a Sentinel-1 SAFE product: a `.SAFE` directory, or a zip holding one.

Recognition is by container shape rather than by file name, since a granule renamed on download is
still readable.
"""
function is_safe_product(path::AbstractString)
    isdir(path) && return isfile(joinpath(path, "manifest.safe"))
    isfile(path) || return false
    # Reading the central directory is enough; the measurement TIFFs are never touched.
    return try
        with_zip(r -> any(endswith("manifest.safe"), zip_names(r)), path)
    catch
        false
    end
end

# The archive is memory-mapped rather than read: a Sentinel-1 zip is a few gigabytes of measurement
# TIFFs around a few hundred kilobytes of annotation XML, and reading it into a buffer to reach the XML
# would cost the whole file's size in memory per call.
function with_zip(f, path::AbstractString)
    return open(path) do io
        f(ZipReader(mmap(io, Vector{UInt8}, filesize(path); grow = false)))
    end
end

zip_names(r) = [zip_name(r, i) for i in 1:zip_nentries(r)]

# An annotation file is `annotation/s1a-iw1-slc-hh-....xml`. The parallel `annotation/calibration/` and
# `annotation/rfi/` trees hold same-named files for products this package does not read, so matching on
# the file name alone would pick one of those up; the parent directory has to be `annotation` itself.
function _is_annotation(name::AbstractString, id_str::AbstractString)
    parts = split(name, '/')
    length(parts) >= 2 || return false
    return parts[end - 1] == "annotation" && occursin(id_str, parts[end])
end

_swath_id(swath::Integer, pol::AbstractString) = string("iw", swath, "-slc-", lowercase(pol))

"""
    annotation_xml(path, swath, polarization) -> String

The text of the annotation XML for one subswath and polarization of a SAFE product.
"""
function annotation_xml(path::AbstractString, swath::Integer, polarization::AbstractString)
    id_str = _swath_id(swath, polarization)
    if isdir(path)
        dir = joinpath(path, "annotation")
        isdir(dir) || throw(ArgumentError("`$path` has no `annotation` directory"))
        hits = filter(f -> occursin(id_str, f) && endswith(f, ".xml"), readdir(dir))
        isempty(hits) && throw(ArgumentError(
            "`$path` has no annotation file for subswath IW$swath polarization " *
            "$(uppercase(polarization))"))
        return read(joinpath(dir, only(sort(hits))), String)
    end
    return with_zip(path) do r
        hits = filter(n -> _is_annotation(n, id_str) && endswith(n, ".xml"), zip_names(r))
        isempty(hits) && throw(ArgumentError(
            "`$path` has no annotation file for subswath IW$swath polarization " *
            "$(uppercase(polarization))"))
        return zip_readentry(r, only(sort(hits)), String)
    end
end

"""
    safe_polarizations(path) -> Vector{String}

The polarizations a SAFE product carries, lowercase, in the order the annotation files sort.
"""
function safe_polarizations(path::AbstractString)
    names = isdir(path) ? (isdir(joinpath(path, "annotation")) ?
                           readdir(joinpath(path, "annotation")) : String[]) :
            with_zip(r -> [split(n, '/')[end] for n in zip_names(r)
                           if _is_annotation(n, "-slc-")], path)
    pols = String[]
    for n in names
        m = match(r"-slc-(hh|hv|vh|vv)-", n)
        m === nothing && continue
        m[1] in pols || push!(pols, m[1])
    end
    isempty(pols) && throw(ArgumentError(
        "`$path` carries no annotation files naming a polarization, so it is not a readable " *
        "Sentinel-1 SLC product"))
    return sort(pols)
end

"""
    default_polarization(path) -> String

The polarization to read when the caller names none.

Co-polarized channels are preferred over cross-polarized ones, and the order otherwise follows
`testGeogrid.py`'s `getPol`, so a dual-polarization product reads the same channel the reference does.
"""
function default_polarization(path::AbstractString)
    available = safe_polarizations(path)
    for pol in ("vv", "vh", "hh", "hv")
        pol in available && return pol
    end
    return first(available)
end

# The annotation records times as ISO-8601 to microseconds, which `DateTime` cannot hold. Times are
# therefore parsed to a `DateTime` plus a residual, and only ever reported as a `Float64` offset from
# the epoch, so no precision is lost to the millisecond truncation a bare `DateTime` would impose.
struct UtcTime
    datetime::DateTime      # truncated to whole seconds
    seconds::Float64        # fractional part, in [0, 1)
end

function parse_utc(s::AbstractString)
    t = strip(_string(s))
    t = replace(t, r"(Z|[+-]\d{2}:?\d{2})$" => "")
    m = match(r"^(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})(?:\.(\d+))?$", t)
    m === nothing && throw(ArgumentError(
        "cannot read a UTC time from \"$t\"; expected `YYYY-MM-DDThh:mm:ss[.ffffff]`"))
    whole = DateTime(m[1])
    frac = m[2] === nothing ? 0.0 : parse(Float64, string("0.", m[2]))
    return UtcTime(whole, frac)
end

# Seconds from one UTC instant to another. The whole-second difference is exact in `Int64`
# milliseconds and the sub-second parts are added afterwards, so the result does not inherit the
# rounding a `DateTime`-only subtraction would introduce.
function seconds_between(from::UtcTime, to::UtcTime)
    whole = Dates.value(to.datetime - from.datetime) / 1000
    return whole + (to.seconds - from.seconds)
end

Base.isless(a::UtcTime, b::UtcTime) =
    a.datetime == b.datetime ? a.seconds < b.seconds : a.datetime < b.datetime

# The epoch a product's times are reported against, truncated to the whole second.
#
# The reference's epoch is the anchor instant less two days exactly, keeping its microseconds; a
# `DateTime` holds only milliseconds, so carrying them here would lose precision silently. Truncating
# to the second instead and leaving the sub-second part in the times keeps `epoch + sensing_start` the
# exact instant the product records. The consequence is that this package's `sensing_start` differs
# from the reference's by the anchor's sub-second part — a constant both report identically once the
# epoch is added back.
epoch_of(t::UtcTime) = t.datetime - Dates.Day(S1_EPOCH_OFFSET_DAYS)

_findtext(node, path) = begin
    hit = findfirst(path, node)
    hit === nothing && throw(ArgumentError("the annotation XML has no `$path` element"))
    strip(nodecontent(hit))
end

_findfloat(node, path) = parse(Float64, _findtext(node, path))
_findint(node, path) = parse(Int, _findtext(node, path))

# One subswath's annotation, reduced to the geometry fields. Everything a burst needs that is shared
# across the subswath's bursts is read once here.
struct SubswathAnnotation
    swath::Int
    polarization::String
    absolute_orbit::Int
    mission::String
    product_type::String
    pass_direction::String
    starting_range::Float64
    range_pixel_spacing::Float64
    wavelength::Float64
    azimuth_time_interval::Float64
    lines_per_burst::Int
    samples_per_burst::Int
    burst_start::Vector{UtcTime}
end

function read_annotation(xml::AbstractString, swath::Integer)
    doc = parsexml(xml)
    r = root(doc)

    range_sampling_rate = _findfloat(r, "generalAnnotation/productInformation/rangeSamplingRate")
    radar_frequency = _findfloat(r, "generalAnnotation/productInformation/radarFrequency")
    slant_range_time = _findfloat(r, "imageAnnotation/imageInformation/slantRangeTime")

    burst_start = UtcTime[]
    for b in findall("swathTiming/burstList/burst", r)
        push!(burst_start, parse_utc(_findtext(b, "azimuthTime")))
    end
    isempty(burst_start) && throw(ArgumentError(
        "the annotation for subswath IW$swath lists no bursts"))

    # These three conversions must not be reordered: they reproduce ISCE3's values to the bit, and a
    # different association of the same operations would change the last digit.
    return SubswathAnnotation(
        Int(swath),
        lowercase(_findtext(r, "adsHeader/polarisation")),
        _findint(r, "adsHeader/absoluteOrbitNumber"),
        _findtext(r, "adsHeader/missionId"),
        _findtext(r, "adsHeader/productType"),
        _findtext(r, "generalAnnotation/productInformation/pass"),
        slant_range_time * SPEED_OF_LIGHT / 2,
        SPEED_OF_LIGHT / (2 * range_sampling_rate),
        SPEED_OF_LIGHT / radar_frequency,
        _findfloat(r, "imageAnnotation/imageInformation/azimuthTimeInterval"),
        _findint(r, "swathTiming/linesPerBurst"),
        _findint(r, "swathTiming/samplesPerBurst"),
        burst_start,
    )
end

read_annotation(b::Sentinel1Backend, swath::Integer) =
    read_annotation(annotation_xml(b.path, swath, b.polarization), swath)

# A burst's last line is a whole number of azimuth lines after its first, rather than a time read from
# the product: the annotation gives only the start of each burst, and every burst carries
# `linesPerBurst` lines.
_burst_stop(a::SubswathAnnotation, start::UtcTime) =
    _advance(start, (a.lines_per_burst - 1) * a.azimuth_time_interval)

_advance(t::UtcTime, seconds::Real) = UtcTime(t.datetime, t.seconds + seconds)

# The instant every time in the result is measured from: the first burst of the lowest-numbered
# subswath asked for. It is that burst rather than the earliest line in the product because the
# reference takes its epoch from `bursts[0].orbit` of `swaths[0]` — and in a TOPS product IW1 does not
# start first, so the two differ by about two seconds.
_anchor(b::Sentinel1Backend, a::SubswathAnnotation) =
    b.burst === nothing ? first(a.burst_start) : a.burst_start[_check_burst(b, a)]

function read_identification(b::Sentinel1Backend)
    a = read_annotation(b, first(b.swaths))
    first_line, last_line = if b.burst === nothing
        # The bounds of what this SLC describes, which for a mosaic spans the subswaths rather than
        # being the one subswath's `productFirst/LastLineUtcTime`.
        _mosaic_window([read_annotation(b, s) for s in b.swaths])
    else
        start = a.burst_start[_check_burst(b, a)]
        start, _burst_stop(a, start)
    end
    return Identification(
        a.mission,
        a.product_type,
        a.absolute_orbit,
        lowercase(a.pass_direction),
        S1_LOOK_SIDE == LookLeft ? "Left" : "Right",
        _utc_string(first_line),
        _utc_string(last_line),
        # The product's footprint is in `manifest.safe` rather than the annotation, and describing a
        # mosaic or a single burst by the whole product's polygon would be wrong for both.
        "",
    )
end

function _utc_string(t::UtcTime)
    # Normalized so a fractional part carried past one second reads as the time it denotes.
    extra = floor(Int, t.seconds)
    frac = t.seconds - extra
    stamp = Dates.format(t.datetime + Dates.Second(extra), "yyyy-mm-ddTHH:MM:SS")
    return string(stamp, ".", lpad(round(Int, frac * 1_000_000), 6, '0'))
end

function _check_burst(b::Sentinel1Backend, a::SubswathAnnotation)
    n = length(a.burst_start)
    i = b.burst
    1 <= i <= n || throw(ArgumentError(
        "`$(b.path)` subswath IW$(a.swath) has $n bursts, so burst $i does not exist"))
    return i
end

# The mosaic's geometry, matching `loadMetadataSlc`. The near-range subswath sets the range origin, the
# azimuth spacing and the wavelength; the sensing window is the union over the subswaths; and the width
# reaches from the near subswath's first sample to the far subswath's last.
function _mosaic_geometry(b::Sentinel1Backend)
    annotations = [read_annotation(b, s) for s in b.swaths]
    near = first(annotations)
    far = last(annotations)

    prf = 1 / near.azimuth_time_interval
    starting_range = near.starting_range
    range_pixel_spacing = near.range_pixel_spacing

    start, stop = _mosaic_window(annotations)

    # IW subswaths are numbered outward from nadir, so the highest-numbered one asked for is the one
    # that sets the far edge. `b.swaths` is sorted, which is what makes `far` that subswath.
    far.starting_range >= starting_range || throw(ArgumentError(
        "`$(b.path)` puts subswath IW$(far.swath) nearer in range than IW$(near.swath), so the " *
        "subswaths are not ordered outward and the mosaic's width cannot be derived"))
    nsamples = round(Int, (far.starting_range - starting_range) / range_pixel_spacing) +
               far.samples_per_burst
    far_range = starting_range + (nsamples - 1.0) * range_pixel_spacing
    nlines = round(Int, seconds_between(start, stop) * prf) + 1

    epoch = epoch_of(_anchor(b, near))
    origin = UtcTime(epoch, 0.0)
    return RadarGeometry(
        starting_range,
        far_range,
        range_pixel_spacing,
        near.wavelength,
        prf,
        seconds_between(origin, start),
        seconds_between(origin, stop),
        nlines,
        nsamples,
        S1_LOOK_SIDE,
        epoch,
    )
end

# The sensing window a mosaic covers: the earliest burst start of any subswath to the latest burst end.
# Every subswath's extent is measured with the near subswath's azimuth spacing, as the reference uses a
# single PRF for the mosaic.
function _mosaic_window(annotations::AbstractVector{SubswathAnnotation})
    near = first(annotations)
    start = first(near.burst_start)
    stop = _advance(last(near.burst_start),
                    (near.lines_per_burst - 1) * near.azimuth_time_interval)
    for a in annotations
        start = min(start, first(a.burst_start))
        stop = max(stop, _advance(last(a.burst_start),
                                  (a.lines_per_burst - 1) * near.azimuth_time_interval))
    end
    return start, stop
end

# One burst's geometry, matching `loadMetadata`: the burst's own start time and the subswath's full
# burst extent, with no mosaicking.
function _burst_geometry(b::Sentinel1Backend)
    a = read_annotation(b, only(b.swaths))
    i = _check_burst(b, a)
    start = a.burst_start[i]
    stop = _burst_stop(a, start)

    prf = 1 / a.azimuth_time_interval
    far_range = a.starting_range + (a.samples_per_burst - 1.0) * a.range_pixel_spacing
    epoch = epoch_of(start)
    origin = UtcTime(epoch, 0.0)
    return RadarGeometry(
        a.starting_range,
        far_range,
        a.range_pixel_spacing,
        a.wavelength,
        prf,
        seconds_between(origin, start),
        seconds_between(origin, stop),
        a.lines_per_burst,
        a.samples_per_burst,
        S1_LOOK_SIDE,
        epoch,
    )
end

read_geometry(b::Sentinel1Backend) =
    b.burst === nothing ? _mosaic_geometry(b) : _burst_geometry(b)

"""
    read_eof_state_vectors(path) -> Vector

Every orbit state vector in a POEORB or RESORB `.EOF` file, as `(UtcTime, position, velocity)`.

The file's `UTC` fields carry a `UTC=` prefix, which is stripped; a reader that kept it would fail to
parse rather than misread, but the times are the only thing tying the orbit to the product's clock.
"""
function read_eof_state_vectors(path::AbstractString)
    doc = readxml(path)
    list = findfirst("//Data_Block/List_of_OSVs", root(doc))
    list === nothing && throw(ArgumentError(
        "`$path` has no `Data_Block/List_of_OSVs` element, so it is not a Sentinel-1 orbit file"))
    out = Tuple{UtcTime,SVector{3,Float64},SVector{3,Float64}}[]
    for osv in eachelement(list)
        time = parse_utc(replace(_findtext(osv, "UTC"), r"^UTC=" => ""))
        pos = SVector{3,Float64}(_findfloat(osv, "X"), _findfloat(osv, "Y"), _findfloat(osv, "Z"))
        vel = SVector{3,Float64}(_findfloat(osv, "VX"), _findfloat(osv, "VY"), _findfloat(osv, "VZ"))
        push!(out, (time, pos, vel))
    end
    isempty(out) && throw(ArgumentError("`$path` lists no orbit state vectors"))
    return out
end

# The state vectors kept are those within `S1_ORBIT_PADDING` of the anchor burst's own window — not of
# the whole mosaic's. `s1reader` builds one orbit per burst and the reference then takes the first
# burst's, so a mosaic's orbit covers that burst rather than the full acquisition; widening it here
# would hand an interpolator a different set of vectors than ISCE3 sees.
function read_orbit(b::Sentinel1Backend)
    a = read_annotation(b, first(b.swaths))
    start = _anchor(b, a)
    stop = _burst_stop(a, start)

    epoch = UtcTime(epoch_of(start), 0.0)
    times = Float64[]
    positions = SVector{3,Float64}[]
    velocities = SVector{3,Float64}[]
    for (t, pos, vel) in read_eof_state_vectors(b.orbit_path)
        dt_start = seconds_between(start, t)
        dt_stop = seconds_between(stop, t)
        dt_stop > S1_ORBIT_PADDING && break
        dt_start > -S1_ORBIT_PADDING || continue
        push!(times, seconds_between(epoch, t))
        push!(positions, pos)
        push!(velocities, vel)
    end
    isempty(times) && throw(ArgumentError(
        "`$(b.orbit_path)` has no state vectors within $(S1_ORBIT_PADDING) s of the acquisition " *
        "window of `$(b.path)`; it is probably the orbit file of a different granule"))

    return StateVectors(
        times, positions, velocities, epoch.datetime,
        # An EOF carries tabulated vectors with no interpolation method of its own. Hermite is what
        # `s1reader` builds an `isce3.core.Orbit` with by default.
        "Hermite",
        _eof_kind(b.orbit_path),
    )
end

# POEORB is the precise orbit, published days later; RESORB is the restituted one available at once.
# Which was used changes the geometry at the metre level, so it is reported rather than dropped.
function _eof_kind(path::AbstractString)
    name = uppercase(basename(path))
    occursin("POEORB", name) && return "POEORB"
    occursin("RESORB", name) && return "RESORB"
    return "Custom"
end

"""
    nbursts(path; swath = 1, polarization = nothing) -> Int

How many bursts one subswath of a Sentinel-1 product carries.

Needed to address bursts by index without opening each in turn, and it reads only the annotation XML.
"""
function nbursts(path::AbstractString; swath::Integer = 1, polarization = nothing)
    pol = polarization === nothing ? default_polarization(path) : lowercase(String(polarization))
    return length(read_annotation(annotation_xml(path, swath, pol), swath).burst_start)
end
