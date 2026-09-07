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
#
# An annotation document is a couple of megabytes and holds a few dozen scalars the reader wants, so the
# subswaths a product is opened for are parsed once into `Sentinel1Product` and every geometry, burst and
# orbit is derived from that. Nothing here re-reads the container.

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

const S1_SUBSWATHS = 1:3

"""
    is_safe_product(path) -> Bool

Whether `path` looks like a Sentinel-1 SAFE product: a `.SAFE` directory, or a zip holding one.

Recognition is by container shape rather than by file name, since a granule renamed on download is
still readable.
"""
function is_safe_product(path::AbstractString)
    isdir(path) && return isfile(joinpath(path, "manifest.safe"))
    isfile(path) || return false
    # Reading the central directory is enough; the measurement TIFFs are never touched. The names are
    # generated rather than collected so the scan stops at the manifest, which is near the front.
    return try
        with_zip(path) do r
            any(i -> endswith(zip_name(r, i), "manifest.safe"), 1:zip_nentries(r))
        end
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

# Every entry name in an archive. The readers filter while listing rather than calling this, so it
# exists for a caller that genuinely wants the whole central directory.
zip_names(r) = [zip_name(r, i) for i in 1:zip_nentries(r)]

# An annotation file is `annotation/s1a-iw1-slc-hh-....xml`. The parallel `annotation/calibration/` and
# `annotation/rfi/` trees hold same-named files for products this package does not read, so matching on
# the file name alone would pick one of those up; the parent directory has to be `annotation` itself.
function _is_annotation(name::AbstractString, id_str::AbstractString)
    endswith(name, ".xml") || return false
    slash = findlast('/', name)
    slash === nothing && return false
    occursin(id_str, SubString(name, slash + 1)) || return false
    rest = SubString(name, 1, slash - 1)
    prev = findlast('/', rest)
    return SubString(rest, something(prev, 0) + 1) == "annotation"
end

# The file-name part of a zip entry. Only for a name `_is_annotation` accepted, which is what guarantees
# the separator is there.
function _zip_basename(name::AbstractString)
    slash = findlast('/', name)
    slash === nothing && throw(ArgumentError("`$name` names no directory, so it has no file part"))
    return SubString(name, slash + 1)
end

_swath_id(swath::Integer, pol::AbstractString) = string("iw", swath, "-slc-", lowercase(pol))

@noinline _no_annotation(path, swath, pol) = throw(ArgumentError(
    "`$path` has no annotation file for subswath IW$swath polarization $(uppercase(pol))"))

"""
    annotation_xml(path, swath, polarization) -> String

The text of the annotation XML for one subswath and polarization of a SAFE product.
"""
annotation_xml(path::AbstractString, swath::Integer, polarization::AbstractString) =
    only(annotation_xml(path, (swath,), polarization))

"""
    AnnotationTree

The `annotation/` entries of an open SAFE container, listed once.

A zip must be memory-mapped and its central directory parsed to reach any entry, and a product's
polarizations and its per-subswath annotations are both answered from that one listing. Holding it makes
opening a product one traversal rather than one per query.
"""
struct AnnotationTree{C}
    path::String
    container::C          # the `ZipReader`, or the `annotation` directory of a `.SAFE`
    names::Vector{String} # entry names as the container addresses them
end

"""
    with_annotations(f, path)

Call `f` on the [`AnnotationTree`](@ref) of the SAFE product at `path`, however it is packaged.
"""
function with_annotations(f, path::AbstractString)
    if isdir(path)
        dir = joinpath(path, "annotation")
        isdir(dir) || throw(ArgumentError("`$path` has no `annotation` directory"))
        names = filter(n -> occursin("-slc-", n) && endswith(n, ".xml"), readdir(dir))
        return f(AnnotationTree(String(path), dir, names))
    end
    return with_zip(path) do r
        names = [zip_name(r, i) for i in 1:zip_nentries(r) if _is_annotation(zip_name(r, i), "-slc-")]
        return f(AnnotationTree(String(path), r, names))
    end
end

# The name of a zip entry or of a file in the `annotation` directory, reduced to the file part the
# subswath and polarization are named in.
_entry_file(r::AnnotationTree{<:AbstractString}, name) = name
_entry_file(::AnnotationTree, name) = _zip_basename(name)

_read_entry(r::AnnotationTree{<:AbstractString}, name) = read(joinpath(r.container, name), String)
_read_entry(r::AnnotationTree, name) = zip_readentry(r.container, name, String)

"""
    read_annotation_entry(r::AnnotationTree, swath, polarization) -> String

The annotation XML of one subswath and polarization, from an already-listed container.
"""
function read_annotation_entry(r::AnnotationTree, swath::Integer, polarization::AbstractString)
    id_str = _swath_id(swath, polarization)
    hits = filter(n -> occursin(id_str, _entry_file(r, n)), r.names)
    isempty(hits) && _no_annotation(r.path, swath, polarization)
    return _read_entry(r, only(sort!(hits)))
end

"""
    annotation_xml(path, swaths, polarization) -> Vector{String}

The annotation XML of several subswaths, read in one pass over the container.
"""
annotation_xml(path::AbstractString, swaths, polarization::AbstractString) =
    with_annotations(r -> annotation_xml(r, swaths, polarization), path)

annotation_xml(r::AnnotationTree, swaths, polarization::AbstractString) =
    [read_annotation_entry(r, s, polarization) for s in swaths]

"""
    safe_polarizations(path) -> Vector{String}
    safe_polarizations(r::AnnotationTree) -> Vector{String}

The polarizations a SAFE product carries, lowercase, in the order the annotation files sort.
"""
safe_polarizations(path::AbstractString) = with_annotations(safe_polarizations, path)

function safe_polarizations(r::AnnotationTree)
    pols = String[]
    for n in r.names
        m = match(r"-slc-(hh|hv|vh|vv)-", _entry_file(r, n))
        m === nothing && continue
        m[1] in pols || push!(pols, m[1])
    end
    isempty(pols) && throw(ArgumentError(
        "`$(r.path)` carries no annotation files naming a polarization, so it is not a readable " *
        "Sentinel-1 SLC product"))
    return sort!(pols)
end

"""
    default_polarization(path) -> String

The polarization to read when the caller names none.

Co-polarized channels are preferred over cross-polarized ones, and the order otherwise follows
`testGeogrid.py`'s `getPol`, so a dual-polarization product reads the same channel the reference does.
"""
default_polarization(path::AbstractString) = _preferred_polarization(safe_polarizations(path))

function _preferred_polarization(available)
    for pol in ("vv", "vh", "hh", "hv")
        pol in available && return pol
    end
    return first(available)
end

@noinline _no_element(path) = throw(ArgumentError("the annotation XML has no `$path` element"))

function _findtext(node, path)
    hit = findfirst(path, node)
    hit === nothing && _no_element(path)
    return strip(nodecontent(hit))
end

_findfloat(node, path) = parse(Float64, _findtext(node, path))
_findint(node, path) = parse(Int, _findtext(node, path))

# A whitespace-separated list of integers, as the per-line valid-sample arrays are written.
_findints(node, path) = [parse(Int, s) for s in eachsplit(_findtext(node, path))]

# The epoch a product's times are reported against, truncated to the whole second.
#
# The reference's epoch is the anchor instant less two days exactly, keeping its microseconds; a
# `DateTime` holds only milliseconds, so carrying them here would lose precision silently. Truncating
# to the second instead and leaving the sub-second part in the times keeps `epoch + sensing_start` the
# exact instant the product records. The consequence is that this package's `sensing_start` differs
# from the reference's by the anchor's sub-second part — a constant both report identically once the
# epoch is added back.
epoch_of(t::UtcTime) = t.datetime - Dates.Day(S1_EPOCH_OFFSET_DAYS)

"""
    SubswathAnnotation

One subswath's annotation, reduced to the geometry fields.

Everything a burst needs that is shared across the subswath's bursts is held once here, so a burst
costs no parsing of its own.
"""
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
    # The valid region of each burst, parallel to `burst_start`. See `read_valid_region`; these are
    # 1-based inclusive bounds into the burst's own `lines_per_burst` x `samples_per_burst` extent.
    first_valid_line::Vector{Int}
    last_valid_line::Vector{Int}
    first_valid_sample::Vector{Int}
    last_valid_sample::Vector{Int}
end

nbursts(a::SubswathAnnotation) = length(a.burst_start)

"""
    valid_lines(a::SubswathAnnotation, burst) -> UnitRange{Int}
    valid_samples(a::SubswathAnnotation, burst) -> UnitRange{Int}

The rows and columns of one burst that carry data, as 1-based inclusive ranges.

A burst's raster is `lines_per_burst` by `samples_per_burst`, but only this sub-rectangle of it was
imaged; the margin outside carries whatever the processor left there, which is not zero. So a reader
combining bursts must clip to these, and a caller reading a burst's raster directly gets the margin
unless it does.
"""
valid_lines(a::SubswathAnnotation, burst::Integer) =
    a.first_valid_line[burst]:a.last_valid_line[burst]
valid_samples(a::SubswathAnnotation, burst::Integer) =
    a.first_valid_sample[burst]:a.last_valid_sample[burst]

# The valid region of one burst, from the per-line `firstValidSample`/`lastValidSample` arrays.
#
# Each array holds one entry per line of the burst, negative where the line carries nothing. This
# reduces them to a rectangle the way `s1reader`'s `burst_from_xml` does: the valid lines are the run
# of non-negative entries, and the sample bounds are taken from the *first and last* of those lines
# only — the widest start and the narrowest end of the two, so the rectangle is inside both.
#
# Reading only two lines would silently widen the rectangle if an interior line were narrower than
# both, so that is checked rather than assumed: the arrays are constant across the valid lines of
# every granule measured, and a product where they are not would otherwise contribute margin samples
# as though they were data.
function read_valid_region(node, swath::Integer, burst::Integer, lines_per_burst::Integer)
    first_sample = _findints(node, "firstValidSample")
    last_sample = _findints(node, "lastValidSample")

    length(first_sample) == lines_per_burst || throw(ArgumentError(
        "burst $burst of subswath IW$swath gives $(length(first_sample)) `firstValidSample` " *
        "entries for $lines_per_burst lines; there must be one per line"))
    length(last_sample) == length(first_sample) || throw(ArgumentError(
        "burst $burst of subswath IW$swath gives $(length(first_sample)) `firstValidSample` " *
        "entries but $(length(last_sample)) `lastValidSample` entries; there must be one of each " *
        "per line"))

    lines = findall(>=(0), first_sample)
    isempty(lines) && throw(ArgumentError(
        "burst $burst of subswath IW$swath has no valid lines, so it carries no data"))
    # A valid region split into several runs is not a rectangle, and reducing it to one would claim
    # lines that carry nothing.
    lines == first(lines):last(lines) || throw(ArgumentError(
        "burst $burst of subswath IW$swath marks valid lines in more than one run " *
        "($(length(lines)) lines between $(first(lines)) and $(last(lines))), so its valid region " *
        "is not a rectangle"))

    lo, hi = first(lines), last(lines)
    first_valid = max(first_sample[lo], first_sample[hi])
    last_valid = min(last_sample[lo], last_sample[hi])
    first_valid <= last_valid || throw(ArgumentError(
        "burst $burst of subswath IW$swath has no sample valid on both its first and last valid " *
        "line, so its valid region is empty"))

    for i in lines
        (first_sample[i] == first_sample[lo] && last_sample[i] == last_sample[lo]) || throw(
            ArgumentError(
                "burst $burst of subswath IW$swath varies its valid sample range from " *
                "$(first_sample[lo])-$(last_sample[lo]) on line $lo to " *
                "$(first_sample[i])-$(last_sample[i]) on line $i. This reader takes a burst's " *
                "valid region to be a rectangle, which it is in every product measured"))
    end

    # `lo` and `hi` index the per-line arrays and so already count from one; the sample bounds are the
    # annotation's own values, which count from zero.
    return lo, hi, first_valid + 1, last_valid + 1
end

function read_annotation(xml::AbstractString, swath::Integer)
    doc = parsexml(xml)
    r = root(doc)

    range_sampling_rate = _findfloat(r, "generalAnnotation/productInformation/rangeSamplingRate")
    radar_frequency = _findfloat(r, "generalAnnotation/productInformation/radarFrequency")
    slant_range_time = _findfloat(r, "imageAnnotation/imageInformation/slantRangeTime")

    bursts = findall("swathTiming/burstList/burst", r)
    isempty(bursts) && throw(ArgumentError(
        "the annotation for subswath IW$swath lists no bursts"))
    lines_per_burst = _findint(r, "swathTiming/linesPerBurst")

    n = length(bursts)
    burst_start = Vector{UtcTime}(undef, n)
    first_valid_line = Vector{Int}(undef, n)
    last_valid_line = Vector{Int}(undef, n)
    first_valid_sample = Vector{Int}(undef, n)
    last_valid_sample = Vector{Int}(undef, n)
    for (i, b) in enumerate(bursts)
        burst_start[i] = parse_utc(_findtext(b, "azimuthTime"))
        first_valid_line[i], last_valid_line[i], first_valid_sample[i], last_valid_sample[i] =
            read_valid_region(b, swath, i, lines_per_burst)
    end

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
        lines_per_burst,
        _findint(r, "swathTiming/samplesPerBurst"),
        burst_start,
        first_valid_line,
        last_valid_line,
        first_valid_sample,
        last_valid_sample,
    )
end

"""
    Sentinel1Product(path; orbit, polarization = nothing, swaths = 1:3)

A Sentinel-1 TOPS product with its annotation already parsed.

`swaths` names the subswaths to read; each is parsed once and held, so every geometry, burst and orbit
derived from this product costs no further parsing of the container. Build the acquisitions with
[`open_slc`](@ref) or [`bursts`](@ref) rather than reaching into the fields.
"""
struct Sentinel1Product
    path::String
    orbit_path::String
    polarization::String
    # Parallel to `swaths`: `annotations[i]` is the annotation of subswath `swaths[i]`. `swaths` is
    # sorted ascending, which is what makes the first entry the near-range subswath of a mosaic.
    swaths::Vector{Int}
    annotations::Vector{SubswathAnnotation}
end

Sentinel1Product(path::AbstractString; orbit::AbstractString,
                 polarization = nothing, swaths = S1_SUBSWATHS) =
    with_annotations(r -> Sentinel1Product(r; orbit, polarization, swaths), path)

# One listing of the container serves the polarization query and every subswath's annotation.
function Sentinel1Product(r::AnnotationTree; orbit::AbstractString,
                          polarization = nothing, swaths = S1_SUBSWATHS)
    # Sorted here so `annotations` is near-range first, which is what makes a mosaic's first entry the
    # subswath that sets its range origin.
    selected = sort!(collect(Int, swaths))
    available = safe_polarizations(r)
    pol = polarization === nothing ? _preferred_polarization(available) :
          lowercase(String(polarization))
    pol in available || _no_polarization(r.path, available, pol)
    annotations = [read_annotation(read_annotation_entry(r, s, pol), s) for s in selected]
    return Sentinel1Product(r.path, String(orbit), pol, selected, annotations)
end

@noinline _no_polarization(path, available, pol) = throw(ArgumentError(
    "`$path` carries polarization$(length(available) == 1 ? " " : "s ")" *
    "$(join(uppercase.(available), ", ")), not $(uppercase(pol))"))

Base.show(io::IO, p::Sentinel1Product) =
    print(io, "Sentinel1Product(", basename(p.path), ", IW",
          join(p.swaths, "+"), ", ", uppercase(p.polarization), ")")

"""
    annotation(p::Sentinel1Product, swath) -> SubswathAnnotation

The parsed annotation of one subswath of `p`, which must be among the ones it was opened for.
"""
function annotation(p::Sentinel1Product, swath::Integer)
    i = findfirst(==(Int(swath)), p.swaths)
    i === nothing && throw(ArgumentError(
        "`$(p.path)` was opened for subswath$(length(p.swaths) == 1 ? " " : "s ")" *
        "$(join(("IW$s" for s in p.swaths), ", ")), not IW$swath"))
    return p.annotations[i]
end

"""
    nbursts(p::Sentinel1Product, swath) -> Int
    nbursts(path; swath = 1, polarization = nothing) -> Int

How many bursts one subswath of a Sentinel-1 product carries.

Needed to address bursts by index without opening each in turn, and it reads only the annotation XML.
"""
nbursts(p::Sentinel1Product, swath::Integer) = nbursts(annotation(p, swath))

function nbursts(path::AbstractString; swath::Integer = 1, polarization = nothing)
    return with_annotations(path) do r
        pol = polarization === nothing ? _preferred_polarization(safe_polarizations(r)) :
              lowercase(String(polarization))
        return nbursts(read_annotation(read_annotation_entry(r, swath, pol), swath))
    end
end

"""
    Sentinel1Backend <: AbstractSLCBackend

Which part of a parsed [`Sentinel1Product`](@ref) an [`SLC`](@ref) describes.

A mosaic spans every subswath the product was opened for, so it has no subswath or burst of its own and
both are `nothing`. An individual burst names the subswath it belongs to and its own 1-based index.
"""
struct Sentinel1Backend <: AbstractSLCBackend
    product::Sentinel1Product
    swath::Union{Nothing,Int}
    burst::Union{Nothing,Int}
end

Sentinel1Backend(p::Sentinel1Product) = Sentinel1Backend(p, nothing, nothing)

_path(b::Sentinel1Backend) = b.product.path

# Whether this describes the mosaic rather than one burst. Both fields are `nothing` together, so either
# answers it; `swath` is the one asked because it is what the mosaic lacks.
_is_mosaic(b::Sentinel1Backend) = b.swath === nothing

# The subswath whose annotation sets the times: the near-range one for a mosaic, the burst's own
# otherwise.
_leading_annotation(b::Sentinel1Backend) =
    _is_mosaic(b) ? first(b.product.annotations) : annotation(b.product, b.swath)

# A burst's last line is a whole number of azimuth lines after its first, rather than a time read from
# the product: the annotation gives only the start of each burst, and every burst carries
# `linesPerBurst` lines.
_burst_stop(a::SubswathAnnotation, start::UtcTime) =
    _advance(start, (a.lines_per_burst - 1) * a.azimuth_time_interval)

# The instant every time in the result is measured from: the first burst of the lowest-numbered
# subswath asked for. It is that burst rather than the earliest line in the product because the
# reference takes its epoch from `bursts[0].orbit` of `swaths[0]` — and in a TOPS product IW1 does not
# start first, so the two differ by about two seconds.
function _anchor(b::Sentinel1Backend)
    a = _leading_annotation(b)
    return _is_mosaic(b) ? first(a.burst_start) : a.burst_start[_check_burst(b, a)]
end

function _check_burst(b::Sentinel1Backend, a::SubswathAnnotation)
    n = nbursts(a)
    i = b.burst
    1 <= i <= n || throw(ArgumentError(
        "`$(_path(b))` subswath IW$(a.swath) has $n bursts, so burst $i does not exist"))
    return i
end

function read_identification(b::Sentinel1Backend)
    a = _leading_annotation(b)
    first_line, last_line = if _is_mosaic(b)
        # The bounds of what this SLC describes, which for a mosaic spans the subswaths rather than
        # being the one subswath's `productFirst/LastLineUtcTime`.
        _mosaic_window(b.product.annotations)
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

# The mosaic's geometry, matching `loadMetadataSlc`. The near-range subswath sets the range origin, the
# azimuth spacing and the wavelength; the sensing window is the union over the subswaths; and the width
# reaches from the near subswath's first sample to the far subswath's last.
function _mosaic_geometry(b::Sentinel1Backend)
    annotations = b.product.annotations
    near = first(annotations)
    far = last(annotations)

    prf = 1 / near.azimuth_time_interval
    starting_range = near.starting_range
    range_pixel_spacing = near.range_pixel_spacing

    start, stop = _mosaic_window(annotations)

    # IW subswaths are numbered outward from nadir, so the highest-numbered one asked for is the one
    # that sets the far edge. `swaths` is sorted, which is what makes `far` that subswath.
    far.starting_range >= starting_range || throw(ArgumentError(
        "`$(_path(b))` puts subswath IW$(far.swath) nearer in range than IW$(near.swath), so the " *
        "subswaths are not ordered outward and the mosaic's width cannot be derived"))
    nsamples = round(Int, (far.starting_range - starting_range) / range_pixel_spacing) +
               far.samples_per_burst
    far_range = starting_range + (nsamples - 1.0) * range_pixel_spacing
    nlines = round(Int, seconds_between(start, stop) * prf) + 1

    epoch = epoch_of(_anchor(b))
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
    a = annotation(b.product, b.swath)
    start = a.burst_start[_check_burst(b, a)]
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
    _is_mosaic(b) ? _mosaic_geometry(b) : _burst_geometry(b)

# The `UTC=` prefix an EOF's time fields carry. A reader that kept it would fail to parse rather than
# misread, but the times are the only thing tying the orbit to the product's clock.
_strip_utc_prefix(s::AbstractString) =
    startswith(s, "UTC=") ? SubString(s, 5) : SubString(s, 1, lastindex(s))

"""
    read_eof_state_vectors(path; from = nothing, to = nothing, padding = 0.0) -> StateVectorTable

The orbit state vectors of a POEORB or RESORB `.EOF` file.

With `from` and `to` given, only the vectors within `padding` seconds of that window are kept, and the
scan stops once past it. A POEORB covers just over a day at one vector every ten seconds while an
acquisition needs a few, so filtering during the scan rather than after it is what keeps reading one
cheap: a record whose time falls outside the window has its six numbers skipped entirely.
"""
function read_eof_state_vectors(path::AbstractString; from::Union{Nothing,UtcTime} = nothing,
                                to::Union{Nothing,UtcTime} = nothing, padding::Real = 0.0)
    doc = readxml(path)
    list = findfirst("//Data_Block/List_of_OSVs", root(doc))
    list === nothing && throw(ArgumentError(
        "`$path` has no `Data_Block/List_of_OSVs` element, so it is not a Sentinel-1 orbit file"))

    times = UtcTime[]
    positions = SVector{3,Float64}[]
    velocities = SVector{3,Float64}[]
    any_record = false
    for osv in eachelement(list)
        any_record = true
        t = _osv_time(osv)
        if to !== nothing && seconds_between(to, t) > padding
            # The vectors are in ascending time order, so nothing later can fall in the window.
            break
        end
        # Strictly inside the leading edge, matching `s1reader`: a vector landing exactly `padding`
        # before the window is not kept, and admitting it would hand an interpolator one more vector
        # than ISCE3 sees.
        from !== nothing && seconds_between(from, t) <= -padding && continue
        pos, vel = _osv_vectors(osv)
        push!(times, t)
        push!(positions, pos)
        push!(velocities, vel)
    end
    any_record || throw(ArgumentError("`$path` lists no orbit state vectors"))
    return StateVectorTable(times, positions, velocities)
end

# One record's time, found by scanning its children rather than by an XPath query: the query builds and
# frees a node set per call, which over the thousands of records in a POEORB dominates the read.
function _osv_time(osv)
    for f in eachelement(osv)
        nodename(f) == "UTC" && return parse_utc(_strip_utc_prefix(strip(nodecontent(f))))
    end
    throw(ArgumentError("an `OSV` record in the orbit file has no `UTC` field"))
end

function _osv_vectors(osv)
    x = y = z = vx = vy = vz = 0.0
    found = 0
    for f in eachelement(osv)
        n = nodename(f)
        # The six names are one or two characters; the record's other fields are longer and are not
        # parsed, since their contents are not numbers.
        length(n) <= 2 || continue
        v = parse(Float64, strip(nodecontent(f)))
        if n == "X"
            x = v
        elseif n == "Y"
            y = v
        elseif n == "Z"
            z = v
        elseif n == "VX"
            vx = v
        elseif n == "VY"
            vy = v
        elseif n == "VZ"
            vz = v
        else
            continue
        end
        found += 1
    end
    found == 6 || throw(ArgumentError(
        "an `OSV` record in the orbit file is missing one of its X/Y/Z/VX/VY/VZ fields"))
    return SVector(x, y, z), SVector(vx, vy, vz)
end

# The state vectors kept are those within `S1_ORBIT_PADDING` of the anchor burst's own window — not of
# the whole mosaic's. `s1reader` builds one orbit per burst and the reference then takes the first
# burst's, so a mosaic's orbit covers that burst rather than the full acquisition; widening it here
# would hand an interpolator a different set of vectors than ISCE3 sees.
function read_orbit(b::Sentinel1Backend)
    a = _leading_annotation(b)
    start = _anchor(b)
    stop = _burst_stop(a, start)
    orbit_path = b.product.orbit_path

    table = read_eof_state_vectors(orbit_path; from = start, to = stop,
                                   padding = S1_ORBIT_PADDING)
    isempty(table.time) && throw(ArgumentError(
        "`$orbit_path` has no state vectors within $(S1_ORBIT_PADDING) s of the acquisition " *
        "window of `$(_path(b))`; it is probably the orbit file of a different granule"))

    epoch = UtcTime(epoch_of(start), 0.0)
    return StateVectors(
        [seconds_between(epoch, t) for t in table.time],
        table.position,
        table.velocity,
        epoch.datetime,
        # An EOF carries tabulated vectors with no interpolation method of its own. Hermite is what
        # `s1reader` builds an `isce3.core.Orbit` with by default.
        "Hermite",
        _eof_kind(orbit_path),
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

# Reaching a subswath's measurement raster.
#
# The raster sits beside the annotation under the same stem: `annotation/s1a-iw2-slc-hh-….xml` names
# `measurement/s1a-iw2-slc-hh-….tiff`. Only a `.SAFE` directory can be indexed, because a zip stores the
# raster deflated — reaching a line means inflating everything before it, and the strip table it would
# need first is at the end of the entry. Metadata reads from a zip are unaffected, so the refusal
# happens here rather than at `open_slc`.

"""
    measurement_path(product::Sentinel1Product, swath) -> String

The measurement raster of one subswath.

Throws for a zipped product, naming what to do instead: the raster is deflated inside the archive, so
its lines are not addressable without inflating the whole entry.
"""
function measurement_path(p::Sentinel1Product, swath::Integer)
    isdir(p.path) || throw(ArgumentError(
        "`$(p.path)` is a zipped Sentinel-1 product, whose measurement raster is deflated inside " *
        "the archive and so cannot be read a window at a time. Unpack it and open the `.SAFE` " *
        "directory, or read the bursts from ASF's burst extractor. Reading the metadata of a zipped " *
        "product needs no unpacking and is unaffected"))

    dir = joinpath(p.path, "measurement")
    isdir(dir) || throw(ArgumentError(
        "`$(p.path)` has no `measurement` directory, so it carries annotation only"))

    id_str = _swath_id(swath, p.polarization)
    hits = filter(n -> occursin(id_str, n) && endswith(n, ".tiff"), readdir(dir))
    isempty(hits) && throw(ArgumentError(
        "`$(p.path)` has no measurement raster for subswath IW$swath polarization " *
        "$(uppercase(p.polarization))"))
    return joinpath(dir, only(sort!(hits)))
end
