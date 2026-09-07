# Folding a subswath's bursts into one acquisition.
#
# The bursts of a subswath are already `SLC`s, so the merge takes them rather than a product and a range:
# whatever vector a caller has — a whole subswath, a slice of one, the result of a `filter` — is the
# selection. What comes back is an `SLC` like any other, so everything downstream of this package reads
# it without knowing it was merged.
#
# The work that is not layout is reconciling the clocks. Each burst reports its times against an epoch
# derived from its own start, so a subswath of nine bursts carries nine different epochs and their
# `sensing_start` values are not on one scale — differencing them across a subswath is out by the whole
# spread of the epochs, some tens of seconds. The merged acquisition adopts the first burst's epoch and
# re-references every other burst onto it, which is the same reconciliation `repeat_interval` performs
# for a pair of products.

"""
    MergedBurstBackend <: AbstractSLCBackend

Several bursts of one Sentinel-1 subswath, described as a single acquisition.

Built by [`merge_bursts`](@ref). Holds the parsed product, the subswath, and the layout the bursts were
placed with, so the geometry and the samples both derive from one grid.
"""
struct MergedBurstBackend <: AbstractSLCBackend
    product::Sentinel1Product
    swath::Int
    bursts::UnitRange{Int}
    grid::BurstGrid
end

_path(b::MergedBurstBackend) = b.product.path
_leading_annotation(b::MergedBurstBackend) = annotation(b.product, b.swath)

# The epoch the merged times are reported against: the first burst placed, as for any single burst. So a
# merge of one burst reports exactly what that burst reports.
_anchor(b::MergedBurstBackend) = _leading_annotation(b).burst_start[first(b.bursts)]

"""
    merge_bursts(bursts::AbstractVector{<:SLC}; tolerance = MAX_BURST_GRID_RESIDUAL) -> SLC

Place consecutive bursts of one subswath on a single uniform azimuth grid.

The bursts must come from one product, subswath and polarization, and be consecutive and in order —
which is what [`bursts`](@ref) returns, so a whole subswath or any slice of one merges directly. The
result is an [`SLC`](@ref) whose `nlines` spans from the first burst's first line to the last burst's
last, with the overlap between neighbours divided at its midpoint rather than repeated.

A line's azimuth time is `sensing_start + line / prf` throughout the result, which is what
[`RadarGeometry`](@ref) means by `prf` and what a stack of whole bursts would break at every seam. The
samples are placed rather than resampled: each burst's start lands on the shared grid within a fraction
of a line, and `tolerance` is how far it may be off before the merge is refused.

Use [`pixels`](@ref) or [`amplitude`](@ref) for the samples, which are read as they are indexed, and
[`validmask`](@ref) for which of them were imaged.

# Examples

```julia
b = bursts(safe; orbit = eof, swath = 2)
ref = merge_bursts(b)          # the whole subswath
ref = merge_bursts(b[3:7])     # five of its bursts
nlines(ref)
amplitude(ref)[1:512, 1:512]
```
"""
function merge_bursts(slcs::AbstractVector{<:AbstractSLC};
                      tolerance::Real = MAX_BURST_GRID_RESIDUAL)
    isempty(slcs) && throw(ArgumentError(
        "a merged acquisition needs at least one burst, but no bursts were given"))

    backends = [_merge_backend(s, i) for (i, s) in enumerate(slcs)]
    product, swath, range = _merge_extent(backends)
    a = annotation(product, swath)
    _check_merge_agreement(slcs, a, swath)

    grid = burst_grid(a, range; tolerance)
    return SLC(MergedBurstBackend(product, swath, range, grid))
end

# The backend of one burst, refused unless it is a burst of a Sentinel-1 product: a mosaic has no
# subswath to place bursts in, and another sensor's product is not made of bursts at all.
function _merge_backend(s::AbstractSLC, i::Integer)
    b = s.backend
    b isa Sentinel1Backend || throw(ArgumentError(
        "acquisition $i is a $(nameof(typeof(b))) rather than a Sentinel-1 burst, and only " *
        "Sentinel-1 bursts merge"))
    b.swath === nothing && throw(ArgumentError(
        "acquisition $i describes the mosaic across `$(_path(b))`'s subswaths rather than a burst " *
        "of one. Merge the bursts of a single subswath, which `bursts(path; swath)` returns"))
    return b
end

# The one product, subswath and consecutive burst range the given bursts amount to.
function _merge_extent(backends::AbstractVector{Sentinel1Backend})
    first_backend = first(backends)
    product = first_backend.product
    swath = first_backend.swath::Int

    for (i, b) in enumerate(backends)
        b.product === product || b.product.path == product.path || throw(ArgumentError(
            "burst $i comes from `$(_path(b))` but the first from `$(product.path)`; bursts of " *
            "different products describe different acquisitions and do not merge"))
        b.swath == swath || throw(ArgumentError(
            "burst $i is in subswath IW$(b.swath) but the first in IW$swath. Subswaths lie at " *
            "different slant ranges, so merging them would need one range origin for two, and " *
            "their samples are not on one grid"))
    end

    idx = [b.burst::Int for b in backends]
    expected = first(idx):(first(idx) + length(idx) - 1)
    idx == collect(expected) || throw(ArgumentError(
        "the bursts given are $idx, which is not a consecutive ascending run. A merged subswath " *
        "spans one unbroken range of bursts, since a gap would leave rows no burst covers and a " *
        "repeat would place one twice"))
    return product, swath, expected
end

# The bursts must agree on what they are, not merely on where they came from. A caller reaching a
# `Sentinel1Product` directly could build bursts from one product with differing polarizations, and the
# geometry would then be a blend of two channels.
function _check_merge_agreement(slcs, a::SubswathAnnotation, swath::Integer)
    ref = first(slcs).identification
    for (i, s) in enumerate(slcs)
        id = s.identification
        id.mission == ref.mission || throw(ArgumentError(
            "burst $i is from $(id.mission) but the first from $(ref.mission)"))
        id.absolute_orbit == ref.absolute_orbit || throw(ArgumentError(
            "burst $i is on absolute orbit $(id.absolute_orbit) but the first on " *
            "$(ref.absolute_orbit); bursts of different orbits are different acquisitions"))
        # Every burst of a subswath shares its range geometry, so a disagreement here means the
        # annotation they were built from is not the one being placed.
        s.geometry.starting_range == a.starting_range || throw(ArgumentError(
            "burst $i starts at range $(s.geometry.starting_range) m but subswath IW$swath at " *
            "$(a.starting_range) m"))
    end
    return nothing
end

# The merged image spans the grid, and its azimuth window runs from the anchor burst's first line to the
# grid's last row — one line per `1/prf` throughout, which is what makes the window and the row count
# consistent with each other.
function read_geometry(b::MergedBurstBackend)
    a = _leading_annotation(b)
    grid = b.grid

    prf = 1 / a.azimuth_time_interval
    start = _anchor(b)
    stop = _advance(start, (grid.nlines - 1) * a.azimuth_time_interval)
    far_range = a.starting_range + (grid.nsamples - 1.0) * a.range_pixel_spacing

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
        grid.nlines,
        grid.nsamples,
        S1_LOOK_SIDE,
        epoch,
    )
end

function read_identification(b::MergedBurstBackend)
    a = _leading_annotation(b)
    start = _anchor(b)
    stop = _advance(start, (b.grid.nlines - 1) * a.azimuth_time_interval)
    return Identification(
        a.mission,
        a.product_type,
        a.absolute_orbit,
        lowercase(a.pass_direction),
        S1_LOOK_SIDE == LookLeft ? "Left" : "Right",
        _utc_string(start),
        _utc_string(stop),
        # As for a burst or a mosaic: the footprint is in `manifest.safe`, and the whole product's
        # polygon would not describe this subswath's merged extent.
        "",
    )
end

# The orbit must cover the merged window rather than one burst's, since a solve anywhere in the image
# has to interpolate rather than extrapolate. So the window widened here is the whole grid's.
function read_orbit(b::MergedBurstBackend)
    a = _leading_annotation(b)
    start = _anchor(b)
    stop = _advance(start, (b.grid.nlines - 1) * a.azimuth_time_interval)
    orbit_path = b.product.orbit_path

    table = read_eof_state_vectors(orbit_path; from = start, to = stop,
                                   padding = S1_ORBIT_PADDING)
    isempty(table.time) && throw(ArgumentError(
        "`$orbit_path` has no state vectors within $(S1_ORBIT_PADDING) s of the merged window of " *
        "`$(_path(b))` subswath IW$(b.swath); it is probably the orbit file of a different granule"))

    epoch = UtcTime(epoch_of(start), 0.0)
    return StateVectors(
        [seconds_between(epoch, t) for t in table.time],
        table.position,
        table.velocity,
        epoch.datetime,
        "Hermite",
        _eof_kind(orbit_path),
    )
end

# A `.SAFE` holds one raster per subswath with every burst stacked in it, so all the placements read from
# one memory-mapped file at a burst-dependent row offset.
function read_pixels(b::MergedBurstBackend)
    a = _leading_annotation(b)
    raster = open_tiff(measurement_path(b.product, b.swath))

    expected = a.lines_per_burst * nbursts(a)
    size(raster, 1) == expected || throw(ArgumentError(
        "`$(path(raster))` has $(size(raster, 1)) lines, but subswath IW$(b.swath) has " *
        "$(nbursts(a)) bursts of $(a.lines_per_burst) lines and so should have $expected. The " *
        "raster and the annotation describe different products"))
    size(raster, 2) == a.samples_per_burst || throw(ArgumentError(
        "`$(path(raster))` is $(size(raster, 2)) samples wide but subswath IW$(b.swath) records " *
        "$(a.samples_per_burst)"))

    sources = [raster for _ in b.grid.placements]
    offsets = [(p.burst - 1) * a.lines_per_burst for p in b.grid.placements]
    return ConcatenatedBursts(sources, b.grid; row_offsets = offsets)
end

"""
    grid(s::SLC) -> BurstGrid

Where each burst sits in a merged acquisition.

Throws for an acquisition that is not a merge of bursts.
"""
grid(s::AbstractSLC) = grid(s.backend)
grid(b::MergedBurstBackend) = b.grid
grid(b::AbstractSLCBackend) = throw(ArgumentError(
    "a $(nameof(typeof(b))) is not a merge of bursts, so it has no burst grid"))
