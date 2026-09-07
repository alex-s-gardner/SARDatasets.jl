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
struct MergedBurstBackend{B<:AbstractSLCBackend} <: AbstractSLCBackend
    # The backends of the bursts merged, in azimuth order and parallel to `grid.placements`. Holding
    # them rather than a product is what lets bursts of a `.SAFE` and bursts delivered one file each
    # merge by the same path: each knows how to reach its own samples.
    sources::Vector{B}
    annotation::SubswathAnnotation
    swath::Int
    bursts::UnitRange{Int}
    grid::BurstGrid
    orbit_path::String
end

_path(b::MergedBurstBackend) = _path(first(b.sources))
_leading_annotation(b::MergedBurstBackend) = b.annotation

# The epoch the merged times are reported against: the first burst placed, as for any single burst. So a
# merge of one burst reports exactly what that burst reports.
_anchor(b::MergedBurstBackend) = b.annotation.burst_start[first(b.bursts)]

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
    swath, range = _merge_extent(backends)
    a = _leading_annotation(first(backends))
    _check_merge_agreement(slcs, a, swath)

    grid = burst_grid(a, range; tolerance)
    return SLC(MergedBurstBackend(backends, a, swath, range, grid,
                                  _merge_orbit_path(backends)))
end

# The backend of one burst, refused unless it is one: a mosaic has no subswath to place bursts in, and
# another sensor's product is not made of bursts at all.
function _merge_backend(s::AbstractSLC, i::Integer)
    b = s.backend
    if b isa Sentinel1Backend
        b.swath === nothing && throw(ArgumentError(
            "acquisition $i describes the mosaic across `$(_path(b))`'s subswaths rather than a " *
            "burst of one. Merge the bursts of a single subswath, which `bursts(path; swath)` " *
            "returns"))
        return b
    end
    b isa AsfBurstBackend && return b
    throw(ArgumentError(
        "acquisition $i is a $(nameof(typeof(b))) rather than a Sentinel-1 burst, and only " *
        "Sentinel-1 bursts merge"))
end

# Which burst of which subswath a backend describes, whatever delivered it, and which product it
# describes it from. `_burst_source` is what a burst must agree on for the merge to be placing bursts of
# one acquisition: the granule for an ASF burst, the container for a `.SAFE` one.
_burst_index(b::Sentinel1Backend) = b.burst::Int
_burst_index(b::AsfBurstBackend) = b.source.burst
_burst_swath(b::Sentinel1Backend) = b.swath::Int
_burst_swath(b::AsfBurstBackend) = b.source.swath
_orbit_path(b::Sentinel1Backend) = b.product.orbit_path
_orbit_path(b::AsfBurstBackend) = b.orbit_path
_burst_source(b::Sentinel1Backend) = b.product.path
_burst_source(b::AsfBurstBackend) = b.source.slc
_burst_polarization(b::Sentinel1Backend) = lowercase(b.product.polarization)
_burst_polarization(b::AsfBurstBackend) = lowercase(b.source.polarization)

# The one subswath and consecutive burst range the given bursts amount to.
function _merge_extent(backends::AbstractVector)
    swath = _burst_swath(first(backends))
    same = typeof(first(backends))
    source = _burst_source(first(backends))
    polarization = _burst_polarization(first(backends))
    for (i, b) in enumerate(backends)
        b isa same || throw(ArgumentError(
            "burst $i is a $(nameof(typeof(b))) but the first a $(nameof(same)); bursts reached " *
            "different ways are not known to come from one product"))
        # Consecutive burst numbers of two products are not consecutive bursts of anything. Two slices
        # of one datatake share an absolute orbit and a range origin, so nothing further downstream
        # would catch it: the grid would place one product's bursts by the other's annotation.
        _burst_source(b) == source || throw(ArgumentError(
            "burst $i comes from `$(_burst_source(b))` but the first from `$source`; bursts of " *
            "different products describe different acquisitions and do not merge"))
        _burst_swath(b) == swath || throw(ArgumentError(
            "burst $i is in subswath IW$(_burst_swath(b)) but the first in IW$swath. Subswaths " *
            "lie at different slant ranges, so merging them would need one range origin for two, " *
            "and their samples are not on one grid"))
        # Two channels of one subswath share their mission, orbit and range geometry, so this is the
        # only check that separates them; merging them would interleave two channels' samples.
        _burst_polarization(b) == polarization || throw(ArgumentError(
            "burst $i is polarization $(uppercase(_burst_polarization(b))) but the first " *
            "$(uppercase(polarization)); a merged subswath carries one channel"))
    end

    idx = [_burst_index(b) for b in backends]
    expected = first(idx):(first(idx) + length(idx) - 1)
    idx == collect(expected) || throw(ArgumentError(
        "the bursts given are $idx, which is not a consecutive ascending run. A merged subswath " *
        "spans one unbroken range of bursts, since a gap would leave rows no burst covers and a " *
        "repeat would place one twice"))
    return swath, expected
end

function _merge_orbit_path(backends::AbstractVector)
    path = _orbit_path(first(backends))
    for (i, b) in enumerate(backends)
        _orbit_path(b) == path || throw(ArgumentError(
            "burst $i takes its state vectors from `$(_orbit_path(b))` but the first from " *
            "`$path`; one acquisition has one orbit"))
    end
    return path
end

# What the bursts report about themselves, against the annotation the grid is built from. `_merge_extent`
# has already established they name one product, subswath and channel; this catches an acquisition whose
# own geometry does not match that annotation, which a caller assembling backends by hand could produce.
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
    orbit_path = b.orbit_path

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

# How the samples are reached depends on how the bursts were delivered, so the merged backend defers to
# its sources: a `.SAFE` holds one raster per subswath with every burst stacked in it, while ASF
# delivers a file per burst.
read_pixels(b::MergedBurstBackend) = _merged_pixels(b, b.sources)

# One memory-mapped file, read at a burst-dependent row offset.
function _merged_pixels(b::MergedBurstBackend, sources::AbstractVector{Sentinel1Backend})
    a = _leading_annotation(b)
    product = first(sources).product
    raster = open_tiff(measurement_path(product, b.swath))

    expected = a.lines_per_burst * nbursts(a)
    size(raster, 1) == expected || throw(ArgumentError(
        "`$(path(raster))` has $(size(raster, 1)) lines, but subswath IW$(b.swath) has " *
        "$(nbursts(a)) bursts of $(a.lines_per_burst) lines and so should have $expected. The " *
        "raster and the annotation describe different products"))
    size(raster, 2) == a.samples_per_burst || throw(ArgumentError(
        "`$(path(raster))` is $(size(raster, 2)) samples wide but subswath IW$(b.swath) records " *
        "$(a.samples_per_burst)"))

    rasters = [raster for _ in b.grid.placements]
    offsets = [(p.burst - 1) * a.lines_per_burst for p in b.grid.placements]
    return ConcatenatedBursts(rasters, b.grid; row_offsets = offsets)
end

# One file per burst, each holding that burst alone, so every placement reads at row offset zero.
_merged_pixels(b::MergedBurstBackend, sources::AbstractVector{<:AsfBurstBackend}) =
    ConcatenatedBursts([read_pixels(s) for s in sources], b.grid)

"""
    grid(s::SLC) -> BurstGrid

Where each burst sits in a merged acquisition.

Throws for an acquisition that is not a merge of bursts.
"""
grid(s::AbstractSLC) = grid(s.backend)
grid(b::MergedBurstBackend) = b.grid
grid(b::AbstractSLCBackend) = throw(ArgumentError(
    "a $(nameof(typeof(b))) is not a merge of bursts, so it has no burst grid"))
