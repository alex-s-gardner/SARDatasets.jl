# Where each burst of a subswath belongs in one image.
#
# A TOPS subswath is a run of bursts overlapping in azimuth by about a tenth of their length, each
# recording the time of its own first line. Two images can be built from them.
#
# One stacks every burst's full extent, which is what a `.SAFE` measurement raster holds: its height is
# `lines_per_burst` times the burst count, the overlap appears twice, and the azimuth time of a line is
# not `sensing_start + line / prf` — the relation jumps at every seam. That image cannot be described by
# a `RadarGeometry`, whose `prf` is exactly the claim that the relation holds.
#
# The other places each burst on one uniform azimuth time axis and splits the overlap between the two
# bursts that cover it, which is what this builds. The line-to-time relation then holds everywhere, so
# the image is one an acquisition can be reported for. It is also exact rather than resampled: a burst's
# start lands on the shared axis to a fraction of a line — measured across the granules this package is
# tested against, within 2.2e-4 of one — so placing it is an integer row offset and a copy. The residual
# is checked rather than assumed, since a product where it did not hold would need interpolation this
# does not do.
#
# Only a burst's valid region is placed. The margin outside it is not zero in the file, so a merge that
# copied whole bursts would admit samples the processor never imaged.

"""
    BurstPlacement

Where one burst's valid region lands in a merged image.

`grid_rows` and `burst_rows` are the same length: row `grid_rows[k]` of the image is row
`burst_rows[k]` of the burst. `cols` indexes both, since bursts of one subswath share a range origin.
"""
struct BurstPlacement
    burst::Int
    grid_rows::UnitRange{Int}
    burst_rows::UnitRange{Int}
    cols::UnitRange{Int}
end

"""
    BurstGrid

The shape of a merged subswath and the placement of every burst in it.

Built by [`burst_grid`](@ref). Holds no pixels and no times: it is the integer layout alone, so it can
be checked against the arithmetic it reproduces without reading a raster.

# Fields
- `nlines`, `nsamples`: the merged image's size.
- `placements`: one per burst merged, in azimuth order, with disjoint `grid_rows`.
- `residual`: how far the worst burst start fell from the shared azimuth axis, in lines. A measure of
  how well the uniform grid describes this product; a copy is exact only because it is small.
"""
struct BurstGrid
    nlines::Int
    nsamples::Int
    placements::Vector{BurstPlacement}
    residual::Float64
end

Base.size(g::BurstGrid) = (g.nlines, g.nsamples)

"""
    MAX_BURST_GRID_RESIDUAL

How far a burst's start may fall from the uniform azimuth axis before a merge is refused, in lines.

A burst is placed at a whole-row offset, so a start lying between rows would shift its samples in
azimuth by the difference. Every granule measured places within 2.2e-4 lines; this bound is two orders
of magnitude above that and still far below the half-line where a placement would round the other way.
"""
const MAX_BURST_GRID_RESIDUAL = 1.0e-2

"""
    burst_grid(a::SubswathAnnotation, bursts = 1:nbursts(a); tolerance = MAX_BURST_GRID_RESIDUAL)
      -> BurstGrid

Lay `bursts` of subswath `a` out on one uniform azimuth time axis.

`bursts` must be consecutive and ascending. Row 1 of the grid is the first burst's first line, whether
or not that line is valid, so the grid's azimuth origin is that burst's own recorded time.
"""
function burst_grid(a::SubswathAnnotation, bursts::AbstractUnitRange{<:Integer} = 1:nbursts(a);
                    tolerance::Real = MAX_BURST_GRID_RESIDUAL)
    n = nbursts(a)
    (1 <= first(bursts) && last(bursts) <= n) || throw(ArgumentError(
        "subswath IW$(a.swath) has $n bursts, so $bursts is not a range of them"))
    isempty(bursts) && throw(ArgumentError("a merged subswath needs at least one burst"))

    idx = collect(bursts)
    dt = a.azimuth_time_interval
    anchor = a.burst_start[first(idx)]

    # Each burst's first line, as a row offset on the axis anchored at the first burst's first line.
    # The offset must be a whole number of rows: it is what makes placement a copy.
    offsets = Vector{Int}(undef, length(idx))
    residual = 0.0
    worst = first(idx)
    for (k, b) in enumerate(idx)
        exact = seconds_between(anchor, a.burst_start[b]) / dt
        offsets[k] = round(Int, exact)
        off_grid = abs(exact - offsets[k])
        if off_grid > residual
            residual = off_grid
            worst = b
        end
    end
    residual <= tolerance || throw(ArgumentError(
        "burst $worst of subswath IW$(a.swath) starts $residual lines from the uniform azimuth " *
        "grid, more than the $tolerance allowed. Placing it at a whole row shifts its samples in " *
        "azimuth by that much, and this reader places bursts rather than resampling them"))

    # The rows of the grid each burst's valid region would occupy before the overlap is divided, as a
    # half-open span: `starts[k]` up to but excluding `stops[k]`.
    starts = [offsets[k] + a.first_valid_line[b] - 1 for (k, b) in enumerate(idx)]
    stops = [starts[k] + length(valid_lines(a, b)) for (k, b) in enumerate(idx)]

    # Each burst's span before the seams are settled, as half-open grid rows and the burst row its
    # first grid row comes from.
    span_lo = Vector{Int}(undef, length(idx))
    span_hi = Vector{Int}(undef, length(idx))
    row_lo = Vector{Int}(undef, length(idx))
    for (k, b) in enumerate(idx)
        # Half of the overlap with each neighbour, so the seam sits at its midpoint. `fld` rather than
        # `÷`: bursts that fall short of overlapping make the count negative, and the two disagree
        # there.
        head = k == 1 ? 0 : fld(stops[k - 1] - starts[k], 2)
        tail = k == length(idx) ? 0 : fld(stops[k] - starts[k + 1], 2)
        span_lo[k] = starts[k] + head
        span_hi[k] = stops[k] - tail
        row_lo[k] = a.first_valid_line[b] + head
    end

    # An overlap of an odd number of rows leaves one row over, which both neighbours' halves reach.
    # It belongs to the later burst — that is the row the reference leaves there, having written the
    # bursts in order — so the earlier one gives it up and the placements stay disjoint.
    for k in 1:(length(idx) - 1)
        span_hi[k] = min(span_hi[k], span_lo[k + 1])
    end

    placements = Vector{BurstPlacement}(undef, length(idx))
    for (k, b) in enumerate(idx)
        rows = span_hi[k] - span_lo[k]
        rows > 0 || throw(ArgumentError(
            "burst $b of subswath IW$(a.swath) is wholly covered by its neighbours' halves of the " *
            "overlap, so it contributes nothing to the merged image"))
        placements[k] = BurstPlacement(b,
                                       (span_lo[k] + 1):span_hi[k],
                                       row_lo[k]:(row_lo[k] + rows - 1),
                                       valid_samples(a, b))
        last(placements[k].burst_rows) <= a.last_valid_line[b] || throw(ArgumentError(
            "burst $b of subswath IW$(a.swath) would contribute rows past its valid region"))
    end

    # The axis runs from the first burst's first line to the last burst's last line, valid or not, so
    # that its length is `1/prf` per row throughout.
    nlines = offsets[end] + a.lines_per_burst
    _check_disjoint(placements, nlines, a.swath)
    return BurstGrid(nlines, a.samples_per_burst, placements, residual)
end

# Two bursts writing the same row would make the later one silently win, and a row inside the span no
# burst reaches would read as zero without saying so. Neither can happen for consecutive bursts of a
# real product, which is why finding one means the placement arithmetic is wrong rather than the
# product unusual.
function _check_disjoint(placements, nlines, swath)
    prev = placements[1]
    (first(prev.grid_rows) >= 1 && last(placements[end].grid_rows) <= nlines) || throw(ArgumentError(
        "subswath IW$swath places bursts in rows " *
        "$(first(prev.grid_rows))-$(last(placements[end].grid_rows)) of a $nlines-row image"))
    for p in @view placements[2:end]
        first(p.grid_rows) == last(prev.grid_rows) + 1 || throw(ArgumentError(
            "burst $(prev.burst) of subswath IW$swath ends at row $(last(prev.grid_rows)) and " *
            "burst $(p.burst) begins at row $(first(p.grid_rows)), so the merged image would " *
            (first(p.grid_rows) > last(prev.grid_rows) + 1 ? "have a gap between them" :
             "have them overlap")))
        prev = p
    end
    return nothing
end

# The placement covering a row, or `nothing` where no burst reaches. Placements are ascending and
# disjoint, so this is a search rather than a scan.
function _placement_at(g::BurstGrid, row::Integer)
    lo, hi = 1, length(g.placements)
    while lo <= hi
        mid = (lo + hi) >>> 1
        p = g.placements[mid]
        row < first(p.grid_rows) ? (hi = mid - 1) :
        row > last(p.grid_rows) ? (lo = mid + 1) :
        return p
    end
    return nothing
end
