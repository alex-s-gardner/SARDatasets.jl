# The merged image itself: a matrix whose samples are read from the bursts when they are asked for.
#
# A subswath is a few gigabytes and a correlator reads it in windows, so nothing is copied up front. The
# array holds a [`BurstGrid`](@ref) and one raster per burst, and a read resolves each row to the burst
# covering it. Rows and samples outside every burst's valid region read as zero — the merged image is
# wider and taller than the union of the regions placed in it, and the corners it leaves are not data.
#
# Reading a window is what this is built for: `A[rows, cols]` groups the rows by the burst they come
# from and takes one strip from each, where scalar indexing pays a row lookup per sample. That is the
# operation `AutoRIFT.jl`'s blocked path uses.
#
# Zero is what the margin reads as because the correlator treats it as absent, and because a filter over
# a window straddling the edge must not see whatever the processor left in the file. It is not a claim
# that no sample there is zero: [`validmask`](@ref) is what says where the data is, and a caller that
# needs the distinction should use it rather than testing for zero.

"""
    ConcatenatedBursts{T,S} <: AbstractMatrix{T}

The bursts of one Sentinel-1 subswath as a single matrix, read on demand.

Rows come from whichever burst covers them under the array's [`BurstGrid`](@ref); rows and samples no
burst reaches read as `zero(T)`. Build one with [`pixels`](@ref) on a merged [`SLC`](@ref) rather than
calling this constructor.

Index with ranges — `A[rows, cols]` — wherever the shape of the read allows it: a window is taken as
one strip per burst it spans, while scalar indexing resolves a row per sample.
"""
struct ConcatenatedBursts{T,S<:AbstractMatrix{T}} <: AbstractMatrix{T}
    # Parallel to `grid.placements`: `sources[k]` is the raster `grid.placements[k]` reads from. One
    # raster per placement rather than one per product, so bursts delivered as separate files and
    # bursts sharing a subswath's raster are the same array.
    sources::Vector{S}
    # Where in `sources[k]` the burst begins: zero when the raster holds that burst alone, and
    # `(burst - 1) * lines_per_burst` when it holds the whole subswath stacked, which is how a `.SAFE`
    # measurement raster stores it. A placement's `burst_rows` count from the burst's own first line,
    # so this is what turns them into rows of the file.
    row_offsets::Vector{Int}
    grid::BurstGrid
end

function ConcatenatedBursts(sources::AbstractVector{<:AbstractMatrix}, grid::BurstGrid;
                            row_offsets::AbstractVector{<:Integer} = zeros(Int, length(sources)))
    length(sources) == length(grid.placements) || throw(ArgumentError(
        "the grid places $(length(grid.placements)) bursts but $(length(sources)) rasters were " *
        "given; there must be one raster per placed burst"))
    length(row_offsets) == length(sources) || throw(ArgumentError(
        "$(length(sources)) rasters were given with $(length(row_offsets)) row offsets; there must " *
        "be one offset per raster"))
    for (k, p) in enumerate(grid.placements)
        src = sources[k]
        off = row_offsets[k]
        # A raster too small for the rows and samples its placement takes would read outside it.
        (off + last(p.burst_rows) <= size(src, 1) && last(p.cols) <= size(src, 2)) ||
            throw(ArgumentError(
                "burst $(p.burst) contributes rows $(off + first(p.burst_rows))-" *
                "$(off + last(p.burst_rows)) and samples $(p.cols) of a raster that is " *
                "$(size(src, 1))x$(size(src, 2))"))
    end
    src = collect(sources)
    return ConcatenatedBursts{eltype(eltype(src)),eltype(src)}(src, collect(Int, row_offsets), grid)
end

Base.size(A::ConcatenatedBursts) = size(A.grid)
Base.IndexStyle(::Type{<:ConcatenatedBursts}) = IndexCartesian()

"""
    grid(A::ConcatenatedBursts) -> BurstGrid

Where each burst sits in the merged image.
"""
grid(A::ConcatenatedBursts) = A.grid

# The placement covering a row and the index of its raster, or `nothing` outside every burst.
function _source_at(A::ConcatenatedBursts, row::Integer)
    ps = A.grid.placements
    lo, hi = 1, length(ps)
    while lo <= hi
        mid = (lo + hi) >>> 1
        p = ps[mid]
        row < first(p.grid_rows) ? (hi = mid - 1) :
        row > last(p.grid_rows) ? (lo = mid + 1) :
        return (p, mid)
    end
    return nothing
end

function Base.getindex(A::ConcatenatedBursts{T}, i::Int, j::Int) where {T}
    @boundscheck checkbounds(A, i, j)
    hit = _source_at(A, i)
    hit === nothing && return zero(T)
    p, k = hit
    j in p.cols || return zero(T)
    row = A.row_offsets[k] + first(p.burst_rows) + (i - first(p.grid_rows))
    return A.sources[k][row, j]
end

function Base.getindex(A::ConcatenatedBursts{T}, rows::AbstractUnitRange{<:Integer},
                       cols::AbstractUnitRange{<:Integer}) where {T}
    @boundscheck checkbounds(A, rows, cols)
    out = Matrix{T}(undef, length(rows), length(cols))
    _read_window!(out, A, rows, cols)
    return out
end

# One strip per burst the window spans, taken from the burst's raster in a single read. The window is
# zeroed first, so the rows and samples no burst covers need no separate pass.
function _read_window!(out, A::ConcatenatedBursts{T}, rows, cols) where {T}
    fill!(out, zero(T))
    (isempty(rows) || isempty(cols)) && return out

    for (k, p) in enumerate(A.grid.placements)
        # The part of this burst's placement the window asks for.
        lo = max(first(rows), first(p.grid_rows))
        hi = min(last(rows), last(p.grid_rows))
        lo <= hi || continue
        cl = max(first(cols), first(p.cols))
        ch = min(last(cols), last(p.cols))
        cl <= ch || continue

        shift = A.row_offsets[k] + first(p.burst_rows) - first(p.grid_rows)
        strip = A.sources[k][(lo + shift):(hi + shift), cl:ch]
        dest = view(out, (lo - first(rows) + 1):(hi - first(rows) + 1),
                    (cl - first(cols) + 1):(ch - first(cols) + 1))
        copyto!(dest, strip)
    end
    return out
end

"""
    BurstValidMask <: AbstractMatrix{Bool}

Which samples of a merged image were imaged, derived from the grid rather than stored.

`true` exactly where [`ConcatenatedBursts`](@ref) reads from a burst's valid region. The image reads
zero outside, but a sample inside can be zero too, so this is what distinguishes absent from dark and
what a correlator should be given rather than a test against zero.
"""
struct BurstValidMask <: AbstractMatrix{Bool}
    grid::BurstGrid
end

Base.size(m::BurstValidMask) = size(m.grid)
Base.IndexStyle(::Type{BurstValidMask}) = IndexCartesian()

function Base.getindex(m::BurstValidMask, i::Int, j::Int)
    @boundscheck checkbounds(m, i, j)
    p = _placement_at(m.grid, i)
    return p !== nothing && j in p.cols
end

function Base.getindex(m::BurstValidMask, rows::AbstractUnitRange{<:Integer},
                       cols::AbstractUnitRange{<:Integer})
    @boundscheck checkbounds(m, rows, cols)
    out = fill(false, length(rows), length(cols))
    (isempty(rows) || isempty(cols)) && return out
    for p in m.grid.placements
        lo = max(first(rows), first(p.grid_rows))
        hi = min(last(rows), last(p.grid_rows))
        lo <= hi || continue
        cl = max(first(cols), first(p.cols))
        ch = min(last(cols), last(p.cols))
        cl <= ch || continue
        out[(lo - first(rows) + 1):(hi - first(rows) + 1),
            (cl - first(cols) + 1):(ch - first(cols) + 1)] .= true
    end
    return out
end

"""
    validmask(A) -> AbstractMatrix{Bool}

Which samples of `A` were imaged.

For a merged subswath this is [`BurstValidMask`](@ref), computed from the burst placements and holding
nothing, so asking for it costs no read of the rasters.
"""
validmask(A::ConcatenatedBursts) = BurstValidMask(A.grid)
