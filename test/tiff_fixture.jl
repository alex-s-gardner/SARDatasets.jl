# Writing a striped TIFF for the reader to read.
#
# The layout a Sentinel-1 measurement raster has is small enough to write directly: a header, a
# directory of the tags the reader checks, and one strip per line. So the reader is tested against
# files built here rather than against a granule, and the departures it must refuse — compression, a
# tiled raster, a BigTIFF, a strip count disagreeing with the line count — are each written on
# purpose.

const TIFF_TAGS = (
    width = 256, length = 257, bits = 258, compression = 259, photometric = 262,
    strip_offsets = 273, samples = 277, rows_per_strip = 278, strip_bytes = 279,
    planar = 284, tile_width = 322, tile_length = 323, sample_format = 339,
)

const TIFF_SHORT = 3
const TIFF_LONG = 4

"""
    write_tiff(path, pixels; kwargs...)

Write `pixels` as a striped TIFF and return `path`.

Each keyword overrides one directory value, so a file that departs from the layout the reader accepts
can be written without a second writer: `compression = 5` writes a file claiming LZW,
`rows_per_strip = 8` one claiming eight lines to a strip, `magic = 43` a BigTIFF version stamp, and
`nstrips` a strip table of the wrong length. `bigendian` writes the whole file `MM`.
"""
function write_tiff(path::AbstractString, pixels::AbstractMatrix{Complex{Int16}};
                    compression::Integer = 1, rows_per_strip::Integer = 1,
                    samples_per_pixel::Integer = 1, bits_per_sample::Integer = 32,
                    sample_format::Integer = 5, planar::Union{Nothing,Integer} = 1,
                    magic::Integer = 42, bigendian::Bool = false,
                    tiled::Bool = false, nstrips::Union{Nothing,Integer} = nothing,
                    strip_bytes::Union{Nothing,Integer} = nothing,
                    pad::Integer = 0)
    nlines, nsamples = size(pixels)
    conv = bigendian ? bswap : identity

    entries = Pair{Int,Tuple{Int,Int,UInt32}}[]      # tag => (type, count, payload)
    put!(tag, typ, count, payload) = push!(entries, tag => (typ, count, UInt32(payload)))

    # The strip table and the raster follow the directory, so their offsets depend on the directory's
    # size. Counting the entries first is what lets the offsets be written in one pass.
    nentries = 10 + (planar === nothing ? 0 : 1) + (tiled ? 2 : 0)
    header = 8
    dirsize = 2 + 12 * nentries + 4
    table_count = nstrips === nothing ? nlines : Int(nstrips)
    offsets_at = header + dirsize
    counts_at = offsets_at + 4 * table_count
    # `pad` shifts the raster off a multiple of the sample size, which is where a real writer leaves it.
    # A reader loading a sample as a word rather than assembling its bytes fails only on such a file.
    raster_at = counts_at + 4 * table_count + Int(pad)

    linebytes = strip_bytes === nothing ? 4 * nsamples : Int(strip_bytes)
    offsets = UInt32[raster_at + (i - 1) * 4 * nsamples for i in 1:table_count]
    counts = UInt32[linebytes for _ in 1:table_count]

    put!(TIFF_TAGS.width, TIFF_SHORT, 1, nsamples)
    put!(TIFF_TAGS.length, TIFF_SHORT, 1, nlines)
    put!(TIFF_TAGS.bits, TIFF_SHORT, 1, bits_per_sample)
    put!(TIFF_TAGS.compression, TIFF_SHORT, 1, compression)
    put!(TIFF_TAGS.photometric, TIFF_SHORT, 1, 1)
    put!(TIFF_TAGS.strip_offsets, TIFF_LONG, table_count, offsets_at)
    put!(TIFF_TAGS.samples, TIFF_SHORT, 1, samples_per_pixel)
    put!(TIFF_TAGS.rows_per_strip, TIFF_SHORT, 1, rows_per_strip)
    put!(TIFF_TAGS.strip_bytes, TIFF_LONG, table_count, counts_at)
    planar === nothing || put!(TIFF_TAGS.planar, TIFF_SHORT, 1, planar)
    if tiled
        put!(TIFF_TAGS.tile_width, TIFF_SHORT, 1, 256)
        put!(TIFF_TAGS.tile_length, TIFF_SHORT, 1, 256)
    end
    put!(TIFF_TAGS.sample_format, TIFF_SHORT, 1, sample_format)

    sort!(entries; by = first)   # a TIFF directory is ordered by tag

    open(path, "w") do io
        write(io, bigendian ? UInt8[0x4d, 0x4d] : UInt8[0x49, 0x49])
        write(io, conv(UInt16(magic)))
        write(io, conv(UInt32(header)))

        write(io, conv(UInt16(length(entries))))
        for (tag, (typ, count, payload)) in entries
            write(io, conv(UInt16(tag)))
            write(io, conv(UInt16(typ)))
            write(io, conv(UInt32(count)))
            # A single SHORT occupies the first two of the payload's four bytes, so on a big-endian
            # file it is written into the high half rather than the low.
            if typ == TIFF_SHORT && count == 1
                bigendian ? write(io, conv(UInt16(payload)), UInt16(0)) :
                            write(io, UInt16(payload), UInt16(0))
            else
                write(io, conv(UInt32(payload)))
            end
        end
        write(io, conv(UInt32(0)))   # no further directory

        for v in offsets; write(io, conv(v)); end
        for v in counts;  write(io, conv(v)); end
        pad > 0 && write(io, zeros(UInt8, Int(pad)))

        # Row-major on disk: line 1's samples, then line 2's.
        for i in 1:nlines, j in 1:nsamples
            p = pixels[i, j]
            write(io, conv(real(p)), conv(imag(p)))
        end
    end
    return path
end

"""
    tiff_pattern(nlines, nsamples) -> Matrix{Complex{Int16}}

A raster whose every sample is a distinct, position-derived value.

Both parts vary with both indices, so a reader that transposed the raster, dropped a line, or read a
neighbouring sample produces a value that cannot be confused with the right one.
"""
function tiff_pattern(nlines::Integer, nsamples::Integer)
    out = Matrix{Complex{Int16}}(undef, nlines, nsamples)
    for i in 1:nlines, j in 1:nsamples
        out[i, j] = Complex{Int16}(Int16(i * 7 + j), Int16(-(i + j * 3)))
    end
    return out
end
