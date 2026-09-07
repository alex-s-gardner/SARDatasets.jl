using SLCDatasets: open_tiff, StripedTiff

@testset "reads a striped raster" begin
    mktempdir() do dir
        pixels = tiff_pattern(9, 13)
        t = open_tiff(write_tiff(joinpath(dir, "p.tiff"), pixels))

        @test size(t) == (9, 13)
        @test eltype(t) == Complex{Int16}
        @test t == pixels

        # Scalar and window indexing must agree, since they take different paths to the same bytes.
        @test t[3, 4] == pixels[3, 4]
        @test t[2:5, 3:9] == pixels[2:5, 3:9]
        @test t[1:1, 1:1] == pixels[1:1, 1:1]
        @test t[9:9, 13:13] == pixels[9:9, 13:13]
        @test collect(t) == pixels
    end
end

@testset "reads a big-endian raster" begin
    mktempdir() do dir
        pixels = tiff_pattern(6, 7)
        t = open_tiff(write_tiff(joinpath(dir, "be.tiff"), pixels; bigendian = true))
        @test size(t) == (6, 7)
        @test t == pixels
        @test t[2:4, 2:6] == pixels[2:4, 2:6]
    end
end

@testset "indexing outside the raster is caught" begin
    mktempdir() do dir
        t = open_tiff(write_tiff(joinpath(dir, "p.tiff"), tiff_pattern(4, 5)))
        @test_throws BoundsError t[5, 1]
        @test_throws BoundsError t[1, 6]
        @test_throws BoundsError t[0, 1]
        @test_throws BoundsError t[1:5, 1:5]
        @test_throws BoundsError t[1:4, 1:6]
    end
end

# Each of these is a layout whose per-line offsets either do not exist or do not mean what the reader
# would take them to mean. Reading one anyway returns plausible garbage, so each must be refused, and
# the message must name the property rather than the byte that failed.
@testset "a raster laid out otherwise is refused" begin
    mktempdir() do dir
        pixels = tiff_pattern(4, 5)

        @test_throws "is compressed" open_tiff(
            write_tiff(joinpath(dir, "lzw.tiff"), pixels; compression = 5))

        @test_throws "one line per strip" open_tiff(
            write_tiff(joinpath(dir, "rows.tiff"), pixels; rows_per_strip = 2))

        @test_throws "tiles rather than strips" open_tiff(
            write_tiff(joinpath(dir, "tiled.tiff"), pixels; tiled = true))

        @test_throws "BigTIFF" open_tiff(
            write_tiff(joinpath(dir, "big.tiff"), pixels; magic = 43))

        @test_throws "one sample per pixel" open_tiff(
            write_tiff(joinpath(dir, "bands.tiff"), pixels; samples_per_pixel = 3))

        @test_throws "complex-integer samples" open_tiff(
            write_tiff(joinpath(dir, "float.tiff"), pixels; sample_format = 3))

        @test_throws "32-bit samples" open_tiff(
            write_tiff(joinpath(dir, "bits.tiff"), pixels; bits_per_sample = 16))

        # A strip table shorter than the raster is tall leaves lines unaddressed.
        @test_throws "one offset per line" open_tiff(
            write_tiff(joinpath(dir, "short.tiff"), pixels; nstrips = 2))

        # A byte count that is not `4 * nsamples` means a line is not the run the tags imply, so
        # indexing would read across line boundaries.
        @test_throws "occupy" open_tiff(
            write_tiff(joinpath(dir, "bytes.tiff"), pixels; strip_bytes = 12))
    end
end

@testset "a file that is not a TIFF is refused" begin
    mktempdir() do dir
        p = joinpath(dir, "not.tiff")
        write(p, "this is not a TIFF file at all, just some bytes")
        @test_throws "byte-order mark" open_tiff(p)

        write(joinpath(dir, "tiny.tiff"), UInt8[0x49, 0x49, 0x2a])
        @test_throws "too short" open_tiff(joinpath(dir, "tiny.tiff"))

        @test_throws "not a file" open_tiff(joinpath(dir, "absent.tiff"))
    end
end

# The one property the whole design rests on: a real measurement raster is uncompressed with one strip
# per line, so its lines are addressable without decoding. Runs only where a granule is available.
if haskey(ENV, "SLCDATASETS_S1_TIFF")
    @testset "a real measurement raster has the layout the reader assumes" begin
        t = open_tiff(ENV["SLCDATASETS_S1_TIFF"])
        @test size(t, 1) > 0
        @test size(t, 2) > 0
        # Window and scalar reads agree deep inside a multi-gigabyte file.
        i, j = size(t, 1) ÷ 2, size(t, 2) ÷ 2
        @test t[i:(i + 3), j:(j + 4)] == [t[a, b] for a in i:(i + 3), b in j:(j + 4)]
    end
end
