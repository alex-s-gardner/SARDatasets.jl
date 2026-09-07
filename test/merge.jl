# Merging a subswath's bursts into one acquisition, and reading its samples.
#
# The property that matters is a time identity: grid row `n` holds the samples recorded at
# `epoch + sensing_start + (n - 1) / prf`, and for the burst that owns the row the same instant is that
# burst's own recorded start plus its own row offset. Asserting the two agree checks the placement
# against the annotation directly, in exact arithmetic over timestamps.
#
# It is checked that way rather than by correlating the overlap between adjacent bursts, which cannot
# work: consecutive TOPS bursts view the same ground about 2.7 s apart at opposite squint, so they are
# fully decorrelated — measured coherence over the overlap is 0.001 — and multilooked amplitude gives a
# correlation surface that is flat in azimuth to four decimal places, because the azimuth profile over
# the overlap is uniform to half a percent. There is no peak to find. The time identity is sharper
# anyway: it holds to under half a microsecond against a line interval of two milliseconds.

using SLCDatasets
using SLCDatasets: MergedBurstBackend, ConcatenatedBursts, BurstValidMask, Amplitude,
                   read_annotation, annotation_xml, burst_grid, measurement_path, open_tiff,
                   Sentinel1Product, UtcTime, seconds_between, valid_lines, valid_samples,
                   MAX_BURST_GRID_RESIDUAL, grid
using Test

# A merged acquisition needs only annotation, so these run against every committed case.
if !isempty(S1_PRODUCTS)
    @testset "a merged subswath spans its bursts" begin
        for p in S1_PRODUCTS
            pol = String(p.gold.polarization)
            for swath in 1:3
                b = bursts(p.safe; orbit = p.eof, swath, polarization = pol)
                m = merge_bursts(b)
                a = read_annotation(annotation_xml(p.safe, swath, pol), swath)
                g = m.geometry

                @test g.nlines == burst_grid(a).nlines
                @test g.nsamples == a.samples_per_burst
                # Shorter than the stack of whole bursts, because the overlap appears once.
                @test g.nlines < a.lines_per_burst * length(b)
                # ...and longer than any single burst.
                @test g.nlines > a.lines_per_burst

                # Range geometry is the subswath's, unchanged by the merge.
                @test g.starting_range == a.starting_range
                @test g.range_pixel_spacing == a.range_pixel_spacing
                @test g.wavelength == a.wavelength
                @test g.prf == 1 / a.azimuth_time_interval

                # The azimuth window and the row count must describe each other: the window is
                # `nlines - 1` line intervals long. This is what `prf` claims.
                @test (g.sensing_stop - g.sensing_start) * g.prf ≈ g.nlines - 1 rtol = 1e-9

                # A merge begins where its first burst begins.
                @test g.epoch == first(b).geometry.epoch
                @test g.sensing_start == first(b).geometry.sensing_start
                @test m.identification.absolute_orbit == first(b).identification.absolute_orbit
                @test m.identification.mission == first(b).identification.mission
            end
        end
    end

    # Every burst's own recorded time must land on the row the grid gives it. A one-row error shows here
    # as a two-millisecond disagreement, four thousand times the tolerance.
    @testset "a row's time is the owning burst's own time" begin
        for p in S1_PRODUCTS
            pol = String(p.gold.polarization)
            for swath in 1:3
                b = bursts(p.safe; orbit = p.eof, swath, polarization = pol)
                m = merge_bursts(b)
                a = read_annotation(annotation_xml(p.safe, swath, pol), swath)
                dt = a.azimuth_time_interval
                anchor = a.burst_start[first(grid(m).placements).burst]

                for pl in grid(m).placements
                    from_burst = seconds_between(anchor, a.burst_start[pl.burst])
                    for (row, burst_row) in ((first(pl.grid_rows), first(pl.burst_rows)),
                                             (last(pl.grid_rows), last(pl.burst_rows)))
                        # The instant the grid puts this row at, and the instant the burst records for
                        # the row it comes from, both as seconds since the anchor.
                        @test (row - 1) * dt ≈ from_burst + (burst_row - 1) * dt atol = 1e-6
                    end
                end
            end
        end
    end

    # Each burst carries its own epoch, so a subswath holds as many epochs as bursts and their
    # `sensing_start` values are not on one scale. The merge puts them on one; without that the times
    # would be out by the whole spread of the epochs.
    @testset "the bursts' epochs are reconciled" begin
        p = first(S1_PRODUCTS)
        pol = String(p.gold.polarization)
        b = bursts(p.safe; orbit = p.eof, swath = 2, polarization = pol)
        m = merge_bursts(b)
        g = m.geometry
        dt = 1 / g.prf

        # The premise: the bursts really do report against different epochs.
        @test length(unique(s.geometry.epoch for s in b)) == length(b)

        # So differencing their `sensing_start` values is meaningless. The last burst starts after the
        # first, but the raw difference is negative — it is the trap this reconciliation avoids.
        raw = last(b).geometry.sensing_start - first(b).geometry.sensing_start
        elapsed = repeat_interval(first(b), last(b))
        @test elapsed > 0
        @test raw < 0
        @test !isapprox(raw, elapsed; atol = 1.0)

        # Every burst's start, as an instant, must be the instant the merged grid puts its first row at.
        origin = UtcTime(g.epoch, g.sensing_start)
        for pl in grid(m).placements
            burst = b[pl.burst]
            start = UtcTime(burst.geometry.epoch, burst.geometry.sensing_start)
            row = first(pl.grid_rows) - (first(pl.burst_rows) - 1)
            @test seconds_between(origin, start) ≈ (row - 1) * dt atol = 1e-6
        end
    end

    # The orbit must bracket the whole merged window, not one burst's: a solve anywhere in the image
    # interpolates rather than extrapolates. This is what the geometry consumer asserts downstream.
    @testset "the orbit covers the merged window" begin
        p = first(S1_PRODUCTS)
        pol = String(p.gold.polarization)
        m = merge_bursts(bursts(p.safe; orbit = p.eof, swath = 2, polarization = pol))
        sv = orbit(m)
        g = m.geometry
        @test sv.epoch == g.epoch
        @test first(sv.time) <= g.sensing_start
        @test g.sensing_stop <= last(sv.time)
    end

    @testset "a subset of a subswath merges" begin
        p = first(S1_PRODUCTS)
        pol = String(p.gold.polarization)
        b = bursts(p.safe; orbit = p.eof, swath = 2, polarization = pol)

        part = merge_bursts(b[3:7])
        @test length(grid(part).placements) == 5
        @test nlines(part) < nlines(merge_bursts(b))
        # A subset is anchored on its own first burst.
        @test part.geometry.epoch == b[3].geometry.epoch
        @test part.geometry.sensing_start == b[3].geometry.sensing_start

        # A merge of one burst spans that burst's own extent.
        one = merge_bursts(b[4:4])
        @test nlines(one) == nlines(b[4])
        @test one.geometry.sensing_start == b[4].geometry.sensing_start
        @test one.geometry.prf == b[4].geometry.prf
    end

    # A merge is over bursts of one subswath in order. Every other arrangement describes something the
    # result could not honestly report a single geometry for.
    @testset "bursts that do not form one subswath are refused" begin
        p = first(S1_PRODUCTS)
        pol = String(p.gold.polarization)
        b2 = bursts(p.safe; orbit = p.eof, swath = 2, polarization = pol)
        b3 = bursts(p.safe; orbit = p.eof, swath = 3, polarization = pol)

        @test_throws "at least one burst" merge_bursts(SLC[])
        @test_throws "different slant ranges" merge_bursts([b2[1], b3[2]])
        @test_throws "consecutive ascending run" merge_bursts([b2[1], b2[3]])
        @test_throws "consecutive ascending run" merge_bursts([b2[3], b2[2]])
        @test_throws "consecutive ascending run" merge_bursts([b2[2], b2[2]])

        # The mosaic across subswaths is not a burst of one.
        mosaic = open_slc(p.safe; orbit = p.eof, polarization = pol)
        @test_throws "rather than a burst" merge_bursts([mosaic])
    end

    @testset "an acquisition with no raster says so" begin
        p = first(S1_PRODUCTS)
        pol = String(p.gold.polarization)
        m = merge_bursts(bursts(p.safe; orbit = p.eof, swath = 2, polarization = pol))
        if !isdir(p.safe)
            # A zipped product: the raster is deflated inside it, so a window is not addressable.
            @test_throws "zipped Sentinel-1 product" pixels(m)
        else
            # The rebuilt fixture carries annotation only.
            @test_throws "annotation only" pixels(m)
        end

        # A backend that reads metadata alone should say so rather than fail on dispatch with a
        # `MethodError` naming an internal type.
        mktempdir() do dir
            nisar = open_slc(write_fixture_product(joinpath(dir, "fixture_rslc.h5")))
            @test_throws "reads metadata only" pixels(nisar)
            @test_throws "not a merge of bursts" grid(nisar)
        end
        # A single burst is not a merge either, so it has no grid.
        @test_throws "not a merge of bursts" grid(
            open_slc(p.safe; orbit = p.eof, swath = 2, burst = 1, polarization = pol))
    end
end

# Reading samples needs an unpacked `.SAFE`, so these run only where one is available.
const MERGE_SAFE = get(ENV, "SLCDATASETS_S1_SAFE", "")
const MERGE_EOF = get(ENV, "SLCDATASETS_S1_SAFE_ORBIT", "")

if !isempty(MERGE_SAFE) && !isempty(MERGE_EOF)
    @testset "the merged samples come from the right burst" begin
        swath = parse(Int, get(ENV, "SLCDATASETS_S1_SAFE_SWATH", "2"))
        pol = lowercase(get(ENV, "SLCDATASETS_S1_SAFE_POL", "hh"))
        b = bursts(MERGE_SAFE; orbit = MERGE_EOF, swath, polarization = pol)
        m = merge_bursts(b)
        a = read_annotation(annotation_xml(MERGE_SAFE, swath, pol), swath)
        raster = open_tiff(measurement_path(
            Sentinel1Product(MERGE_SAFE; orbit = MERGE_EOF, polarization = pol), swath))

        px = pixels(m)
        mask = validmask(m)
        @test size(px) == (m.geometry.nlines, m.geometry.nsamples)
        @test size(mask) == size(px)
        @test eltype(px) == Complex{Int16}

        # A `.SAFE` raster holds every burst stacked, so the sample a placed row reads must be the one
        # at that burst's own offset in the file.
        L = a.lines_per_burst
        for pl in grid(m).placements
            for (row, burst_row) in ((first(pl.grid_rows), first(pl.burst_rows)),
                                     (last(pl.grid_rows), last(pl.burst_rows)))
                file_row = (pl.burst - 1) * L + burst_row
                for col in (first(pl.cols), first(pl.cols) + 137, last(pl.cols))
                    @test px[row, col] == raster[file_row, col]
                    @test mask[row, col]
                end
            end
        end
    end

    @testset "samples no burst imaged read as zero" begin
        swath = parse(Int, get(ENV, "SLCDATASETS_S1_SAFE_SWATH", "2"))
        pol = lowercase(get(ENV, "SLCDATASETS_S1_SAFE_POL", "hh"))
        m = merge_bursts(bursts(MERGE_SAFE; orbit = MERGE_EOF, swath, polarization = pol))
        px = pixels(m)
        mask = validmask(m)
        placements = grid(m).placements

        # The margin before the first placed row belongs to no burst. The raster has data there — the
        # invalid margin of a burst is not zeroed in the file — so this asserts the clip, not the file.
        top = first(placements[1].grid_rows)
        if top > 1
            @test all(iszero, px[1:(top - 1), 1000:1010])
            @test !any(mask[1:(top - 1), 1000:1010])
        end

        # Samples outside the valid range of the burst owning a row.
        pl = placements[1]
        row = first(pl.grid_rows) + 10
        if first(pl.cols) > 1
            @test iszero(px[row, first(pl.cols) - 1])
            @test !mask[row, first(pl.cols) - 1]
        end
        if last(pl.cols) < size(px, 2)
            @test iszero(px[row, last(pl.cols) + 1])
            @test !mask[row, last(pl.cols) + 1]
        end
        @test mask[row, first(pl.cols)]
        @test mask[row, last(pl.cols)]
    end

    # Window reads are the path a correlator takes, so they must agree with scalar reads everywhere —
    # including across a seam, where a window spans two bursts.
    @testset "a window agrees with scalar reads" begin
        swath = parse(Int, get(ENV, "SLCDATASETS_S1_SAFE_SWATH", "2"))
        pol = lowercase(get(ENV, "SLCDATASETS_S1_SAFE_POL", "hh"))
        m = merge_bursts(bursts(MERGE_SAFE; orbit = MERGE_EOF, swath, polarization = pol))
        px = pixels(m)
        mask = validmask(m)
        amp = amplitude(m)
        placements = grid(m).placements

        seam = first(placements[2].grid_rows)
        for (rows, cols) in (((seam - 3):(seam + 3), 12000:12009),      # across a seam
                             (1:6, 1:8),                                # the top-left margin
                             (100:107, 1000:1007),                      # inside one burst
                             ((size(px, 1) - 5):size(px, 1), (size(px, 2) - 7):size(px, 2)))
            win = px[rows, cols]
            mwin = mask[rows, cols]
            awin = amp[rows, cols]
            @test size(win) == (length(rows), length(cols))
            for (di, i) in enumerate(rows), (dj, j) in enumerate(cols)
                @test win[di, dj] == px[i, j]
                @test mwin[di, dj] == mask[i, j]
                @test awin[di, dj] == Float32(abs(px[i, j]))
            end
        end

        @test eltype(amp) == Float32
        @test size(amp) == size(px)
        # `amplitude` wraps the samples rather than copying them, so its parent is the same kind of
        # array `pixels` returns and reads the same values.
        @test parent(amp) isa typeof(px)
        @test parent(amplitude(px)) === px
        @test validmask(amp) == mask
    end
end
