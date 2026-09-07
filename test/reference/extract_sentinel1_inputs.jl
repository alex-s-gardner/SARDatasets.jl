# Reduce the real granules the golden values were dumped from to the inputs the reader actually reads,
# so `test/sentinel1.jl` can rebuild an equivalent product instead of needing 4.8 GB of granules.
#
# Run against a directory holding the granules and orbit files named in `sentinel1_metadata*.json`:
#
#     SLCDATASETS_S1_DIR=<dir> julia --project=test test/reference/extract_sentinel1_inputs.jl
#
# Every scalar is stored as the annotation's own text rather than as a parsed number: the reader parses
# that text, so keeping it verbatim makes the rebuilt product yield bit-identical values without this
# script having to round-trip a float correctly.

using JSON3
# Through the package rather than as a test dependency: this is a developer script run by hand against
# granules that are not in the repository, not part of the suite.
using SLCDatasets: annotation_xml, parse_utc, seconds_between, UtcTime, S1_ORBIT_PADDING,
                   _advance, _findtext, _strip_utc_prefix,
                   parsexml, readxml, root, findfirst, findall, eachelement

const GOLDEN = ["sentinel1_metadata.json", "sentinel1_metadata_s1b.json"]

# The valid region of one burst, as the four numbers the per-line arrays reduce to.
#
# Storing the arrays verbatim would be some forty thousand integers per granule to record a rectangle.
# They are constant across a burst's valid lines in every granule measured — `read_valid_region`
# checks that rather than assuming it — so the rectangle is what is kept, and the fixture writes arrays
# back out from it.
function burst_valid_region(node)
    first_sample = [parse(Int, s) for s in eachsplit(_findtext(node, "firstValidSample"))]
    last_sample = [parse(Int, s) for s in eachsplit(_findtext(node, "lastValidSample"))]
    lines = findall(>=(0), first_sample)
    lo, hi = first(lines), last(lines)
    return Dict{String,Any}(
        # 1-based line indices into the burst, as the arrays are indexed here.
        "firstValidLine" => lo,
        "lastValidLine" => hi,
        # The annotation's own 0-based sample bounds.
        "firstValidSample" => first_sample[lo],
        "lastValidSample" => last_sample[lo],
    )
end

# The annotation fields `read_annotation` reads, and nothing else.
function subswath_inputs(xml::AbstractString)
    r = root(parsexml(xml))
    burst_nodes = findall("swathTiming/burstList/burst", r)
    bursts = [_findtext(b, "azimuthTime") for b in burst_nodes]
    return Dict{String,Any}(
        "burstValidRegions" => [burst_valid_region(b) for b in burst_nodes],
        "polarisation" => _findtext(r, "adsHeader/polarisation"),
        "absoluteOrbitNumber" => _findtext(r, "adsHeader/absoluteOrbitNumber"),
        "missionId" => _findtext(r, "adsHeader/missionId"),
        "productType" => _findtext(r, "adsHeader/productType"),
        "pass" => _findtext(r, "generalAnnotation/productInformation/pass"),
        "rangeSamplingRate" => _findtext(r, "generalAnnotation/productInformation/rangeSamplingRate"),
        "radarFrequency" => _findtext(r, "generalAnnotation/productInformation/radarFrequency"),
        "slantRangeTime" => _findtext(r, "imageAnnotation/imageInformation/slantRangeTime"),
        "azimuthTimeInterval" => _findtext(r, "imageAnnotation/imageInformation/azimuthTimeInterval"),
        "linesPerBurst" => _findtext(r, "swathTiming/linesPerBurst"),
        "samplesPerBurst" => _findtext(r, "swathTiming/samplesPerBurst"),
        "burstAzimuthTimes" => bursts,
    )
end

# The state vectors any window in the product could ask for: the acquisition's own span widened by the
# padding, so every per-burst orbit the tests build finds the same vectors the real file holds.
function orbit_inputs(eof::AbstractString, first_burst::UtcTime, last_stop::UtcTime)
    list = findfirst("//Data_Block/List_of_OSVs", root(readxml(eof)))
    out = Vector{Dict{String,String}}()
    for osv in eachelement(list)
        stamp = _findtext(osv, "UTC")
        t = parse_utc(_strip_utc_prefix(stamp))
        seconds_between(first_burst, t) < -2 * S1_ORBIT_PADDING && continue
        seconds_between(last_stop, t) > 2 * S1_ORBIT_PADDING && break
        push!(out, Dict(f => String(_findtext(osv, f))
                        for f in ("UTC", "X", "Y", "Z", "VX", "VY", "VZ")))
    end
    return out
end

dir = get(ENV, "SLCDATASETS_S1_DIR", "")
isempty(dir) && error("set SLCDATASETS_S1_DIR to the directory holding the granules")

for name in GOLDEN
    gold = JSON3.read(read(joinpath(@__DIR__, name), String))
    safe = joinpath(dir, String(gold.safe))
    eof = joinpath(dir, String(gold.orbit_file))
    (isfile(safe) && isfile(eof)) || error("`$safe` or `$eof` is missing")
    pol = String(gold.polarization)

    swaths = Dict{String,Any}()
    starts = UtcTime[]
    stops = UtcTime[]
    for swath in 1:3
        s = subswath_inputs(annotation_xml(safe, swath, pol))
        swaths[string(swath)] = s
        times = parse_utc.(s["burstAzimuthTimes"])
        push!(starts, minimum(times))
        push!(stops, _advance(maximum(times),
                             (parse(Int, s["linesPerBurst"]) - 1) *
                             parse(Float64, s["azimuthTimeInterval"])))
    end

    out = Dict{String,Any}(
        "safe" => String(gold.safe),
        "orbit_file" => String(gold.orbit_file),
        "polarization" => pol,
        "swaths" => swaths,
        "state_vectors" => orbit_inputs(eof, minimum(starts), maximum(stops)),
    )
    dest = joinpath(@__DIR__, replace(name, "_metadata" => "_inputs"))
    open(dest, "w") do io
        JSON3.pretty(io, out)
    end
    println("wrote ", dest, " (", round(filesize(dest) / 1024; digits = 1), " KiB, ",
            length(out["state_vectors"]), " state vectors)")
end
