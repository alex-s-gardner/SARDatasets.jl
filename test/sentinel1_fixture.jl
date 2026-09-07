# Rebuilding a Sentinel-1 SAFE product and its orbit file from the committed inputs.
#
# `reference/sentinel1_inputs*.json` holds the annotation fields and state vectors the reader reads,
# taken verbatim from the granules the golden values were dumped from — the text as the annotation
# writes it, not a re-rounded float. Writing them back into the group and element layout a real product
# uses gives the reader a product to open, so the whole Sentinel-1 path is exercised against ISCE3's
# values without the 4.8 GB of granules those values came from.
#
# `reference/extract_sentinel1_inputs.jl` is what produced these files.

using JSON3

# The golden-value files, each naming the granule its values were dumped from.
const S1_GOLDEN = ["sentinel1_metadata.json", "sentinel1_metadata_s1b.json"]

# The annotation elements `read_annotation` looks for, in the paths it looks for them. Fields a real
# annotation carries but this reader ignores are omitted rather than filled with stand-ins, so what the
# reader depends on is visible here.
function s1_annotation_document(sw::AbstractDict)
    bursts = join(("<burst><azimuthTime>$t</azimuthTime></burst>"
                   for t in sw["burstAzimuthTimes"]))
    return """<?xml version="1.0" encoding="UTF-8"?>
    <product>
      <adsHeader>
        <missionId>$(sw["missionId"])</missionId>
        <productType>$(sw["productType"])</productType>
        <polarisation>$(sw["polarisation"])</polarisation>
        <absoluteOrbitNumber>$(sw["absoluteOrbitNumber"])</absoluteOrbitNumber>
      </adsHeader>
      <generalAnnotation>
        <productInformation>
          <pass>$(sw["pass"])</pass>
          <rangeSamplingRate>$(sw["rangeSamplingRate"])</rangeSamplingRate>
          <radarFrequency>$(sw["radarFrequency"])</radarFrequency>
        </productInformation>
      </generalAnnotation>
      <imageAnnotation>
        <imageInformation>
          <azimuthTimeInterval>$(sw["azimuthTimeInterval"])</azimuthTimeInterval>
          <slantRangeTime>$(sw["slantRangeTime"])</slantRangeTime>
        </imageInformation>
      </imageAnnotation>
      <swathTiming>
        <linesPerBurst>$(sw["linesPerBurst"])</linesPerBurst>
        <samplesPerBurst>$(sw["samplesPerBurst"])</samplesPerBurst>
        <burstList count="$(length(sw["burstAzimuthTimes"]))">$bursts</burstList>
      </swathTiming>
    </product>
    """
end

function s1_eof_document(vectors)
    osvs = join(("<OSV>" * join(("<$f>$(v[f])</$f>" for f in ("UTC", "X", "Y", "Z", "VX", "VY", "VZ"))) *
                 "</OSV>" for v in vectors))
    return """<?xml version="1.0"?>
    <Earth_Explorer_File><Data_Block type="xml"><List_of_OSVs count="$(length(vectors))">$osvs</List_of_OSVs></Data_Block></Earth_Explorer_File>
    """
end

"""
    write_s1_fixture(dir, inputs) -> (safe, eof)

Write `inputs` as a `.SAFE` directory and a POEORB `.EOF` under `dir`.

The orbit file keeps the name from the granule it was extracted from, because the reader reports
POEORB or RESORB by matching that name.
"""
function write_s1_fixture(dir::AbstractString, inputs)
    pol = String(inputs["polarization"])
    mission = lowercase(String(inputs["swaths"]["1"]["missionId"]))
    safe = joinpath(dir, "fixture.SAFE")
    mkpath(joinpath(safe, "annotation", "calibration"))
    write(joinpath(safe, "manifest.safe"), "<xfdu:XFDU/>")
    for swath in 1:3
        sw = inputs["swaths"][string(swath)]
        name = "$mission-iw$swath-slc-$pol-20000101t000000-20000101t000030-000001-000001-00$swath.xml"
        write(joinpath(safe, "annotation", name), s1_annotation_document(sw))
        # The parallel `calibration` tree holds same-named files the reader must not pick up.
        write(joinpath(safe, "annotation", "calibration", name), "<product/>")
    end
    eof = joinpath(dir, String(inputs["orbit_file"]))
    write(eof, s1_eof_document(inputs["state_vectors"]))
    return safe, eof
end

"""
    s1_fixtures() -> Vector{NamedTuple}

The committed Sentinel-1 cases: each golden-value file paired with the inputs to rebuild its product.
"""
function s1_fixtures()
    out = NamedTuple[]
    for name in S1_GOLDEN
        gold_file = joinpath(@__DIR__, "reference", name)
        input_file = joinpath(@__DIR__, "reference", replace(name, "_metadata" => "_inputs"))
        (isfile(gold_file) && isfile(input_file)) || continue
        push!(out, (; name,
                    gold = JSON3.read(read(gold_file, String)),
                    inputs = JSON3.read(read(input_file, String), Dict{String,Any})))
    end
    return out
end
