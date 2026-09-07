# Synthetic products at the scale of real ones, for benchmarking the readers without a multi-gigabyte
# granule. The annotation XML carries the bulk a real IW annotation does — geolocation grid, antenna
# pattern, Doppler and FM-rate lists — because the readers parse the whole document to reach a few dozen
# scalars, so the filler is what the parse actually costs.

using Dates: DateTime
using Printf
using SLCDatasets: UtcTime, _advance, _utc_string
using ZipArchives: ZipWriter, zip_newfile

# The instant the generated products' times are measured from.
const EPOCH = UtcTime(DateTime(2019, 1, 1), 0.0)

const N_BURSTS = 9
const LINES_PER_BURST = 1500
const SAMPLES_PER_BURST = 22000

function burst_xml(i::Int, t0::Float64, az_interval::Float64)
    az = t0 + (i - 1) * LINES_PER_BURST * az_interval
    stamp = utc_stamp(az)
    fvs = join((j < 20 || j > LINES_PER_BURST - 20 ? -1 : 300 for j in 1:LINES_PER_BURST), " ")
    lvs = join((j < 20 || j > LINES_PER_BURST - 20 ? -1 : 21000 for j in 1:LINES_PER_BURST), " ")
    return """
      <burst>
        <azimuthTime>$stamp</azimuthTime>
        <azimuthAnxTime>$(@sprintf("%.6f", 2000.0 + i))</azimuthAnxTime>
        <sensingTime>$stamp</sensingTime>
        <byteOffset>$(i * 1000000)</byteOffset>
        <firstValidSample count="$LINES_PER_BURST">$fvs</firstValidSample>
        <lastValidSample count="$LINES_PER_BURST">$lvs</lastValidSample>
      </burst>"""
end

# 2019-01-01T00:00:00 plus `s` seconds, in the form the annotation writes. The package's own formatter
# is used so the generated products carry exactly the timestamps its parser round-trips.
utc_stamp(s::Real) = _utc_string(_advance(EPOCH, s))

# The geolocation grid a real annotation carries: ~10 azimuth lines x 21 range samples per burst.
function geolocation_grid(t0::Float64, az_interval::Float64)
    io = IOBuffer()
    n = N_BURSTS * 10 * 21
    print(io, "<geolocationGrid><geolocationGridPointList count=\"$n\">")
    for b in 1:N_BURSTS, l in 1:10, p in 1:21
        line = (b - 1) * LINES_PER_BURST + (l - 1) * 150
        px = (p - 1) * 1100
        print(io, "<geolocationGridPoint><azimuthTime>",
              utc_stamp(t0 + line * az_interval), "</azimuthTime>",
              "<line>", line, "</line><pixel>", px, "</pixel>",
              "<latitude>", @sprintf("%.10e", 6.5e1 + 1e-4 * line), "</latitude>",
              "<longitude>", @sprintf("%.10e", -4.5e1 + 1e-4 * px), "</longitude>",
              "<height>", @sprintf("%.6e", 1.0e2), "</height>",
              "<incidenceAngle>", @sprintf("%.10e", 3.0e1 + 1e-4 * px), "</incidenceAngle>",
              "<elevationAngle>", @sprintf("%.10e", 2.7e1), "</elevationAngle>",
              "<slantRangeTime>", @sprintf("%.10e", 5.3e-3), "</slantRangeTime>",
              "</geolocationGridPoint>")
    end
    print(io, "</geolocationGridPointList></geolocationGrid>")
    return String(take!(io))
end

# The antenna pattern list: one record per burst, each a few thousand space-separated values.
function antenna_pattern(t0::Float64)
    io = IOBuffer()
    print(io, "<antennaPattern><antennaPatternList count=\"$N_BURSTS\">")
    vals = join((@sprintf("%.6e", 1.0 + 1e-6 * k) for k in 1:2000), " ")
    for b in 1:N_BURSTS
        print(io, "<antennaPattern><azimuthTime>", utc_stamp(t0 + 3.0 * b), "</azimuthTime>",
              "<slantRangeTime count=\"2000\">", vals, "</slantRangeTime>",
              "<elevationAngle count=\"2000\">", vals, "</elevationAngle>",
              "<elevationPattern count=\"4000\">", vals, " ", vals, "</elevationPattern>",
              "<incidenceAngle count=\"2000\">", vals, "</incidenceAngle>",
              "<terrainHeight count=\"2000\">", vals, "</terrainHeight>",
              "<roll>0.0</roll></antennaPattern>")
    end
    print(io, "</antennaPatternList></antennaPattern>")
    return String(take!(io))
end

function poly_list(tag, inner_tag, t0::Float64, n::Int)
    io = IOBuffer()
    print(io, "<$tag count=\"$n\">")
    for k in 1:n
        print(io, "<$inner_tag><azimuthTime>", utc_stamp(t0 + 3.0 * k), "</azimuthTime>",
              "<t0>5.3e-3</t0>",
              "<", tag == "azimuthFmRateList" ? "azimuthFmRatePolynomial" : "dataDcPolynomial",
              " count=\"3\">-2.7e3 4.1e5 -1.5e8</",
              tag == "azimuthFmRateList" ? "azimuthFmRatePolynomial" : "dataDcPolynomial", ">",
              "</$inner_tag>")
    end
    print(io, "</$tag>")
    return String(take!(io))
end

"""
    annotation_document(swath) -> String

One subswath's annotation XML, at the size and shape of a real IW SLC annotation.
"""
function annotation_document(swath::Int; polarization = "vv", mission = "S1A")
    az_interval = 2.0555563e-3
    # IW subswaths sit at increasing slant range; the times stagger so IW1 does not start first.
    t0 = 43200.0 + 0.8 * (swath - 1)
    slant_range_time = 5.3e-3 + 1.2e-3 * (swath - 1)
    return string(
        """<?xml version="1.0" encoding="UTF-8"?>
        <product>
          <adsHeader>
            <missionId>$mission</missionId>
            <productType>SLC</productType>
            <polarisation>$(uppercase(polarization))</polarisation>
            <mode>IW</mode>
            <swath>IW$swath</swath>
            <startTime>$(utc_stamp(t0))</startTime>
            <stopTime>$(utc_stamp(t0 + 28.0))</stopTime>
            <absoluteOrbitNumber>25514</absoluteOrbitNumber>
            <missionDataTakeId>123456</missionDataTakeId>
            <imageNumber>00$swath</imageNumber>
          </adsHeader>
          <generalAnnotation>
            <productInformation>
              <pass>DESCENDING</pass>
              <timelinessCategory>Fast-24h</timelinessCategory>
              <platformHeading>-1.6e2</platformHeading>
              <projection>Slant Range</projection>
              <rangeSamplingRate>6.4345238e7</rangeSamplingRate>
              <radarFrequency>5.40500045e9</radarFrequency>
              <azimuthSteeringRate>1.590368784</azimuthSteeringRate>
            </productInformation>
          </generalAnnotation>
          <imageAnnotation>
            <imageInformation>
              <productFirstLineUtcTime>$(utc_stamp(t0))</productFirstLineUtcTime>
              <productLastLineUtcTime>$(utc_stamp(t0 + 28.0))</productLastLineUtcTime>
              <ascendingNodeTime>$(utc_stamp(t0 - 2000.0))</ascendingNodeTime>
              <azimuthTimeInterval>$(@sprintf("%.10e", az_interval))</azimuthTimeInterval>
              <slantRangeTime>$(@sprintf("%.10e", slant_range_time))</slantRangeTime>
              <rangePixelSpacing>2.329562</rangePixelSpacing>
              <numberOfLines>$(N_BURSTS * LINES_PER_BURST)</numberOfLines>
              <numberOfSamples>$SAMPLES_PER_BURST</numberOfSamples>
            </imageInformation>
          </imageAnnotation>
          <swathTiming>
            <linesPerBurst>$LINES_PER_BURST</linesPerBurst>
            <samplesPerBurst>$SAMPLES_PER_BURST</samplesPerBurst>
            <burstList count="$N_BURSTS">""",
        join((burst_xml(i, t0, az_interval) for i in 1:N_BURSTS), "\n"),
        "</burstList></swathTiming>",
        geolocation_grid(t0, az_interval),
        antenna_pattern(t0),
        "<generalAnnotation>",
        poly_list("azimuthFmRateList", "azimuthFmRate", t0, 40),
        "</generalAnnotation><dopplerCentroid>",
        poly_list("dcEstimateList", "dcEstimate", t0, 40),
        "</dopplerCentroid></product>",
    )
end

"""
    write_safe_zip(path) -> String

A zipped `.SAFE` product carrying three subswaths of annotation and a stand-in measurement entry.
"""
function write_safe_zip(path::AbstractString; polarization = "vv", mission = "S1A")
    ZipWriter(path) do w
        zip_newfile(w, "PRODUCT.SAFE/manifest.safe")
        write(w, "<xfdu:XFDU/>")
        for swath in 1:3
            name = "PRODUCT.SAFE/annotation/$(lowercase(mission))-iw$swath-slc-" *
                   "$(lowercase(polarization))-20190101t120000-20190101t120028-025514-02d2e5-00$swath.xml"
            zip_newfile(w, name)
            write(w, annotation_document(swath; polarization, mission))
            # The calibration and RFI trees hold same-named files, which the reader must not pick up.
            for sub in ("calibration", "rfi")
                zip_newfile(w, replace(name, "annotation/" => "annotation/$sub/"))
                write(w, "<product/>")
            end
        end
        zip_newfile(w, "PRODUCT.SAFE/measurement/dummy.tiff")
        write(w, zeros(UInt8, 4096))
    end
    return path
end

"""
    write_safe_dir(path) -> String

The same product as an unzipped `.SAFE` directory.
"""
function write_safe_dir(path::AbstractString; polarization = "vv", mission = "S1A")
    mkpath(joinpath(path, "annotation"))
    write(joinpath(path, "manifest.safe"), "<xfdu:XFDU/>")
    for swath in 1:3
        name = "$(lowercase(mission))-iw$swath-slc-$(lowercase(polarization))-" *
               "20190101t120000-20190101t120028-025514-02d2e5-00$swath.xml"
        write(joinpath(path, "annotation", name),
              annotation_document(swath; polarization, mission))
    end
    return path
end

"""
    write_eof(path; n = 9361) -> String

A POEORB-shaped orbit file. `n` defaults to the count a real one carries: one state vector every ten
seconds over just past a day.
"""
function write_eof(path::AbstractString; n::Int = 9361)
    io = IOBuffer()
    print(io, """<?xml version="1.0"?>
    <Earth_Explorer_File><Data_Block type="xml"><List_of_OSVs count="$n">""")
    for k in 1:n
        t = utc_stamp(10.0 * (k - 1))
        print(io, "<OSV><TAI>TAI=", t, "</TAI><UTC>UTC=", t, "</UTC><UT1>UT1=", t, "</UT1>",
              "<Absolute_Orbit>+25514</Absolute_Orbit>",
              "<X unit=\"m\">", @sprintf("%.6f", 3.0e6 + 1e3 * k), "</X>",
              "<Y unit=\"m\">", @sprintf("%.6f", -4.0e6 + 2e3 * k), "</Y>",
              "<Z unit=\"m\">", @sprintf("%.6f", 5.0e6 + 3e3 * k), "</Z>",
              "<VX unit=\"m/s\">", @sprintf("%.6f", 1.0e3), "</VX>",
              "<VY unit=\"m/s\">", @sprintf("%.6f", 2.0e3), "</VY>",
              "<VZ unit=\"m/s\">", @sprintf("%.6f", -6.0e3), "</VZ>",
              "<Quality>NOMINAL</Quality></OSV>")
    end
    print(io, "</List_of_OSVs></Data_Block></Earth_Explorer_File>")
    write(path, take!(io))
    return path
end
