# Where the time in a read goes: the pieces of one open, measured separately.

using BenchmarkTools
using Printf
using SLCDatasets
using SLCDatasets: Sentinel1Backend, read_annotation, read_identification, read_geometry,
                   read_orbit, read_eof_state_vectors, annotation_xml, safe_polarizations,
                   default_polarization, is_safe_product, parse_utc, with_zip, zip_names
using EzXML: readxml, parsexml, root

include(joinpath(@__DIR__, "gen_synthetic.jl"))

const WORK = get(ENV, "SLC_BENCH_DIR", "/tmp/slcbench")
mkpath(WORK)
const ZIP = joinpath(WORK, "S1A_synthetic.zip")
const DIR = joinpath(WORK, "S1A_synthetic.SAFE")
const EOF_ = joinpath(WORK, "S1A_OPER_AUX_POEORB_synthetic.EOF")
isfile(ZIP) || write_safe_zip(ZIP)
isdir(DIR) || write_safe_dir(DIR)
isfile(EOF_) || write_eof(EOF_)

const XML1 = annotation_xml(DIR, 1, "vv")
const EOFTEXT = read(EOF_, String)

b_mosaic = Sentinel1Backend(Sentinel1Product(DIR; orbit = EOF_))
b_burst = Sentinel1Backend(Sentinel1Product(DIR; orbit = EOF_, swaths = [2]), 2, 5)

function show_case(name, f)
    f()
    t = @benchmark $f() samples = 15 evals = 1 seconds = 25
    @printf("%-42s %10s %10d %10s\n", name,
            BenchmarkTools.prettytime(minimum(t.times)), minimum(t.allocs),
            BenchmarkTools.prettymemory(minimum(t.memory)))
end

@printf("%-42s %10s %10s %10s\n", "piece", "time", "allocs", "memory")
println("-"^76)

show_case("read the annotation text (dir)", () -> annotation_xml(DIR, 1, "vv"))
show_case("read the annotation text (zip)", () -> annotation_xml(ZIP, 1, "vv"))
show_case("libxml2 parse only, annotation", () -> parsexml(XML1))
show_case("read_annotation (text -> struct)", () -> read_annotation(XML1, 1))
println()
show_case("is_safe_product (zip)", () -> is_safe_product(ZIP))
show_case("safe_polarizations (zip)", () -> safe_polarizations(ZIP))
show_case("zip central directory listing", () -> with_zip(zip_names, ZIP))
println()
show_case("read_identification (mosaic, dir)", () -> read_identification(b_mosaic))
show_case("read_geometry (mosaic, dir)", () -> read_geometry(b_mosaic))
show_case("read_identification (burst, dir)", () -> read_identification(b_burst))
show_case("read_geometry (burst, dir)", () -> read_geometry(b_burst))
println()
show_case("libxml2 parse only, EOF", () -> parsexml(EOFTEXT))
show_case("read_eof_state_vectors (all 9361)", () -> read_eof_state_vectors(EOF_))
show_case("read_orbit (keeps ~20)", () -> read_orbit(b_mosaic))
show_case("parse_utc x1000", () -> for _ in 1:1000; parse_utc("2019-01-01T12:00:00.123456"); end)
