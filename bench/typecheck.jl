# Type stability and dispatch of the reader entry points.
#
# `@report_opt` finds runtime dispatch and `@report_call` finds inference errors; `Base.infer_return_type`
# checks each entry point resolves to a concrete type. Run as
# `julia --project=bench bench/typecheck.jl`.

using JET
using InteractiveUtils
using SLCDatasets
using SLCDatasets: orbit, Sentinel1Backend, NisarBackend, read_identification, read_geometry,
                   read_orbit, read_annotation, read_eof_state_vectors, annotation_xml,
                   parse_utc, seconds_between, UtcTime, safe_polarizations, is_safe_product,
                   annotation, _utc_string, parse_cf_epoch, _truncated_datetime

include(joinpath(@__DIR__, "gen_synthetic.jl"))

const WORK = get(ENV, "SLC_BENCH_DIR", "/tmp/slcbench")
mkpath(WORK)
const ZIP = joinpath(WORK, "S1A_synthetic.zip")
const DIR = joinpath(WORK, "S1A_synthetic.SAFE")
const EOF_ = joinpath(WORK, "S1A_OPER_AUX_POEORB_synthetic.EOF")
const H5 = joinpath(WORK, "nisar_fixture.h5")
isfile(ZIP) || write_safe_zip(ZIP)
isdir(DIR) || write_safe_dir(DIR)
isfile(EOF_) || write_eof(EOF_)

isfile(H5) || error("run bench/bench_nisar.jl first to write $H5")

const XML1 = annotation_xml(DIR, 1, "vv")
const PROD = Sentinel1Product(DIR; orbit = EOF_)
const B_MOSAIC = Sentinel1Backend(PROD)
const B_BURST = Sentinel1Backend(PROD, 2, 5)
const B_NISAR = NisarBackend(H5, "LSAR", "RSLC", "A")

println("="^78)
println("inferred return types")
println("="^78)
for (label, f, types) in [
        ("read_annotation", read_annotation, (String, Int)),
        ("read_identification(S1)", read_identification, (Sentinel1Backend,)),
        ("read_geometry(S1)", read_geometry, (Sentinel1Backend,)),
        ("read_orbit(S1)", read_orbit, (Sentinel1Backend,)),
        ("read_identification(NISAR)", read_identification, (NisarBackend,)),
        ("read_geometry(NISAR)", read_geometry, (NisarBackend,)),
        ("read_orbit(NISAR)", read_orbit, (NisarBackend,)),
        ("read_eof_state_vectors", read_eof_state_vectors, (String,)),
        ("parse_utc", parse_utc, (String,)),
        ("parse_utc(::SubString)", parse_utc, (SubString{String},)),
        ("seconds_between", seconds_between, (UtcTime, UtcTime)),
        ("_utc_string", _utc_string, (UtcTime,)),
        ("parse_cf_epoch", parse_cf_epoch, (String,)),
        ("_truncated_datetime", _truncated_datetime, (String,)),
        ("annotation", annotation, (Sentinel1Product, Int)),
        ("nbursts(product)", nbursts, (Sentinel1Product, Int)),
        ("orbit(::SLC)", orbit, (SLC{Sentinel1Backend},)),
        ("epoch_offset", epoch_offset, (RadarGeometry,)),
        ("safe_polarizations", safe_polarizations, (String,)),
        ("is_safe_product", is_safe_product, (String,)),
    ]
    t = Base.infer_return_type(f, types)
    concrete = isconcretetype(t) || t === Union{} ? "" : "   <-- not concrete"
    println(rpad(label, 30), t, concrete)
end

# The lazy burst vector: indexing and iteration must infer, or a consumer's loop over it boxes.
series = bursts(DIR; orbit = EOF_, swath = 2)
println()
println(rpad("eltype(bursts(...))", 30), eltype(series))
println(rpad("getindex", 30), Base.infer_return_type(getindex, (typeof(series), Int)))
println(rpad("iterate", 30), Base.infer_return_type(iterate, (typeof(series),)))
println(rpad("collect", 30), Base.infer_return_type(collect, (typeof(series),)))
println(rpad("first", 30), Base.infer_return_type(first, (typeof(series),)))

println()
println("="^78)
println("@report_opt: runtime dispatch")
println("="^78)
# Scoped to this package: EzXML and HDF5 dispatch dynamically inside their own wrappers, and Base's
# `show` error paths are reachable from any `string(...)` in a message. Neither is this package's to fix,
# and either buries what is.
for (label, thunk) in [
        ("read_annotation", () -> @report_opt target_modules = (SLCDatasets,) read_annotation(XML1, 1)),
        ("read_geometry(S1 mosaic)", () -> @report_opt target_modules = (SLCDatasets,) read_geometry(B_MOSAIC)),
        ("read_geometry(S1 burst)", () -> @report_opt target_modules = (SLCDatasets,) read_geometry(B_BURST)),
        ("read_identification(S1)", () -> @report_opt target_modules = (SLCDatasets,) read_identification(B_MOSAIC)),
        ("read_orbit(S1)", () -> @report_opt target_modules = (SLCDatasets,) read_orbit(B_MOSAIC)),
        ("read_eof_state_vectors", () -> @report_opt target_modules = (SLCDatasets,) read_eof_state_vectors(EOF_)),
        ("parse_utc", () -> @report_opt target_modules = (SLCDatasets,) parse_utc("2019-01-01T12:00:00.123456")),
        ("read_identification(NISAR)", () -> @report_opt target_modules = (SLCDatasets,) read_identification(B_NISAR)),
        ("read_geometry(NISAR)", () -> @report_opt target_modules = (SLCDatasets,) read_geometry(B_NISAR)),
        ("read_orbit(NISAR)", () -> @report_opt target_modules = (SLCDatasets,) read_orbit(B_NISAR)),
        ("getindex(bursts)", () -> @report_opt target_modules = (SLCDatasets,) series[3]),
        ("orbit(::SLC)", () -> @report_opt target_modules = (SLCDatasets,) orbit(series[3])),
    ]
    println("\n--- ", label, " ---")
    show(stdout, thunk())
    println()
end

println()
println("="^78)
println("@report_call: inference errors")
println("="^78)
for (label, thunk) in [
        ("open_slc(S1 dir)", () -> @report_call target_modules = (SLCDatasets,) open_slc(DIR; orbit = EOF_)),
        ("open_slc(NISAR)", () -> @report_call target_modules = (SLCDatasets,) open_slc(H5)),
        ("bursts", () -> @report_call target_modules = (SLCDatasets,) bursts(DIR; orbit = EOF_, swath = 2)),
        ("read_annotation", () -> @report_call target_modules = (SLCDatasets,) read_annotation(XML1, 1)),
        ("read_orbit(S1)", () -> @report_call target_modules = (SLCDatasets,) read_orbit(B_MOSAIC)),
        ("read_orbit(NISAR)", () -> @report_call target_modules = (SLCDatasets,) read_orbit(B_NISAR)),
        ("parse_utc", () -> @report_call target_modules = (SLCDatasets,) parse_utc("2019-01-01T12:00:00.123456")),
    ]
    println("\n--- ", label, " ---")
    show(stdout, thunk())
    println()
end
