# The reader benchmark: opening a product, and reading its orbit, for each format.
#
# Run as `julia --project=bench bench/bench.jl`. The products are synthetic but sized like real ones,
# so the numbers track what a real granule costs in parse work; they do not track disk or network.

using BenchmarkTools
using Printf
using SLCDatasets
using SLCDatasets: orbit

include(joinpath(@__DIR__, "gen_synthetic.jl"))

const WORK = get(ENV, "SLC_BENCH_DIR", mktempdir(; prefix = "slcbench_"))

function products()
    mkpath(WORK)
    zip = joinpath(WORK, "S1A_synthetic.zip")
    dir = joinpath(WORK, "S1A_synthetic.SAFE")
    eof = joinpath(WORK, "S1A_OPER_AUX_POEORB_synthetic.EOF")
    isfile(zip) || write_safe_zip(zip)
    isdir(dir) || write_safe_dir(dir)
    isfile(eof) || write_eof(eof)
    return (; zip, dir, eof)
end

const P = products()

@printf("annotation XML : %.1f MB per subswath\n",
        length(annotation_document(1)) / 1e6)
@printf("zip            : %.1f MB\n", filesize(P.zip) / 1e6)
@printf("EOF            : %.1f MB\n", filesize(P.eof) / 1e6)
println()

const CASES = [
    "S1 mosaic, zip" =>
        () -> open_slc(P.zip; orbit = P.eof),
    "S1 mosaic, dir" =>
        () -> open_slc(P.dir; orbit = P.eof),
    "S1 burst, zip" =>
        () -> open_slc(P.zip; orbit = P.eof, swath = 2, burst = 5),
    "S1 all 9 bursts of IW2, zip" =>
        () -> [open_slc(P.zip; orbit = P.eof, swath = 2, burst = i) for i in 1:N_BURSTS],
    "S1 all 9 bursts via bursts()" =>
        () -> collect(bursts(P.zip; orbit = P.eof, swath = 2)),
    "S1 nbursts" =>
        () -> nbursts(P.zip; swath = 2),
    "S1 orbit (mosaic)" =>
        () -> orbit(open_slc(P.zip; orbit = P.eof)),
]

@printf("%-30s %10s %12s %10s\n", "case", "time", "allocs", "memory")
println("-"^66)
for (name, f) in CASES
    f()  # warm up compilation before measuring
    t = @benchmark $f() samples = 20 evals = 1 seconds = 30
    @printf("%-30s %10s %12d %10s\n", name,
            BenchmarkTools.prettytime(minimum(t.times)),
            minimum(t.allocs),
            BenchmarkTools.prettymemory(minimum(t.memory)))
end
