# Every field of every reader output, before and after the restructure, compared bitwise.
#
# Run as two passes over the same synthetic products: `julia --project=bench bench/equiv.jl dump <out>`
# under each version, then `diff` the two dumps. The scalars are written as hex float literals so the
# comparison is bit-exact rather than to a printed precision.

using Printf
using SLCDatasets
using SLCDatasets: orbit

include(joinpath(@__DIR__, "gen_synthetic.jl"))

const WORK = get(ENV, "SLC_BENCH_DIR", "/tmp/slcbench")
mkpath(WORK)
const ZIP = joinpath(WORK, "S1A_synthetic.zip")
const DIR = joinpath(WORK, "S1A_synthetic.SAFE")
const EOF_ = joinpath(WORK, "S1A_OPER_AUX_POEORB_synthetic.EOF")
isfile(ZIP) || write_safe_zip(ZIP)
isdir(DIR) || write_safe_dir(DIR)
isfile(EOF_) || write_eof(EOF_)

h(x::Float64) = @sprintf("%a", x)
h(x) = string(x)

function dump_slc(io, label, s)
    id, g = s.identification, s.geometry
    for (k, v) in [("mission", id.mission), ("product_type", id.product_type),
                   ("absolute_orbit", id.absolute_orbit), ("pass", id.pass_direction),
                   ("look_direction", id.look_direction), ("start_time", id.start_time),
                   ("stop_time", id.stop_time), ("polygon", id.bounding_polygon),
                   ("starting_range", g.starting_range), ("far_range", g.far_range),
                   ("range_pixel_spacing", g.range_pixel_spacing), ("wavelength", g.wavelength),
                   ("prf", g.prf), ("sensing_start", g.sensing_start),
                   ("sensing_stop", g.sensing_stop), ("nlines", g.nlines),
                   ("nsamples", g.nsamples), ("look_side", g.look_side), ("epoch", g.epoch)]
        println(io, label, "\t", k, "\t", h(v))
    end
    o = orbit(s)
    println(io, label, "\torbit.n\t", length(o.time))
    println(io, label, "\torbit.epoch\t", o.epoch)
    println(io, label, "\torbit.interp\t", o.interp_method)
    println(io, label, "\torbit.kind\t", o.kind)
    for i in eachindex(o.time)
        println(io, label, "\torbit.t[", i, "]\t", h(o.time[i]))
        println(io, label, "\torbit.p[", i, "]\t", join(h.(o.position[i]), ","))
        println(io, label, "\torbit.v[", i, "]\t", join(h.(o.velocity[i]), ","))
    end
end

out = length(ARGS) >= 2 ? ARGS[2] : "/tmp/equiv.txt"
open(out, "w") do io
    dump_slc(io, "mosaic_zip", open_slc(ZIP; orbit = EOF_))
    dump_slc(io, "mosaic_dir", open_slc(DIR; orbit = EOF_))
    for sw in 1:3
        dump_slc(io, "swath$sw", open_slc(ZIP; orbit = EOF_, swath = sw))
        for i in 1:N_BURSTS
            dump_slc(io, "burst_$(sw)_$i", open_slc(ZIP; orbit = EOF_, swath = sw, burst = i))
        end
    end
    dump_slc(io, "swaths12", open_slc(ZIP; orbit = EOF_, swaths = [1, 2]))
    dump_slc(io, "swaths23", open_slc(ZIP; orbit = EOF_, swaths = [2, 3]))
    for sw in 1:3
        println(io, "nbursts\t$sw\t", nbursts(ZIP; swath = sw))
    end
    println(io, "polarizations\t-\t", join(SLCDatasets.safe_polarizations(ZIP), ","))
    println(io, "default_pol\t-\t", SLCDatasets.default_polarization(ZIP))
    println(io, "is_safe_zip\t-\t", SLCDatasets.is_safe_product(ZIP))
    println(io, "is_safe_dir\t-\t", SLCDatasets.is_safe_product(DIR))
end
println("wrote $out")
