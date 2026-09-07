# Benchmarks

Synthetic products sized like real granules, so the readers can be measured without a multi-gigabyte
file or a network. `gen_synthetic.jl` writes a three-subswath SAFE product — as a zip and as a
directory, with a 2.3 MB annotation document per subswath — and a POEORB-shaped `.EOF` carrying the
9,361 state vectors a real one does.

```
julia --project=bench bench/bench.jl         # the reader, per format and per view
julia --project=bench bench/bench_nisar.jl   # the NISAR path, piece by piece
julia --project=bench bench/profile.jl       # where one Sentinel-1 read spends its time
```

`SLC_BENCH_DIR` sets where the products are written, so a run reuses them instead of regenerating.

## Checking a change did not alter what is read

`equiv.jl` dumps every field of every reader output — both mosaics, all three subswaths, all 27
bursts, and the state vectors — with the floats as hex literals, so a comparison is bit-exact rather
than to a printed precision. Run it before and after a change and diff:

```
git worktree add /tmp/slc_baseline HEAD
julia --project=bench bench/equiv.jl dump /tmp/before.txt   # in the worktree
julia --project=bench bench/equiv.jl dump /tmp/after.txt    # here
diff /tmp/before.txt /tmp/after.txt
```

This holds the reader against its own previous output over a wider spread of views than the test suite
covers — every burst of every subswath, both mosaics, and subswath subsets. Agreeing with yesterday's
output is not agreeing with ISCE3, so it complements rather than replaces `test/sentinel1.jl`, whose
golden values are ISCE3's and run without any granule.
