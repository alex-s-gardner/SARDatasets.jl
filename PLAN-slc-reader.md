# Complexity Assessment: Pure-Julia SLC Reader (NISAR + Sentinel-1)

## Context

SAR.jl aims to provide SAR geometry and metadata handling in pure Julia, with no dependency on
a compiled ISCE3 install. The first deliverable is the SLC metadata reader: **a Julia `SLC` type
covering both NISAR RSLC and Sentinel-1, translated from the ISCE3 ecosystem's Python and C++
sources.**

The reader comes first because it is the prerequisite for everything else. Tracing ISCE3's
`bin/geo2rdr.py` down to its numerics shows the geometry kernel itself is small — roughly 250
lines of scalar Newton–Raphson on a Doppler cost function, with no linear algebra beyond
3-vectors. The cost and the risk live in the input plumbing: versioned NISAR HDF5 layouts and
IPF-dependent Sentinel-1 annotation XML. That plumbing is this document's scope; the geometry
kernel is deferred to a later stage.

Scope decisions:
- **Sentinel-1 source of truth:** `isce-framework/s1-reader`, which produces ISCE3 objects
  natively and so translates with the least friction.
- **Depth:** geometry only — orbit, radar grid, Doppler, burst extents. No radiometric
  calibration, thermal noise, EAP, RFI, or ETAD.
- **Type design:** `abstract type AbstractSLC` with `NisarSLC` and `Sentinel1BurstSLC` subtypes
  sharing one accessor interface.

### Where Sentinel-1 support lives in the ISCE framework

The `isce3` repository proper contains no Sentinel-1 reader — its only matches for "sentinel"
are the programming sense of *sentinel value*, a bibliography entry citing the S1 IPF algorithm
spec, and an ETAD documentation URL. S1 support lives in sibling repositories instead:
**ISCE2** carries a mature TOPS reader (`components/isceobj/Sensor/TOPS/Sentinel1.py`,
`BurstSLC.py`), and **`isce-framework/s1-reader`** is a first-party repo in the same
organization built specifically to feed S1 into ISCE3.

This shapes the port: NISAR translates from `isce3` proper, while Sentinel-1 translates from a
separate repository with a genuinely different data model — bursts and subswaths rather than a
single swath grid.

## Verdict

| Scope | Julia LOC | Effort | Risk |
|---|---|---|---|
| NISAR RSLC reader alone | ~250–300 | 3–5 days | moderate (schema drift) |
| Sentinel-1 geometry-only reader alone | ~450–550 | 1–2 weeks | moderate (XML/IPF variance) |
| Shared type + accessor interface | ~120 | 1–2 days | low |
| **Unified reader, both sensors** | **~850–1000** | **2–3 weeks** | **moderate** |

The two halves share almost nothing at the I/O layer — HDF5 datasets vs. zipped XML — and
converge only at the `AbstractSLC` interface. **Sentinel-1 is roughly 2× the NISAR effort**,
inverting the naive expectation that the NASA-native format is the harder one.

## Type Design

```julia
abstract type AbstractSLC end

struct Orbit
    time::Vector{Float64}                      # seconds since refepoch
    position::Vector{SVector{3,Float64}}       # ECEF metres
    velocity::Vector{SVector{3,Float64}}       # ECEF m/s
    refepoch::DateTime
    interp_method::Symbol                      # :hermite | :legendre
end

struct RadarGrid
    sensing_start::Float64                     # s since refepoch
    prf::Float64                               # Hz
    wavelength::Float64                        # m
    starting_range::Float64                    # m
    range_pixel_spacing::Float64               # m
    length::Int                                # azimuth lines
    width::Int                                 # range samples
    look_side::Symbol                          # :left | :right
    refepoch::DateTime
end
```

`Doppler` is a small union: `ZeroDoppler()` singleton (NISAR `geo2rdr` path, and the honest
default) or `LUT2dDoppler` (regular-grid bilinear). Keeping the zero case a distinct singleton
type — not a zero-filled LUT — keeps the geometry inner loop allocation-free and type-stable.

Required accessor interface, implemented by both subtypes:

```julia
orbit(::AbstractSLC)::Orbit
radar_grid(::AbstractSLC)::RadarGrid
doppler(::AbstractSLC)
wavelength(::AbstractSLC)::Float64
look_side(::AbstractSLC)::Symbol
```

`NisarSLC` additionally carries `filename`, `frequency`, `polarizations`, `root_path`.
`Sentinel1BurstSLC` adds `burst_index`, `subswath`, `polarization`, `tiff_path`,
`first_valid_sample`/`last_valid_sample`/`first_valid_line`/`last_valid_line`,
`azimuth_steer_rate`, `azimuth_fm_rate`, `border`/`center`.

**`azimuth_steer_rate` and `azimuth_fm_rate` are not optional extras for S1** — TOPS mode
sweeps the beam along-track, so any correct S1 geometry or resampling needs them. This is a
real structural asymmetry, not a cosmetic one, and it is the main argument for the abstract-type
design over one flattened struct.

## NISAR Side (~250–300 LOC)

Translate from `isce3`: `nisar/products/readers/SLC/{SLC,RSLC,SLCBase}.py`,
`Base/{Base,Identification}.py`, and the C++ HDF5 layout in `cxx/isce3/core/Serialization.h` +
`cxx/isce3/product/Serialization.h`.

Everything needed, via `HDF5.jl`:

| Quantity | HDF5 path | Notes |
|---|---|---|
| root path | probe `/science/LSAR`, then `/science/SSAR` | `get_hdf5_file_root_path` |
| product group | `{root}/RSLC`, fall back to `{root}/SLC` | early products used `SLC` |
| orbit | `{root}/metadata/orbit/{time,position,velocity}` | position/velocity are N×3 |
| orbit refepoch | **`units` attribute** of `time` | regex out the ISO datetime |
| orbit interp | `{root}/metadata/orbit/interpMethod` | optional, default `"Hermite"` |
| zeroDopplerTime | `{root}/swaths/zeroDopplerTime` | + its `units` attr for refepoch |
| slantRange | `{root}/swaths/frequency{A,B}/slantRange` | |
| az time spacing | `{root}/swaths/zeroDopplerTimeSpacing` | `prf = 1/spacing` |
| centre frequency | `{root}/swaths/.../processedCenterFrequency` | `wavelength = c/f` |
| look side | `{root}/identification/lookDirection` | string "Left"/"Right" |

Derived exactly as ISCE3's `RadarGridParameters(Swath, LookSide)` does:
`starting_range = slantRange[1]`, `range_pixel_spacing = slantRange[2]-slantRange[1]`,
`length = length(zeroDopplerTime)`, `width = length(slantRange)`,
`sensing_start = zeroDopplerTime[1]`.

**Primary risk: documented schema drift.** This is not speculative. `RSLC.py`'s
`getNoiseEquivalentBackscatter` carries an explicit three-version layout table keyed
`v1.2.0` / `v1.0.0` / `v0.0.0`, iterating until a path exists; `Base.py` has two-location
fallbacks for `zeroDopplerTime`/`slantRange`; `ProductPath` probes `RSLC` then `SLC`;
`Identification.py` handles `diagnosticModeFlag` as string, bool-array, *or* integer across
spec generations. Build the same probe-and-fall-back structure in Julia rather than pinning one
layout, and **write it against a real granule, not from the spec.**

**Silent-failure trap:** the reference epoch lives in a `units` *attribute*, not a dataset, and
may be `bytes` or `str` (`RSLC.py` explicitly handles both). Getting it wrong shifts every
azimuth time by a constant with no error raised. Also note ISCE3 never reconciles the orbit's
epoch against the swath's — assert equality (or rebase) instead of inheriting that assumption.

## Sentinel-1 Side (~450–550 LOC)

Translate the geometry subset of `s1-reader`: `s1_reader.py` (1401 LOC) and the geometry fields
of `s1_burst_slc.py` (1192 LOC). The chosen depth excludes most of `s1_annotation.py`
(1338 LOC) — `BurstNoise`, `BurstCalibration`, `BurstEAP`, `AuxCal`, `SwathRfiInfo` are all out
of scope, which is what keeps this tractable.

Needed pieces:

1. **Container access** (~60 LOC). Read either a `.SAFE` directory or the zip directly.
   `ZipArchives.jl` is already in the local depot. Locate `annotation/s1a-iw{N}-slc-{pol}-*.xml`
   while excluding `annotation/calibration/` and `annotation/rfi/` — mirror
   `_is_zip_annotation_xml`.

2. **Annotation XML parsing** (~150 LOC). `EzXML.jl` (also already present). Pull from
   `generalAnnotation/productInformation` (`azimuthSteeringRate` — **note it's converted to
   radians**, `radarFrequency`, `rangeSamplingRate`, `pass`),
   `imageAnnotation/imageInformation` (`azimuthTimeInterval`, `slantRangeTime`,
   `ascendingNodeTime`, `productFirst/LastLineUtcTime`), `swathTiming`
   (`linesPerBurst`, `samplesPerBurst`, `burstList`), and `adsHeader/absoluteOrbitNumber`.

   Derived: `wavelength = c/radar_freq`, `starting_range = slant_range_time*c/2`,
   `range_pixel_spacing = c/(2*range_sampling_rate)`.

3. **Polynomial lists** (~80 LOC). `azimuthFmRateList` and `dopplerCentroid/dcEstimateList`,
   each a time-tagged `Poly1d(coeffs, r0, half_c)`. Then `get_nearest_polynomial` selects by
   burst mid-time. **Format trap, explicitly flagged in the source:** the FM-rate polynomial
   layout changed somewhere between IPF 2.36 and 2.82 — newer files have a named
   `azimuthFmRatePolynomial` element with space-separated text, older ones use positional
   children `elem[2:]`. Must handle both.

4. **Orbit from external EOF** (~90 LOC). S1 orbits are *not* in the product. Parse
   `Data_Block/List_of_OSVs` from a separate POEORB/RESORB `.EOF` XML, filter to
   `[sensing_start - pad, sensing_stop + pad]`, and note the OSV timestamps need the leading
   4 chars stripped (`osv[1].text[4:]` — a `UTC=` prefix). Optionally merge multiple EOFs.
   **`s1-reader` sets `ref_epoch = sensing_start - 2 days`** and `doppler_poly1d_to_lut2d`
   hardcodes a matching `offset_ref_epoch = 2*24*3600`; these two must stay consistent or every
   azimuth time is off by exactly two days.

5. **Per-burst assembly** (~100 LOC). Loop `swathTiming/burstList`. Each burst gets its own
   `azimuthTime` (→ `sensing_start`), `sensingTime`, and valid-sample window. The
   valid-line logic is fiddly and worth transcribing exactly:
   `first_valid_line = first index where firstValidSample >= 0`,
   `n_valid_lines = count of firstValidSample >= 0`, then
   `first_valid_sample = max(fvs[first], fvs[last])` and `last_sample = min(lvs[first], lvs[last])`.

6. **Doppler as LUT2d** (~40 LOC). `doppler_poly1d_to_lut2d` evaluates the range polynomial
   across all samples and `vstack`s it into a 2-row LUT — no azimuth dependence, but a 2D LUT
   is required by the interface. Straightforward, just preserve the epoch offset above.

**Deliberately dropped from scope**, and worth stating so the omission is visible rather than
silent: `S1BurstId` / ESA burst-ID scheme (144 LOC + a track/burst lookup table),
`get_ascending_node_time_orbit` ANX refinement, burst centre/border polygons
(`get_burst_centers_and_boundaries`), `iw2_mid_range` (only used by downstream calibration),
and the entire correction stack. If ESA-compatible burst IDs are needed later for
cross-referencing or stack organisation, add `s1_burst_id.py` — it's self-contained.

## Recommended Approach

Add as a module in the existing `NISAR.jl` (already depends on `HDF5`, `GDAL`, `Proj`,
`Rasters`; verified no radar-geometry code present, so greenfield) or a new `SLCReaders.jl`.
Local depot already has `EzXML`, `XMLDict`, `ZipArchives`, `ZipFile`, `HDF5`. `AutoRIFT.jl`
(7411 LOC, pure Julia, extension-based optional deps) is the right in-house style reference.

Files, buildable and testable in this order:

1. `types.jl` — `AbstractSLC`, `Orbit`, `RadarGrid`, `ZeroDoppler`/`LUT2dDoppler`, accessors.
2. `epoch.jl` — the `units`-attribute / ISO-8601 datetime parsing shared by both sensors.
   Small, and the single highest-leverage correctness item in the whole reader.
3. `nisar.jl` — `NisarSLC` + `read_nisar_slc(path; frequency='A')`, with probe/fallback chains.
4. `sentinel1.jl` — `Sentinel1BurstSLC` + `load_bursts(path, orbit_path, subswath; pol=:vv)`.
5. `orbit_eof.jl` — POEORB/RESORB parsing, kept separate since it's S1-specific and reusable.

Per CLAUDE.md: `SVector{3,Float64}` throughout for type stability, descriptive geospatial names
(`lon`, `lat`, `slant_range`, `az_time` — not `i`/`j`), explicit type annotations on numeric
signatures, and document CRS/epoch conventions in comments since they aren't evident from code.
Prefer `Dates.DateTime` for epochs and plain `Float64` seconds-since-epoch on the hot path.

## Verification

- **Round-trip invariants, no reference implementation needed.** Assert
  `slantRange[end] == starting_range + (width-1)*range_pixel_spacing`; assert
  `zeroDopplerTime` spacing matches `1/prf` to tolerance; assert orbit times are monotonic and
  bracket the swath's sensing window with padding.
- **Epoch cross-check.** Independently derive sensing start from `identification/zeroDopplerStartTime`
  (NISAR) or `productFirstLineUtcTime` (S1) and confirm it agrees with
  `refepoch + sensing_start`. This is the specific test that catches the two-day-offset and
  `units`-attribute bugs, which are otherwise silent.
- **Golden-value comparison against Python, needs the reference stack once.** On a machine with
  ISCE3 + `s1-reader`, dump every `RadarGrid` and `Orbit` field to JSON for a real RSLC and a
  real S1 SAFE, then assert the Julia reader matches to `1e-9` relative. Cheap to produce,
  and it is the only test that catches a systematically misread field that still passes the
  internal consistency checks above.
- **Cross-sensor interface test.** A single generic function taking `::AbstractSLC` and printing
  `orbit`/`radar_grid`/`wavelength`/`look_side` must run unmodified on both subtypes. Add
  `@inferred` assertions on the accessors to lock in type stability.
- **Schema-variance coverage.** Test NISAR against both an `RSLC`- and an `SLC`-group granule if
  available, and S1 against pre- and post-IPF-2.82 products, since both fallback paths exist
  precisely because real files differ. Note in the test file which variants were *not* obtainable
  rather than leaving the gap implicit.

## Bottom Line

~850–1000 LOC, 2–3 weeks, for a unified reader over both sensors at geometry-only depth. The
work is dominated by format archaeology — versioned NISAR HDF5 layouts and IPF-dependent S1
XML — not by algorithms; there is essentially no math here. Sentinel-1 is the larger half at
roughly 2:1, because orbits live in separate EOF files, the container is zipped XML, the data
model is per-burst, and TOPS carries extra geometry state (`azimuth_steer_rate`,
`azimuth_fm_rate`) that NISAR simply doesn't have.

The riskiest items are all silent-failure modes around time: the NISAR `units`-attribute epoch,
the S1 hardcoded two-day epoch offset, and the unreconciled orbit-vs-swath epochs that ISCE3
itself assumes away. Test those explicitly and first.
