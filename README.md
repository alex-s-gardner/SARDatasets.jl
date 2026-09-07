# SLCDatasets

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://alex-s-gardner.github.io/SLCDatasets.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alex-s-gardner.github.io/SLCDatasets.jl/dev/)
[![Build Status](https://github.com/alex-s-gardner/SLCDatasets.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alex-s-gardner/SLCDatasets.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/alex-s-gardner/SLCDatasets.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/alex-s-gardner/SLCDatasets.jl)

Read single-look complex SAR products into one type, whatever the sensor.

```julia
using SLCDatasets

s = open_slc("NISAR_L1_PR_RSLC_....h5")
s.geometry.starting_range   # 895255.2277025
s.geometry.prf              # 1520.0
s.geometry.look_side        # LookLeft
orbit(s)                    # state vectors, read on first access
```

`open_slc` returns a `SLC` carrying an `Identification` — what the acquisition is — and a
`RadarGeometry`: the slant-range/azimuth geometry in the units a geometry kernel wants. Both are parsed
when the product is opened, since a product whose metadata cannot be read is not usable. The state
vectors are deferred, and the sensor is chosen by inspecting the product rather than by the caller
naming it.

The field set is sensor-neutral by construction. `look_side` and `prf` are fields rather than constants
because the two reference loaders they follow disagree on both: NISAR looks left with the pulse
repetition frequency from the zero-Doppler time spacing, Sentinel-1 looks right with it from the azimuth
time interval.

## Reading a remote product without transferring it

A NISAR RSLC is tens of gigabytes, nearly all of it image samples. Its metadata is tens of kilobytes
lying within the first few megabytes — measured by logging every byte range the HDF5 library requests,
one granule needs 26 KiB spread over the first 4.79 MB.

So a remote product is opened by fetching its head with one ranged request into a sparse local file and
reading that. The tail is never fetched and occupies no disk.

```julia
s = open_slc(RemoteHTTP(url))                  # 8 MB of an 11.77 GB granule
s = open_slc(RemoteS3("s3://bucket/key.h5"))   # inside the bucket's region
```

Earthdata URLs authenticate from `~/.netrc` and follow the redirect to the signed data URL, so a
`machine urs.earthdata.nasa.gov` entry is the only setup. A prefetch window too small to hold the
metadata fails with the knob to raise rather than an HDF5 stack trace — a truncated read is never
returned as though it were complete.

## Feeding a geometry kernel

This package reads products and says nothing about geometry. Converting an acquisition into a geometry
package's types belongs to that package, since they are its types:
[ImagePairGeometry.jl](https://github.com/alex-s-gardner/ImagePairGeometry.jl) extends its own
constructors over a `SLC`, so loading both is all that is needed.

```julia
using ImagePairGeometry, SLCDatasets

pair = CoregisteredPair(open_slc(url1), open_slc(url2))
pair.coordinate    # a RadarCoordinate, incidence angle included
pair.dt / 86400    # the repeat interval in days
```

What does live here is what a pair of *products* means: `repeat_interval` spans both epochs, and
`epoch_offset` gives the seconds between a product's epoch and midnight — the constant a kernel indexing
azimuth lines against midnight needs.

Verified against the reference implementations on real granules. For NISAR, against `h5py` and `isce3`:
all 26 geometry and identification values agree bitwise, and the scene-center range, azimuth time,
position and velocity agree to 0 ULP. For Sentinel-1, against `isce3` + `s1reader` on two granules from
different missions — every range, spacing, wavelength, PRF and image dimension of the mosaic and of all
27 bursts agrees bitwise, state vectors agree bitwise, and azimuth times agree to a nanosecond.

`test/reference/dump_sentinel1.py` regenerates the golden values.

## Sentinel-1

An IW SLC is read from a `.SAFE` directory or the zip of one. The state vectors come from a POEORB or
RESORB `.EOF` file, which the product does not contain, so `orbit` is required — an acquisition whose
orbit is unavailable is refused rather than returned half-read.

```julia
safe = "S1A_IW_SLC__1SSH_20151120T080202_....zip"
eof  = "S1A_OPER_AUX_POEORB_OPOD_....EOF"

mosaic = open_slc(safe; orbit = eof)                       # across subswaths 1-3
burst  = open_slc(safe; orbit = eof, swath = 2, burst = 3)  # one burst
nbursts(safe; swath = 2)                                    # 9
```

A TOPS product is three subswaths of bursts at different slant ranges rather than one image, so there
are two defensible geometries and both are available: the mosaic spanning the subswaths, and a single
burst. Only the annotation XML is read, so the measurement TIFFs cost nothing and a zip need not be
unpacked.

To work through a subswath's bursts, `bursts` parses the annotation once and returns them as an
`AbstractVector` of `SLC`s, where a call to `open_slc` per burst would re-read the product each time.

```julia
b = bursts(safe; orbit = eof, swath = 2)   # an AbstractVector{SLC}, one parse
length(b)                                  # 9
b[3].geometry.starting_range               # built on indexing
[nlines(s) for s in b]
```

A `Sentinel1Product` is the parsed product itself, when several subswaths are wanted from one read:

```julia
p = Sentinel1Product(safe; orbit = eof)    # subswaths 1-3, parsed once
nbursts(p, 2)
bursts(p, 3)
```

Times are reported against an `epoch` two days before the anchor burst's sensing start, matching the
offset `s1reader` and ISCE3's Doppler LUTs use; `epoch + sensing_start` is the instant the product
records.

## Scope

SLCs in radar geometry: NISAR-format HDF5 (RSLC) and Sentinel-1 IW SLC.

Other SAR products are out of scope rather than unimplemented. An interferogram (RIFG, RUNW) or a
covariance product (GCOV) carries a multilooked grid, and a geocoded SLC (GSLC) carries map coordinates
rather than a slant-range axis. None is described by `RadarGeometry`, which is the type this package
exists to produce, so a geocoded product is refused by name with that reason.
