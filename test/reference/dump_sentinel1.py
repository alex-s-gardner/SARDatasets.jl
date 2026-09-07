#!/usr/bin/env python
"""Dump Sentinel-1 geometry from isce3 + s1reader as the golden values the Julia reader is tested against.

Mirrors hyp3-autorift's `vend/testGeogrid.py`: `loadMetadataSlc` for the subswath mosaic and
`loadMetadata` for a single burst. Every float is written with its hexadecimal literal beside its
decimal so the comparison is bit-exact rather than to a printed precision.

Run in an environment carrying isce3 and s1reader:

    micromamba run -p /tmp/s1ref python dump_sentinel1.py \
        --safe S1A_..._105E.zip --orbit S1A_OPER_AUX_POEORB_..._005943.EOF \
        --out sentinel1_metadata.json
"""

import argparse
import json
import os
from datetime import timedelta

import isce3
from s1reader import load_bursts

# `loadMetadata` and `loadMetadataSlc` both hardcode this, since Sentinel-1 is a right-looking mission.
LOOK_SIDE = 'Right'


def fx(value):
    """A float as both its decimal and its hexadecimal literal, so a reader can compare bitwise."""
    v = float(value)
    return {'dec': v, 'hex': v.hex()}


def seconds_since(epoch, when):
    """Seconds from an isce3 reference epoch to a datetime, the way the reference computes `aztime`."""
    return float((isce3.core.DateTime(when) - epoch).total_seconds())


def orbit_record(orbit):
    return {
        'n': len(orbit.time),
        'epoch': str(orbit.reference_epoch),
        'time': [fx(t) for t in orbit.time],
        'position': [[fx(c) for c in sv] for sv in orbit.position],
        'velocity': [[fx(c) for c in sv] for sv in orbit.velocity],
    }


def burst_geometry(burst, orbit):
    """One burst's geometry, as `loadMetadata` derives it."""
    prf = 1 / burst.azimuth_time_interval
    nlines, nsamples = burst.shape
    starting_range = burst.starting_range
    far_range = starting_range + (nsamples - 1.0) * burst.range_pixel_spacing
    sensing_start = seconds_since(orbit.reference_epoch, burst.sensing_start)
    return {
        'subswath': burst.burst_id.subswath,
        'polarization': burst.polarization,
        # `sensing_start` is an offset from this epoch, which is the subswath's first burst less two
        # days and so differs from subswath to subswath.
        'epoch': str(orbit.reference_epoch),
        'starting_range': fx(starting_range),
        'far_range': fx(far_range),
        'range_pixel_spacing': fx(burst.range_pixel_spacing),
        'wavelength': fx(burst.wavelength),
        'prf': fx(prf),
        'sensing_start': fx(sensing_start),
        # The reference advances by whole lines rather than reading a stop time from the product.
        'sensing_stop': fx(sensing_start + (nlines - 1.0) / prf),
        'nlines': nlines,
        'nsamples': nsamples,
        'look_side': LOOK_SIDE,
        'sensing_start_utc': burst.sensing_start.isoformat(),
        'azimuth_time_interval': fx(burst.azimuth_time_interval),
    }


def mosaic_geometry(safe, orbit_path, swaths, pol):
    """The geometry spanning several subswaths, as `loadMetadataSlc` derives it."""
    by_swath = {s: load_bursts(safe, orbit_path, s, pol) for s in swaths}
    # The epoch and the near-range parameters come from the lowest-numbered subswath, and the
    # sensing window is the union across all of them.
    first = by_swath[min(swaths)]
    orbit = first[0].orbit

    prf = 1 / first[0].azimuth_time_interval
    starting_range = first[0].starting_range
    range_pixel_spacing = first[0].range_pixel_spacing

    sensing_start_utc = first[0].sensing_start
    sensing_stop_utc = None
    for swath in swaths:
        bursts = by_swath[swath]
        dt = first[0].azimuth_time_interval
        stop = bursts[-1].sensing_start + timedelta(seconds=(bursts[-1].shape[0] - 1) * dt)
        sensing_start_utc = min(sensing_start_utc, bursts[0].sensing_start)
        sensing_stop_utc = stop if sensing_stop_utc is None else max(sensing_stop_utc, stop)

    # Width spans the subswaths: the far-range swath's own width plus its range offset from the near one.
    last = by_swath[max(swaths)]
    nsamples = round((last[-1].starting_range - starting_range) / range_pixel_spacing) + last[-1].shape[1]

    far_range = starting_range + (nsamples - 1.0) * range_pixel_spacing
    sensing_start = seconds_since(orbit.reference_epoch, sensing_start_utc)
    sensing_stop = seconds_since(orbit.reference_epoch, sensing_stop_utc)
    nlines = round((sensing_stop_utc - sensing_start_utc).total_seconds() * prf) + 1

    return {
        'swaths': list(swaths),
        'polarization': pol,
        'starting_range': fx(starting_range),
        'far_range': fx(far_range),
        'range_pixel_spacing': fx(range_pixel_spacing),
        'wavelength': fx(first[0].wavelength),
        'prf': fx(prf),
        'sensing_start': fx(sensing_start),
        'sensing_stop': fx(sensing_stop),
        'nlines': nlines,
        'nsamples': nsamples,
        'look_side': LOOK_SIDE,
        'sensing_start_utc': sensing_start_utc.isoformat(),
        'sensing_stop_utc': sensing_stop_utc.isoformat(),
    }, orbit, by_swath


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--safe', required=True)
    p.add_argument('--orbit', required=True)
    p.add_argument('--pol', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()

    swaths = [1, 2, 3]
    mosaic, orbit, by_swath = mosaic_geometry(args.safe, args.orbit, swaths, args.pol)

    record = {
        # Base names only: the golden values travel with the repository but the products do not, so a
        # test identifies the granule by name and locates it under its own data directory.
        'safe': os.path.basename(args.safe),
        'orbit_file': os.path.basename(args.orbit),
        'polarization': args.pol,
        'isce3_version': isce3.__version__,
        'mosaic': mosaic,
        'orbit': orbit_record(orbit),
        # Every burst of every subswath, so the per-burst path is checked in full rather than at
        # whichever index the reference's own loop happens to land on.
        'bursts': {
            str(s): [burst_geometry(b, by_swath[s][0].orbit) for b in by_swath[s]] for s in swaths
        },
    }

    with open(args.out, 'w') as f:
        json.dump(record, f, indent=1)
    print(f'wrote {args.out}')


if __name__ == '__main__':
    main()
