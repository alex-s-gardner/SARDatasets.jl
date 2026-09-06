#!/usr/bin/env python
"""Harvest a NISAR product's geometry metadata as a committed fixture.

Reads the datasets `SAR.jl`'s NISAR backend reads, through `h5py` rather than through Julia, and
records every float as both a decimal and a hex literal so the Julia side can assert bit-exact
agreement rather than agreement to some printed precision.

The fixture is a few kilobytes, so the test suite runs with no product, no network and no HDF5
library beyond the one Julia ships.

Run in the reference environment:

    micromamba run -n geogrid-ref python test/reference/gen_nisar_metadata.py <product.h5> [out.json]
"""
import json
import os
import sys

import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
SPEED_OF_LIGHT = 299792458.0


def h(x):
    return {'dec': float(x), 'hex': float(x).hex()}


def s(x):
    return x.decode().rstrip('\x00') if isinstance(x, bytes) else str(x)


def band_of(f):
    for b in ('LSAR', 'SSAR'):
        if b in f['science']:
            return b
    raise SystemExit('no LSAR or SSAR group')


def harvest(path):
    f = h5py.File(path, 'r')
    band = band_of(f)
    ident = f[f'science/{band}/identification']
    declared = s(ident['productType'][()])
    product = 'RSLC' if declared == 'SLC' and 'RSLC' in f[f'science/{band}'] else declared
    freq = s(ident['listOfFrequencies'][:][0])

    p = f'science/{band}/{product}'
    sw = f[f'{p}/swaths']
    fr = f[f'{p}/swaths/frequency{freq}']
    orb = f[f'{p}/metadata/orbit']

    sr = fr['slantRange']
    zdt = sw['zeroDopplerTime']
    zdt_spacing = float(sw['zeroDopplerTimeSpacing'][()])
    fc = float(fr['processedCenterFrequency'][()])

    t = orb['time'][:]
    pos = orb['position'][:]
    vel = orb['velocity'][:]

    return {
        'source': os.path.basename(path),
        'band': band,
        'product_type': product,
        'frequency': freq,
        'identification': {
            'mission': s(ident['missionId'][()]),
            'product_type': declared,
            'absolute_orbit': int(ident['absoluteOrbitNumber'][()]),
            'pass_direction': s(ident['orbitPassDirection'][()]),
            'look_direction': s(ident['lookDirection'][()]),
            'start_time': s(ident['zeroDopplerStartTime'][()]),
            'stop_time': s(ident['zeroDopplerEndTime'][()]),
        },
        'geometry': {
            'starting_range': h(sr[0]),
            'far_range': h(sr[-1]),
            'range_pixel_spacing': h(fr['slantRangeSpacing'][()]),
            'wavelength': h(SPEED_OF_LIGHT / fc),
            'prf': h(1.0 / zdt_spacing),
            'sensing_start': h(zdt[0]),
            'sensing_stop': h(zdt[-1]),
            'nlines': int(zdt.shape[0]),
            'nsamples': int(sr.shape[0]),
            'epoch': s(zdt.attrs['units']),
        },
        'orbit': {
            'n': int(t.shape[0]),
            'epoch': s(orb['time'].attrs['units']),
            'interp_method': s(orb['interpMethod'][()]) if 'interpMethod' in orb else 'Hermite',
            'kind': s(orb['orbitType'][()]) if 'orbitType' in orb else 'Custom',
            'time': [h(v) for v in t],
            'position': [[h(v) for v in row] for row in pos],
            'velocity': [[h(v) for v in row] for row in vel],
        },
        'versions': {'h5py': h5py.version.version, 'hdf5': h5py.version.hdf5_version},
    }


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'nisar_metadata.json')
    json.dump(harvest(sys.argv[1]), open(out, 'w'), indent=1)
    print('wrote', out)
