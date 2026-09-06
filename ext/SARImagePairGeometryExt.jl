module SARImagePairGeometryExt

# Turning a read acquisition into the geometry types ImagePairGeometry computes with.
#
# The conversion lives here rather than in ImagePairGeometry because that package deliberately depends
# on no IO stack: its radar path takes a `RadarCoordinate` and asks where the numbers came from.
# Reading a product is this package's job, so the bridge belongs on this side of the boundary.
#
# Three things the conversion has to get right, each of which would produce a plausible-looking wrong
# answer if it were assumed instead of checked:
#
# Uniform spacing. `Orbit` interpolates against a uniform axis and rejects a non-uniform one. Products
# have supplied uniform state vectors everywhere measured, but a product that did not would otherwise
# be reported by `Orbit`'s own error, which names state vector indices rather than the product.
#
# The two clocks. `RadarCoordinate` measures the azimuth index against seconds-since-midnight and
# interpolates the orbit against the orbit's own epoch, carrying `orbit_epoch_offset` between them.
# `RadarGeometry` reports both times against one epoch, so the offset is that epoch's seconds past
# midnight — computed here rather than assumed zero, which it happens to be for NISAR only because its
# epoch *is* midnight of the acquisition day.
#
# Coverage. A solve at a time the orbit does not bracket extrapolates rather than failing, so the
# bracket is checked before the coordinate is built.

using SAR
using SAR: SAR, Radar, RadarGeometry, StateVectors, LookSide, LookLeft, LookRight,
           radar_coordinate, image_pair, repeat_interval
using Dates: DateTime, Date, Millisecond, value
using ImagePairGeometry: ImagePairGeometry, RadarCoordinate, CoregisteredPair, incidence_angle
import ImagePairGeometry

# The enums are distinct types with the same meaning; neither package should import the other's.
_look(side::LookSide) = side == LookLeft ? ImagePairGeometry.LookLeft : ImagePairGeometry.LookRight

# Seconds from midnight of the epoch's own day to the epoch. `RadarCoordinate` wants azimuth times on
# the seconds-since-midnight scale the reference measures the azimuth index against
# (`geogridRadar.cpp:972`), and the orbit on its own scale, so this is the constant between them.
function _epoch_offset(epoch::DateTime)
    midnight = DateTime(Date(epoch))
    return value(Millisecond(epoch - midnight)) / 1000
end

"""
    ImagePairGeometry.Orbit(sv::StateVectors) -> Orbit

The state vectors as an interpolating orbit.

Throws if the times are not uniformly spaced, which the interpolant requires.
"""
function ImagePairGeometry.Orbit(sv::StateVectors)
    n = length(sv.time)
    n >= 4 || throw(ArgumentError(
        "an interpolating orbit needs at least 4 state vectors, but the product supplies $n"))
    return ImagePairGeometry.Orbit(; time = sv.time, position = sv.position,
                                  velocity = sv.velocity)
end

"""
    radar_coordinate(s::Radar; zrange = nothing, chebyshev = false) -> RadarCoordinate

The acquisition as an [`ImagePairGeometry.RadarCoordinate`](@ref).

The scene-center incidence angle is computed here because the type stores it: the reference computes
it before running the geometry, so it is an input rather than something derived on demand.

`chebyshev` swaps the orbit interpolant for `ImagePairGeometry.chebyshev_orbit`, which is faster and
not bitwise identical to the default. `zrange` overrides the elevation pair the incidence angle is
averaged over.

# Examples

```julia
using SAR, ImagePairGeometry
s = open_sar("NISAR_L1_PR_RSLC_....h5")
coord = radar_coordinate(s)
```
"""
function SAR.radar_coordinate(s::Radar; zrange = nothing, chebyshev::Bool = false)
    g = s.geometry
    sv = orbit(s)

    g.epoch == sv.epoch || throw(ArgumentError(
        "the azimuth times are measured against $(g.epoch) but the state vectors against " *
        "$(sv.epoch). Converting between them is not implemented, because every product measured " *
        "puts both on one epoch and a product that does not may differ in more than this"))

    # An out-of-range solve extrapolates rather than failing, so a product whose orbit does not span
    # its own acquisition is refused here.
    first(sv.time) <= g.sensing_start && g.sensing_stop <= last(sv.time) || throw(ArgumentError(
        "the state vectors span $(first(sv.time))–$(last(sv.time)) s but the acquisition runs " *
        "$(g.sensing_start)–$(g.sensing_stop) s, so the orbit does not cover it"))

    orb = ImagePairGeometry.Orbit(sv)
    chebyshev && (orb = ImagePairGeometry.chebyshev_orbit(orb))

    offset = _epoch_offset(g.epoch)
    side = _look(g.look_side)
    kwargs = (; orbit = orb, starting_range = g.starting_range, dr = g.range_pixel_spacing,
              sensing_start = g.sensing_start, prf = g.prf, nsamples = g.nsamples,
              nlines = g.nlines, look_side = side, wavelength = g.wavelength,
              orbit_epoch_offset = offset)
    ia = zrange === nothing ? incidence_angle(; kwargs...) :
         incidence_angle(; kwargs..., zrange)
    return RadarCoordinate(; kwargs..., incidence_angle = ia)
end

"""
    image_pair(reference::Radar, secondary::Radar; kwargs...) -> CoregisteredPair

The two acquisitions as a pair, for the geometry of the reference and the interval between them.

Radar geometry comes from the reference alone — `testGeogrid.py:427-470` takes every radar parameter
from image 1 and the secondary only for the repeat interval — so `secondary` contributes its sensing
time and nothing else. `kwargs` are forwarded to [`radar_coordinate`](@ref).

# Examples

```julia
pair = image_pair(open_sar(url1), open_sar(url2))
pair.dt / 86400   # the repeat interval in days
```
"""
function SAR.image_pair(reference::Radar, secondary::Radar; kwargs...)
    dt = repeat_interval(reference, secondary)
    dt > 0 || throw(ArgumentError(
        "the secondary acquisition starts $(-dt) s before the reference; pass them in acquisition " *
        "order"))
    return CoregisteredPair(radar_coordinate(reference; kwargs...); dt)
end

end
