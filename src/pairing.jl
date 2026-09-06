# Relating an acquisition's clock to another scale, and two acquisitions to each other.
#
# Both are properties of the products, so they live here. Turning a pair into a geometry package's own
# types does not — those are its types, so it extends its own constructors over a `SLC`. Load
# ImagePairGeometry alongside this package to get that.

"""
    repeat_interval(reference::SLC, secondary::SLC) -> Float64

Seconds from the reference acquisition's start to the secondary's.

Measured between the two epochs as well as the two azimuth times, so a pair of products carrying
different epochs — which two acquisitions weeks apart always do — is differenced on one scale rather
than silently on two.

# Examples

```julia
repeat_interval(a, b) / 86400   # the interval in days
```
"""
function repeat_interval(reference::AbstractSLC, secondary::AbstractSLC)
    a, b = reference.geometry, secondary.geometry
    between_epochs = a.epoch == b.epoch ? 0.0 :
                     Dates.value(Dates.Millisecond(b.epoch - a.epoch)) / 1000
    return between_epochs + b.sensing_start - a.sensing_start
end

"""
    epoch_offset(g::RadarGeometry) -> Float64
    epoch_offset(s::SLC) -> Float64

Seconds from midnight of the epoch's own day to the epoch.

[`RadarGeometry`](@ref) reports its azimuth times against `epoch`. A consumer measuring them against
midnight instead — which is the scale a radar geometry kernel indexes azimuth lines on — needs this
constant between the two.

It is zero for a NISAR product, whose epoch *is* midnight of the acquisition day, but that is a property
of the format rather than a general one.
"""
epoch_offset(g::RadarGeometry) =
    Dates.value(Dates.Millisecond(g.epoch - Dates.DateTime(Dates.Date(g.epoch)))) / 1000
epoch_offset(s::AbstractSLC) = epoch_offset(s.geometry)
