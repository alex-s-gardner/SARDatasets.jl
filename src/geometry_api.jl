# The geometry conversions, as stubs an extension fills in.
#
# Turning an acquisition into a geometry package's own types needs that package loaded, so the methods
# live in an extension. The functions are declared here so they can be exported and documented from
# one place, and so calling one without the extension loaded says what to load rather than reporting an
# undefined name.

"""
    radar_coordinate(s::Radar; kwargs...)

The acquisition as a geometry package's radar coordinate.

Provided by an extension. Load ImagePairGeometry alongside this package to get it:

```julia
using SARDatasets, ImagePairGeometry
coord = radar_coordinate(open_sar(path))
```
"""
function radar_coordinate(s::AbstractRadar; kwargs...)
    throw(ArgumentError(
        "`radar_coordinate` needs a geometry package loaded to know what to build: " *
        "`using ImagePairGeometry` alongside this one enables it"))
end

"""
    image_pair(reference::Radar, secondary::Radar; kwargs...)

The two acquisitions as a coregistered pair.

Provided by an extension, as [`radar_coordinate`](@ref) is.
"""
function image_pair(reference::AbstractRadar, secondary::AbstractRadar; kwargs...)
    throw(ArgumentError(
        "`image_pair` needs a geometry package loaded to know what to build: " *
        "`using ImagePairGeometry` alongside this one enables it"))
end

"""
    repeat_interval(reference::Radar, secondary::Radar) -> Float64

Seconds from the reference acquisition's start to the secondary's.

Measured between the two epochs as well as the two azimuth times, so a pair of products carrying
different epochs — which two dates always do — is differenced on one scale rather than silently on
two.
"""
function repeat_interval(reference::AbstractRadar, secondary::AbstractRadar)
    a, b = reference.geometry, secondary.geometry
    between_epochs = (a.epoch == b.epoch) ? 0.0 :
                     Dates.value(Dates.Millisecond(b.epoch - a.epoch)) / 1000
    return between_epochs + b.sensing_start - a.sensing_start
end
