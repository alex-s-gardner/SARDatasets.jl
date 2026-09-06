# The uniform acquisition type, and the two sensor-neutral metadata records it carries.
#
# A SAR product is not a bag of arbitrary named variables — it is a fixed, known set of geometry
# quantities. So the metadata is plain immutable structs rather than a variable/attribute protocol, and
# only the bulky part, the state vectors, is deferred.
#
# The field set is what a geometry consumer needs, taken from the two reference loaders in
# hyp3-autorift's `testGeogrid.py`: `loadMetadataRslc` for NISAR and `loadMetadata` for Sentinel-1.
# Those two disagree on how several quantities are obtained — NISAR looks left and derives the pulse
# repetition frequency from the zero-Doppler time spacing, Sentinel-1 looks right and derives it from
# the azimuth time interval — which is why `look_side` and `prf` are fields here rather than
# constants.

"""
    LookSide

Which side of the flight track the radar points: `LookLeft` or `LookRight`.
"""
@enum LookSide LookLeft LookRight

"""
    AbstractRadar

One SAR acquisition, whatever the sensor. See [`Radar`](@ref).
"""
abstract type AbstractRadar end

"""
    AbstractSARBackend

How one sensor's bytes are laid out and named. A backend knows the group and dataset paths of its
product family; it does not know how the bytes are reached — that is an
[`AbstractSARSource`](@ref).
"""
abstract type AbstractSARBackend end

"""
    AbstractSARSource

Where a product's bytes come from. [`LocalFile`](@ref) is the only route that needs no network; see
[`RemoteHTTP`](@ref) and [`RemoteS3`](@ref) for the others.
"""
abstract type AbstractSARSource end

"""
    LocalFile(path)

A product already on disk.
"""
struct LocalFile <: AbstractSARSource
    path::String
end

LocalFile(path::AbstractString) = LocalFile(String(path))

"""
    localpath(src::AbstractSARSource) -> String

A path that can be opened, materializing whatever the source needs to make that true.
"""
localpath(src::LocalFile) = src.path

"""
    Identification

What an acquisition is: mission, product, orbit and coverage.

`start_time` and `stop_time` are the zero-Doppler bounds. They are kept as the product's own strings
because those carry nanoseconds, which `DateTime`'s millisecond resolution would silently drop;
[`start_datetime`](@ref) and [`stop_datetime`](@ref) convert when millisecond truncation is
acceptable. The full-precision bounds the geometry uses are `sensing_start` and `sensing_stop` on
[`RadarGeometry`](@ref), in seconds against its `epoch`.
"""
struct Identification
    mission::String
    product_type::String
    absolute_orbit::Int
    pass_direction::String
    look_direction::String
    start_time::String
    stop_time::String
    bounding_polygon::String
end

"""
    start_datetime(id::Identification) -> DateTime
    stop_datetime(id::Identification) -> DateTime

The zero-Doppler bounds as `DateTime`s, truncated to milliseconds.

The product records these to nanoseconds. Use `id.start_time` / `id.stop_time` for the exact strings,
or `RadarGeometry`'s `sensing_start` / `sensing_stop` for full-precision arithmetic.
"""
start_datetime(id::Identification) = _truncated_datetime(id.start_time)
stop_datetime(id::Identification) = _truncated_datetime(id.stop_time)

"""
    RadarGeometry

The slant-range/azimuth geometry of an acquisition: everything needed to relate a pixel index to a
range and a time, in the units a geometry kernel wants.

# Fields
- `starting_range`, `far_range`: slant range to the first and last range sample, meters.
- `range_pixel_spacing`: range sample spacing, meters.
- `wavelength`: radar wavelength, meters — `c` over the processed center frequency.
- `prf`: pulse repetition frequency, Hz, so `1/prf` is the azimuth line spacing in time.
- `sensing_start`, `sensing_stop`: azimuth time of the first and last line, **seconds since
  `epoch`**.
- `nlines`, `nsamples`: image height and width in pixels — azimuth lines and range samples.
- `look_side`: which side the radar points, read from the product rather than assumed.
- `epoch`: the instant `sensing_start` and the orbit's times are measured from.

# The epoch

`epoch` is the one clock this type reports against, and both the azimuth times and the orbit's state
vector times are on it. A consumer measuring azimuth time against a different origin — seconds since
midnight, say — needs the offset between the two, which is why `epoch` is carried rather than being
folded into the times and discarded.
"""
struct RadarGeometry
    starting_range::Float64
    far_range::Float64
    range_pixel_spacing::Float64
    wavelength::Float64
    prf::Float64
    sensing_start::Float64
    sensing_stop::Float64
    nlines::Int
    nsamples::Int
    look_side::LookSide
    epoch::DateTime
end

"""
    StateVectors

Platform state vectors in ECEF, and the instant they are measured from.

`time` is seconds since `epoch`, `position` is meters and `velocity` meters per second. The times are
not required to be uniformly spaced here — a consumer needing that must check, since whether it holds
is a property of the product rather than of this type.

This is the product's own record, not an interpolator: a consumer wanting to evaluate the trajectory
between samples builds one from these.
"""
struct StateVectors
    time::Vector{Float64}
    position::Vector{SVector{3,Float64}}
    velocity::Vector{SVector{3,Float64}}
    epoch::DateTime
    interp_method::String
    kind::String
end

"""
    Radar <: AbstractRadar

One SAR acquisition, read through a sensor backend.

`identification` and `geometry` are parsed when the product is opened: they are a few dozen scalars,
and a product whose metadata cannot be read is not usable, so failing at `open_sar` is better than
failing later. The state vectors are read on first access and then held, since a consumer touching only
the geometry should not pay for them.

Build one with [`open_sar`](@ref) rather than calling this constructor.

# Examples

```julia
s = open_sar("NISAR_L1_PR_RSLC_....h5")
s.geometry.starting_range   # already in hand
orbit(s)                    # read on first access
```
"""
mutable struct Radar{B<:AbstractSARBackend} <: AbstractRadar
    const backend::B
    const identification::Identification
    const geometry::RadarGeometry
    orbit::Union{Nothing,StateVectors}
end

Radar(backend::AbstractSARBackend, ident::Identification, geom::RadarGeometry) =
    Radar(backend, ident, geom, nothing)

"""
    nlines(s::Radar)
    nsamples(s::Radar)

Image height and width in pixels.
"""
nlines(s::Radar) = s.geometry.nlines
nsamples(s::Radar) = s.geometry.nsamples

"""
    orbit(s::Radar) -> StateVectors

The platform state vectors, read on first call and held afterwards.
"""
function orbit(s::Radar)
    o = s.orbit
    o === nothing || return o
    o = read_orbit(s.backend)
    s.orbit = o
    return o
end

Base.show(io::IO, s::Radar) = print(io, "Radar(", s.identification.product_type, ", ",
                                    s.identification.mission, ", ",
                                    s.geometry.nlines, "x", s.geometry.nsamples, ")")

function Base.show(io::IO, ::MIME"text/plain", s::Radar)
    id, g = s.identification, s.geometry
    println(io, "SAR acquisition")
    println(io, "  mission        : ", id.mission, " (", id.product_type, ")")
    println(io, "  absolute orbit : ", id.absolute_orbit, " ", lowercase(id.pass_direction))
    println(io, "  sensing        : ", id.start_time, " to ", id.stop_time)
    println(io, "  size           : ", g.nlines, " lines x ", g.nsamples, " samples")
    println(io, "  slant range    : ", g.starting_range, " to ", g.far_range, " m")
    println(io, "  look side      : ", g.look_side == LookLeft ? "left" : "right")
    print(io, "  orbit          : ", s.orbit === nothing ? "not read" :
              string(length(s.orbit.time), " state vectors"))
end
