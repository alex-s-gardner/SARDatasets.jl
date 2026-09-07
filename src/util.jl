# Every scalar string in a NISAR product is stored as fixed-length bytes, and HDF5.jl surfaces those
# variously as a `String` or as a byte vector depending on the dataset. Trailing NULs are padding.
#
# The return type is annotated because the callers pass the result of an HDF5 read, whose type is only
# known once the file is open: without it every string field of a product would infer as `Any` and each
# use of one would dispatch at runtime.
_string(x::String)::String = String(rstrip(x, '\0'))
_string(x::AbstractString)::String = _string(String(x))
_string(x::AbstractVector{UInt8})::String = _string(String(x))
_string(x::AbstractArray)::String = _string(only(x))
