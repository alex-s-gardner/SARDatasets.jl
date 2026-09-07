# Reading a remote product's metadata without transferring the product.
#
# A NISAR RSLC is tens of gigabytes, essentially all of it image samples. Its metadata is tens of
# kilobytes, and — measured by logging every range libhdf5 asks for — lies within the first few
# megabytes: the granule this package's fixture comes from needs 26 KiB spread over the first 4.79 MB.
#
# So the head of the file is fetched with one ranged request into a *sparse* local file truncated to
# the product's full length, and that file is opened normally. The tail is never fetched and never
# occupies disk; the filesystem reports holes as zeros, which the metadata reads never reach.
#
# Two things this deliberately does not do. It does not implement an HDF5 virtual file driver: HDF5.jl
# exposes no Julia-callable one, and `h5open` takes a path or a byte image rather than an `IO`. It does
# not guess whether the fetched window was enough — after opening, `open_slc` reads every dataset it
# needs, and a read that falls in the hole fails. `RemoteHTTP` turns that failure into a message naming
# the knob to raise, because a silently truncated read would be worse than a slow one.

"""
    DEFAULT_PREFETCH

Bytes of a remote product fetched up front, 8 MiB.

Large enough for every NISAR product measured, and small enough to be one request. Raise it with
`RemoteHTTP(url; prefetch = ...)` if a product's metadata lies further in.
"""
const DEFAULT_PREFETCH = 8 * 1024 * 1024

"""
    RemoteHTTP(url; prefetch = DEFAULT_PREFETCH, dir = nothing)

A product reachable over HTTP from a server honoring range requests.

Only the first `prefetch` bytes are fetched. `dir` is where the sparse local file goes, defaulting to
a temporary directory; pass one to keep it across sessions.

Earthdata URLs are authenticated from `~/.netrc` and the redirect to the signed data URL is followed,
so a `urs.earthdata.nasa.gov` entry there is all the setup needed.
"""
struct RemoteHTTP <: AbstractSLCSource
    url::String
    prefetch::Int
    dir::Union{Nothing,String}
end

function RemoteHTTP(url::AbstractString; prefetch::Integer = DEFAULT_PREFETCH, dir = nothing)
    prefetch > 0 || throw(ArgumentError("prefetch must be positive, got $prefetch bytes"))
    return RemoteHTTP(String(url), Int(prefetch), dir === nothing ? nothing : String(dir))
end

"""
    RemoteS3(uri; prefetch = DEFAULT_PREFETCH, dir = nothing)

A product in S3, addressed as `s3://bucket/key`.

Reads the same way [`RemoteHTTP`](@ref) does, through the bucket's HTTPS endpoint, and so needs
credentials in the environment (`AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, and
`AWS_SESSION_TOKEN` for the temporary credentials a DAAC issues) and to be running in the bucket's
region. Requester-pays and cross-region access are refused by the endpoint, not here.
"""
struct RemoteS3 <: AbstractSLCSource
    uri::String
    prefetch::Int
    dir::Union{Nothing,String}
end

function RemoteS3(uri::AbstractString; prefetch::Integer = DEFAULT_PREFETCH, dir = nothing)
    startswith(uri, "s3://") || throw(ArgumentError(
        "an S3 source must be addressed as `s3://bucket/key`, got `$uri`"))
    prefetch > 0 || throw(ArgumentError("prefetch must be positive, got $prefetch bytes"))
    return RemoteS3(String(uri), Int(prefetch), dir === nothing ? nothing : String(dir))
end

"""
    s3_https_url(uri; region = nothing) -> String

The virtual-hosted HTTPS endpoint of an `s3://bucket/key` URI.

`region` defaults to `AWS_REGION`, then `AWS_DEFAULT_REGION`. A bucket reached through the wrong
region's endpoint is refused by S3 with a redirect rather than served, so the region has to be right.
"""
function s3_https_url(uri::AbstractString; region = nothing)
    rest = uri[(length("s3://") + 1):end]
    slash = findfirst('/', rest)
    slash === nothing && throw(ArgumentError(
        "`$uri` names a bucket but no key; expected `s3://bucket/key`"))
    bucket = rest[1:(slash - 1)]
    key = rest[(slash + 1):end]
    isempty(bucket) && throw(ArgumentError("`$uri` has an empty bucket name"))
    isempty(key) && throw(ArgumentError("`$uri` has an empty key"))
    reg = region !== nothing ? String(region) :
          get(ENV, "AWS_REGION", get(ENV, "AWS_DEFAULT_REGION", ""))
    isempty(reg) && throw(ArgumentError(
        "reading `$uri` needs a region: set AWS_REGION, or pass `region` to `s3_https_url`"))
    return "https://$bucket.s3.$reg.amazonaws.com/$key"
end

# `curl` is used rather than Downloads.jl because this needs two things Downloads does not expose
# together: a `Range` header on a request that follows redirects, and `.netrc` credentials that are
# *not* replayed to the redirect target. Earthdata redirects to a pre-signed URL which rejects an
# Authorization header, so the credentials must be dropped at the hop — `--location-trusted` would
# keep them and fail.
function _curl_range(url::AbstractString, first_byte::Integer, last_byte::Integer, dest;
                     netrc::Bool = true, headers = String[])
    cookies = tempname()
    args = ["--silent", "--show-error", "--location", "--fail",
            "--range", "$first_byte-$last_byte",
            "--cookie", cookies, "--cookie-jar", cookies,
            "--write-out", "%{http_code} %{size_download}"]
    netrc && push!(args, "--netrc")
    for h in headers
        append!(args, ["--header", h])
    end
    out = IOBuffer()
    err = IOBuffer()
    cmd = pipeline(`curl $args --output $dest $url`; stdout = out, stderr = err)
    try
        run(cmd)
    catch e
        rm(cookies; force = true)
        throw(ErrorException(
            "fetching bytes $first_byte-$last_byte of `$url` failed: " *
            strip(String(take!(err))) * " (curl exited nonzero)"))
    end
    rm(cookies; force = true)
    fields = split(strip(String(take!(out))))
    code = length(fields) >= 1 ? fields[1] : "?"
    got = length(fields) >= 2 ? parse(Int, fields[2]) : 0
    return code, got
end

# The total length comes from `Content-Range` on a one-byte request, not from a `HEAD`: the CDN
# Earthdata redirects to answers HEAD with `Content-Length: 0`, which would size the sparse file to
# nothing.
function _remote_size(url::AbstractString; netrc::Bool = true, headers = String[])
    cookies = tempname()
    args = ["--silent", "--show-error", "--location", "--fail", "--dump-header", "-",
            "--range", "0-0", "--output", "/dev/null",
            "--cookie", cookies, "--cookie-jar", cookies]
    netrc && push!(args, "--netrc")
    for h in headers
        append!(args, ["--header", h])
    end
    hdrs = try
        read(`curl $args $url`, String)
    catch
        rm(cookies; force = true)
        throw(ErrorException(
            "cannot reach `$url`. Check the URL, and for an Earthdata host that `~/.netrc` has a " *
            "`machine urs.earthdata.nasa.gov` entry"))
    end
    rm(cookies; force = true)
    m = match(r"(?im)^content-range:\s*bytes\s+\d+-\d+/(\d+)"m, hdrs)
    m === nothing && throw(ErrorException(
        "`$url` did not answer a range request with a Content-Range header, so its length is " *
        "unknown and it cannot be read without transferring it whole. The server must support " *
        "byte ranges"))
    # The group is not optional, so a match always captured it.
    return parse(Int, m[1]::AbstractString)
end

# A sparse file: the head holds real bytes, the rest is a hole the filesystem reports as zeros. Only
# the head occupies disk, and only the head was transferred.
function _write_sparse(dest::AbstractString, head::Vector{UInt8}, total::Integer)
    open(dest, "w") do io
        write(io, head)
        truncate(io, total)
    end
    return dest
end

function _materialize(url::AbstractString, prefetch::Int, dir::Union{Nothing,String},
                      label::AbstractString; netrc::Bool = true, headers = String[])
    total = _remote_size(url; netrc, headers)
    want = min(prefetch, total)
    into = dir === nothing ? mktempdir(; prefix = "sar_") : dir
    isdir(into) || mkpath(into)
    dest = joinpath(into, basename(split(label, '?')[1]))
    if isfile(dest) && filesize(dest) == total
        # A previous run of the same prefetch left this behind; reuse it rather than refetch.
        return dest
    end
    chunk = tempname()
    code, got = _curl_range(url, 0, want - 1, chunk; netrc, headers)
    head = read(chunk)
    rm(chunk; force = true)
    length(head) == want || throw(ErrorException(
        "requested the first $want bytes of `$url` but received $(length(head)) (HTTP $code); the " *
        "server may not honor byte ranges"))
    return _write_sparse(dest, head, total)
end

localpath(src::RemoteHTTP) = _materialize(src.url, src.prefetch, src.dir, src.url)

function localpath(src::RemoteS3)
    url = s3_https_url(src.uri)
    # S3 authenticates with SigV4 headers, which are not `.netrc` credentials; an anonymous or
    # role-based read carries none. Signing is left to the environment so no credential handling
    # lives here.
    return _materialize(url, src.prefetch, src.dir, src.uri; netrc = false)
end

# What a truncating prefetch looks like when it was too small: HDF5 reads into the hole and finds
# zeros where a group or dataset header should be. The error names the knob rather than leaving a
# caller to infer it from an HDF5 stack trace.
function _prefetch_error(src::Union{RemoteHTTP,RemoteS3}, path::AbstractString, e)
    return ErrorException(
        "read the first $(src.prefetch) bytes of `$(_source_label(src))` but its metadata does not " *
        "lie within them, so the product cannot be opened. Retry with a larger prefetch, e.g. " *
        "`prefetch = $(2 * src.prefetch)`. The underlying failure was: " *
        sprint(showerror, e))
end

_source_label(src::RemoteHTTP) = src.url
_source_label(src::RemoteS3) = src.uri
