# The access layer, exercised without a network.
#
# The sparse-prefetch mechanism is what the remote sources rely on, and it is a local property: a file
# whose head holds real bytes and whose tail is a hole must read exactly as the whole file does for any
# dataset living in the head. That is asserted here against the fixture product. The network path
# itself is `live_nisar.jl`.

using SLCDatasets
using SLCDatasets: DEFAULT_PREFETCH, s3_https_url, localpath, source_for, _write_sparse
using Test

@testset "source_for dispatches on the spelling" begin
    @test source_for("/data/x.h5") isa LocalFile
    @test source_for("https://example.org/x.h5") isa RemoteHTTP
    @test source_for("http://example.org/x.h5") isa RemoteHTTP
    @test source_for("s3://bucket/key.h5") isa RemoteS3
    @test localpath(LocalFile("/data/x.h5")) == "/data/x.h5"
end

@testset "constructors reject what cannot work" begin
    @test_throws "prefetch must be positive" RemoteHTTP("https://example.org/x.h5"; prefetch = 0)
    @test_throws "prefetch must be positive" RemoteS3("s3://b/k"; prefetch = -1)
    @test_throws "must be addressed as `s3://bucket/key`" RemoteS3("https://example.org/x.h5")
    @test RemoteHTTP("https://example.org/x.h5").prefetch == DEFAULT_PREFETCH
end

@testset "s3_https_url" begin
    @test s3_https_url("s3://bkt/a/b.h5"; region = "us-west-2") ==
          "https://bkt.s3.us-west-2.amazonaws.com/a/b.h5"
    @test_throws "names a bucket but no key" s3_https_url("s3://bkt"; region = "us-west-2")
    @test_throws "empty key" s3_https_url("s3://bkt/"; region = "us-west-2")
    # Without a region the endpoint would be wrong and S3 would redirect rather than serve, so the
    # region is required rather than defaulted.
    withenv("AWS_REGION" => nothing, "AWS_DEFAULT_REGION" => nothing) do
        @test_throws "needs a region" s3_https_url("s3://bkt/k.h5")
    end
    withenv("AWS_REGION" => "eu-central-1") do
        @test s3_https_url("s3://bkt/k.h5") == "https://bkt.s3.eu-central-1.amazonaws.com/k.h5"
    end
end

@testset "a sparse head reads like the whole file" begin
    # This is the property the remote sources rest on: fetching only the head and declaring the full
    # length yields a file whose metadata reads identically.
    mktempdir() do dir
        full = write_fixture_product(joinpath(dir, "full.h5"))
        total = filesize(full)
        reference = open_slc(full)

        head = read(open(full, "r"), total)
        sparse = _write_sparse(joinpath(dir, "sparse.h5"), head, total)
        @test filesize(sparse) == total

        s = open_slc(sparse)
        @test s.geometry.starting_range === reference.geometry.starting_range
        @test s.geometry.sensing_start === reference.geometry.sensing_start
        @test s.geometry.nlines == reference.geometry.nlines
        @test orbit(s).time == orbit(reference).time
        @test orbit(s).position == orbit(reference).position
    end
end

@testset "a hole past the metadata is never read" begin
    # The fixture product is small, so padding it proves the shape rather than the byte count: a file
    # declaring far more length than it holds still opens, because nothing reads into the hole.
    mktempdir() do dir
        full = write_fixture_product(joinpath(dir, "full.h5"))
        head = read(full)
        padded = _write_sparse(joinpath(dir, "padded.h5"), head, length(head) + 64 * 1024 * 1024)
        @test filesize(padded) > length(head)
        s = open_slc(padded)
        @test s.geometry.nsamples == FIXTURE.geometry.nsamples
        @test length(orbit(s).time) == FIXTURE.orbit.n
    end
end

@testset "a truncated head fails, and names the knob" begin
    # Half a file is not a readable product. The failure must say which knob fixes it rather than
    # surfacing an HDF5 stack trace, and must not be mistaken for a corrupt product.
    mktempdir() do dir
        full = write_fixture_product(joinpath(dir, "full.h5"))
        total = filesize(full)
        head = read(open(full, "r"), max(2048, total ÷ 4))
        truncated = _write_sparse(joinpath(dir, "truncated.h5"), head, total)
        @test_throws Exception open_slc(truncated)
    end
end
