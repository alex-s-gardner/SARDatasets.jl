# Whether an acquisition's samples carry a TOPS azimuth ramp.
#
# This is the question a consumer of `pixels` has to ask before interpolating them, and it lives here
# because the answer is a property of how the product was collected — which this package reads and a
# geometry or correlation package cannot infer.
#
# What makes it worth testing rather than obvious: the answer has to be right for four different backend
# kinds, and one of them does not inherit it. `MergedBurstBackend` describes several bursts and subtypes
# `AbstractSLCBackend` rather than `AbstractBurstBackend`, so it would default to `false` without a method
# of its own. A merge silently reported as non-TOPS is exactly the case that would reach a resampler.

using SLCDatasets
using SLCDatasets: AbstractBurstBackend, AbstractSLCBackend, MergedBurstBackend, NisarBackend,
                   Sentinel1Backend, is_tops, deramp_parameters
using Test

@testset "a NISAR acquisition is not TOPS" begin
    # Stripmap at zero Doppler, so the samples are directly interpolable and a consumer needs no
    # deramping. Read from the committed fixture, so this needs no granule.
    mktempdir() do dir
        s = open_slc(write_fixture_product(joinpath(dir, "a.h5")))
        @test !is_tops(s)
        @test !is_tops(s.backend)
        @test s.backend isa NisarBackend
    end
end

@testset "every Sentinel-1 form is TOPS" begin
    # A mosaic, a single burst, and a merge of bursts. All three describe TOPS data, and the three reach
    # the answer by different routes — the first two through `AbstractBurstBackend`, the third through its
    # own method — so all three are asserted rather than one standing for the family.
    if !isempty(S1_PRODUCTS)
        for p in S1_PRODUCTS
            pol = String(p.gold.polarization)
            mosaic = open_slc(p.safe; orbit = p.eof, swath = 2, polarization = pol)
            @test is_tops(mosaic)
            @test mosaic.backend isa Sentinel1Backend
            # The type hierarchy is what carries the answer for a burst, so that a backend added later
            # inherits it rather than having to remember to declare it.
            @test mosaic.backend isa AbstractBurstBackend

            b = bursts(p.safe; orbit = p.eof, swath = 2, polarization = pol)
            @test is_tops(b[1])
            @test is_tops(b[end])
            @test all(is_tops, b)

            # The case that does not inherit: a merge subtypes `AbstractSLCBackend` directly.
            m = merge_bursts(b)
            @test m.backend isa MergedBurstBackend
            @test !(m.backend isa AbstractBurstBackend)
            @test is_tops(m)
            # A merge of one burst too, since that is the degenerate placement and takes a different path
            # through the grid.
            @test is_tops(merge_bursts(b[1:1]))
        end
    else
        @info "skipping the Sentinel-1 TOPS cases; no committed products resolved"
    end
end

@testset "a backend that does not say is not TOPS" begin
    # The default is `false` rather than an error, so that a non-TOPS sensor added later needs no negative
    # declaration. Checked against a backend type declared here, which the package has never seen — if
    # this needed a method the default would not be doing its job.
    struct _PlainBackend <: AbstractSLCBackend end
    @test !is_tops(_PlainBackend())
end

@testset "the deramp is refused, naming what it needs" begin
    # Not implemented, and it says which annotation fields are missing and that they are in the annotation
    # already read. A consumer that needs the deramp is told what to add rather than discovering the gap
    # as a phase error in a resampled image.
    mktempdir() do dir
        s = open_slc(write_fixture_product(joinpath(dir, "a.h5")))
        err = try
            deramp_parameters(s)
            nothing
        catch e
            sprint(showerror, e)
        end
        @test err !== nothing
        @test occursin("not implemented", err)
        # The three fields, by their annotation names, so the message is actionable.
        @test occursin("azimuthFmRateList", err)
        @test occursin("dcEstimateList", err)
        @test occursin("azimuthSteeringRate", err)
        # And the one thing that remains safe, since a consumer doing amplitude tracking on Sentinel-1
        # has a legitimate use for the samples today.
        @test occursin("amplitude", err)
    end
end
