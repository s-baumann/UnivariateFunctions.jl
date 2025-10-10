using Test

@testset "Isotonic Regression Tests" begin
    using UnivariateFunctions
    using Dates, Random, DataFrames

    tol = 10*eps()

    obs = 1000
    twister = MersenneTwister(3)
    dd = DataFrame(x = rand(twister, obs) .* 2)
    # decreasing data.
    dd[!, :y] = -1 *  dd.x.^2 .+ rand(twister, Normal(),obs) .+ 1
    iso_increasing = UnivariateFunctions.isotonic_regression(dd, :x, :y; increasing = true)
    plot(iso_increasing, dd)
    @test iso_increasing isa Piecewise_Function
    @test iso_increasing(1.9) >= iso_increasing(0.1)
    iso_decreasing = UnivariateFunctions.isotonic_regression(dd, :x, :y; increasing = false)
    plot(iso_decreasing, dd)
    @test iso_decreasing isa Piecewise_Function
    @test iso_decreasing(1.9) < iso_decreasing(0.1)

    # decreasing data.
    dd[!, :y] = -1 *  dd.x.^2 .+ rand(twister, Normal(),obs) .+ 1
    mono_increasing = UnivariateFunctions.monotonic_regression(dd, :x, :y; increasing = true)
    plot(mono_increasing, dd)
    @test mono_increasing isa Piecewise_Function
    @test mono_increasing(1.9) >= mono_increasing(0.1)
    mono_decreasing = UnivariateFunctions.monotonic_regression(dd, :x, :y; nbins=15, equally_spaced_bins=false, increasing = false)
    plot(mono_decreasing, dd)
    @test mono_decreasing isa Piecewise_Function
    @test mono_decreasing(1.9) < mono_decreasing(0.1)


    # increasing data.
    dd = DataFrame(x = rand(twister, obs) .* 2)
    dd[!, :y] = 1 *  dd.x.^2 .+ rand(twister, Normal(),obs) .+ 1
    iso_increasing = UnivariateFunctions.isotonic_regression(dd, :x, :y; increasing = true)
    plot(iso_increasing, dd)
    @test iso_increasing isa Piecewise_Function
    @test iso_increasing(1.9) > iso_increasing(0.1)
    iso_decreasing = UnivariateFunctions.isotonic_regression(dd, :x, :y; increasing = false)
    plot(iso_decreasing, dd)
    @test iso_decreasing isa Piecewise_Function
    @test iso_decreasing(1.9) <= iso_decreasing(0.1)

    # increasing data.
    dd = DataFrame(x = rand(twister, obs) .* 2)
    dd[!, :y] = 1 *  dd.x.^2 .+ rand(twister, Normal(),obs) .+ 1
    mono_increasing = UnivariateFunctions.monotonic_regression(dd, :x, :y; nbins=10, equally_spaced_bins=true, increasing=true)
    plot(mono_increasing, dd)
    @test mono_increasing isa Piecewise_Function
    @test mono_increasing(1.9) > mono_increasing(0.1)
    mono_decreasing = UnivariateFunctions.monotonic_regression(dd, :x, :y; nbins=5, equally_spaced_bins=false, increasing = false)
    plot(mono_decreasing, dd)
    @test mono_decreasing isa Piecewise_Function
    @test mono_decreasing(1.9) <= mono_decreasing(0.1)
end
