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

    plot([iso_increasing, mono_increasing], 0.0, 2.0)

    # ========== Weighted isotonic regression ==========
    twister2 = MersenneTwister(10)
    obs2 = 200
    x_w = rand(twister2, obs2) .* 2
    y_w = x_w.^2 .+ rand(twister2, Normal(), obs2) .* 0.5

    # Equal weights should match unweighted
    iso_unweighted = isotonic_regression(x_w, y_w; increasing=true)
    iso_equal_w = isotonic_regression(x_w, y_w; increasing=true, weights=ones(obs2))
    @test abs(iso_unweighted(1.0) - iso_equal_w(1.0)) < 1e-8
    @test abs(iso_unweighted(0.5) - iso_equal_w(0.5)) < 1e-8

    # Non-uniform weights should produce a valid result
    w = rand(twister2, obs2) .+ 0.1
    iso_weighted = isotonic_regression(x_w, y_w; increasing=true, weights=w)
    @test iso_weighted isa Piecewise_Function
    @test iso_weighted(1.9) >= iso_weighted(0.1)

    # Decreasing with weights
    iso_weighted_dec = isotonic_regression(x_w, y_w; increasing=false, weights=w)
    @test iso_weighted_dec isa Piecewise_Function
    @test iso_weighted_dec(1.9) <= iso_weighted_dec(0.1)

    # ========== Weighted monotonic regression ==========

    # Equal weights should match unweighted
    mono_unweighted = monotonic_regression(x_w, y_w; increasing=true, nbins=10)
    mono_equal_w = monotonic_regression(x_w, y_w; increasing=true, nbins=10, weights=ones(obs2))
    @test abs(mono_unweighted(1.0) - mono_equal_w(1.0)) < 1e-6
    @test abs(mono_unweighted(0.5) - mono_equal_w(0.5)) < 1e-6

    # Non-uniform weights should produce a valid result
    mono_weighted = monotonic_regression(x_w, y_w; increasing=true, weights=w)
    @test mono_weighted isa Piecewise_Function
    @test mono_weighted(1.9) >= mono_weighted(0.1)

    # Decreasing with weights
    mono_weighted_dec = monotonic_regression(x_w, y_w; increasing=false, weights=w)
    @test mono_weighted_dec isa Piecewise_Function
    @test mono_weighted_dec(1.9) <= mono_weighted_dec(0.1)

    # High weight on one region should pull the fit toward that region
    # Create data where low-x points suggest y=0 and high-x points suggest y=10
    x_pull = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    y_pull = [0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
    # Heavy weight on high values should raise the fit at the midpoint
    w_high = [0.1, 0.1, 0.1, 0.1, 0.1, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
    w_low  = [10.0, 10.0, 10.0, 10.0, 10.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    mono_w_high = monotonic_regression(x_pull, y_pull; increasing=true, weights=w_high, nbins=5)
    mono_w_low  = monotonic_regression(x_pull, y_pull; increasing=true, weights=w_low, nbins=5)
    # With high weights on the upper region, the midpoint fit should be higher
    @test mono_w_high(2.5) >= mono_w_low(2.5) - 0.1

    # DataFrame interface with weights
    dd_w = DataFrame(x = x_w, y = y_w)
    iso_df_w = isotonic_regression(dd_w, :x, :y; increasing=true, weights=w)
    @test iso_df_w isa Piecewise_Function
    mono_df_w = monotonic_regression(dd_w, :x, :y; increasing=true, weights=w)
    @test mono_df_w isa Piecewise_Function

end
