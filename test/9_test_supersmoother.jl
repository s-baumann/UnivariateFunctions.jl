using Test

@testset "SuperSmoother Tests" begin
    using UnivariateFunctions
    using Dates, Random, DataFrames, Distributions

    tol = 10*eps()

    obs = 200
    twister = MersenneTwister(42)

    # Test with sinusoidal data
    x = collect(range(0, 2π, length=obs))
    y = sin.(x) .+ rand(twister, Normal(0, 0.2), obs)

    ss_fit = supersmoother(x, y)
    @test ss_fit isa Piecewise_Function

    # Check that smoothed values are reasonable (close to true sin curve)
    mid_point = π
    @test abs(ss_fit(mid_point) - sin(mid_point)) < 0.5  # allowing for noise

    # Test with DataFrame interface
    dd = DataFrame(x = x, y = y)
    ss_fit_df = supersmoother(dd, :x, :y)
    @test ss_fit_df isa Piecewise_Function

    # Test with different spans
    ss_fit_smooth = supersmoother(x, y; spans=[0.1, 0.3, 0.6])
    @test ss_fit_smooth isa Piecewise_Function

    # Test with bass enhancement (should give smoother result)
    ss_fit_bass = supersmoother(x, y; bass=5.0)
    @test ss_fit_bass isa Piecewise_Function

    # Test that derivative works
    deriv = derivative(ss_fit)
    @test deriv isa Piecewise_Function

    # Test that integral works
    integ = indefinite_integral(ss_fit)
    @test integ isa Piecewise_Function

    # Test with monotonic data - should still work
    x_mono = collect(range(0, 5, length=100))
    y_mono = 2.0 .* x_mono .+ rand(twister, Normal(0, 0.5), 100)
    ss_mono = supersmoother(x_mono, y_mono)
    @test ss_mono isa Piecewise_Function
    @test ss_mono(4.0) > ss_mono(1.0)  # should be increasing

    # Test definite integral
    area = evaluate_integral(ss_fit, 0.0, 2π)
    @test area isa Real
    # For sin(x) from 0 to 2π, integral should be close to 0
    @test abs(area) < 2.0  # allowing for noise

    # Test evaluation at multiple points (broadcasting)
    x_eval = collect(0.0:0.5:2π)
    y_eval = ss_fit.(x_eval)
    @test length(y_eval) == length(x_eval)
    @test all(y -> y isa Real, y_eval)

    # Test arithmetic operations on fitted function
    scaled = 2.0 * ss_fit
    @test scaled isa Piecewise_Function
    @test scaled(1.0) ≈ 2.0 * ss_fit(1.0)

    shifted = ss_fit + PE_Function(1.0)
    @test shifted isa Piecewise_Function
    @test shifted(1.0) ≈ ss_fit(1.0) + 1.0

    # Test with smaller dataset (edge case)
    x_small = collect(range(0, 1, length=20))
    y_small = x_small.^2 .+ rand(twister, Normal(0, 0.1), 20)
    ss_small = supersmoother(x_small, y_small)
    @test ss_small isa Piecewise_Function

    println("SuperSmoother tests passed.")
end
