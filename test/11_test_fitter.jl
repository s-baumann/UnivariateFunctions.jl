using Test

@testset "Simplify Tests" begin
    using UnivariateFunctions

    tol = 10 * eps()

    # Build a piecewise linear function from known data
    x = collect(range(0.0, 10.0, length=11))
    y = sin.(x)
    pw = create_linear_interpolation(x, y)

    # simplify with same number of points and linear method should recover the same values
    pw_simplified = simplify(pw, 11, 0.0, 10.0, :linear)
    for xi in x
        @test abs(evaluate(pw, xi) - evaluate(pw_simplified, xi)) < tol
    end

    # simplify with :quadratic method
    pw_quad = simplify(pw, 11, 0.0, 10.0, :quadratic)
    for xi in x
        @test abs(evaluate(pw, xi) - evaluate(pw_quad, xi)) < 0.01
    end

    # simplify with :constant_right
    pw_cr = simplify(pw, 11, 0.0, 10.0, :constant_right)
    for xi in x
        @test abs(evaluate(pw, xi) - evaluate(pw_cr, xi)) < 1.0  # constants are coarser
    end

    # simplify with :constant_left
    pw_cl = simplify(pw, 11, 0.0, 10.0, :constant_left)
    @test pw_cl isa Piecewise_Function

    # Error on n_points < 2
    @test_throws ErrorException simplify(pw, 1, 0.0, 10.0)

    # Error on unknown method
    @test_throws ErrorException simplify(pw, 10, 0.0, 10.0, :unknown)

    # simplify with fewer points produces a coarser but valid function
    pw_coarse = simplify(pw, 5, 0.0, 10.0, :linear)
    @test pw_coarse isa Piecewise_Function
    @test abs(evaluate(pw_coarse, 0.0) - evaluate(pw, 0.0)) < tol
    @test abs(evaluate(pw_coarse, 10.0) - evaluate(pw, 10.0)) < tol
end

@testset "UnivariateFitter Constructor Tests" begin
    using UnivariateFunctions

    tol = 10 * eps()

    # Keyword constructor with defaults
    fitter = UnivariateFitter(:increasing)
    @test fitter.method == :increasing
    @test fitter.times_through == 0
    @test fitter.nbins == 10
    @test fitter.equally_spaced_bins == true
    @test fitter.weight_on_new == 1.0
    @test fitter.simplification_frequency == 0
    @test fitter.left_ == -Inf
    @test fitter.right_ == Inf
    # Initial function is zero
    @test abs(fitter(0.5)) < tol

    # Keyword constructor with custom values
    fitter2 = UnivariateFitter(:supersmoother; nbins=20, weight_on_new=0.5,
                                simplification_frequency=3, left=0.0, right=10.0)
    @test fitter2.method == :supersmoother
    @test fitter2.nbins == 20
    @test fitter2.weight_on_new == 0.5
    @test fitter2.simplification_frequency == 3
    @test fitter2.left_ == 0.0
    @test fitter2.right_ == 10.0
end

@testset "UnivariateFitter fit! Tests" begin
    using UnivariateFunctions
    using Random

    Random.seed!(42)
    tol = 0.5  # tolerance for regression approximation

    n = 200
    x = sort(rand(n) .* 10.0)
    y_increasing = 2.0 .* x .+ randn(n) .* 0.1

    # Test :increasing method
    fitter = UnivariateFitter(:increasing; nbins=15, left=0.0, right=10.0)
    UnivariateFunctions.fit!(fitter, x, y_increasing)
    @test fitter.times_through == 1
    @test fitter.fun isa UnivariateFunction
    # Check it's roughly increasing: evaluate at a few points
    vals = [evaluate(fitter, xi) for xi in 0.0:2.0:10.0]
    for i in 2:length(vals)
        @test vals[i] >= vals[i-1] - tol
    end

    # Test :decreasing method
    y_decreasing = -1.5 .* x .+ randn(n) .* 0.1
    fitter_dec = UnivariateFitter(:decreasing; nbins=15, left=0.0, right=10.0)
    UnivariateFunctions.fit!(fitter_dec, x, y_decreasing)
    @test fitter_dec.times_through == 1
    vals_dec = [evaluate(fitter_dec, xi) for xi in 0.0:2.0:10.0]
    for i in 2:length(vals_dec)
        @test vals_dec[i] <= vals_dec[i-1] + tol
    end

    # Test :supersmoother method
    y_smooth = sin.(x) .+ randn(n) .* 0.05
    fitter_ss = UnivariateFitter(:supersmoother; left=0.0, right=10.0)
    UnivariateFunctions.fit!(fitter_ss, x, y_smooth)
    @test fitter_ss.times_through == 1
    @test fitter_ss.fun isa UnivariateFunction
    # Check at a known point
    @test abs(evaluate(fitter_ss, 0.0) - sin(0.0)) < 1.0

    # Test :convex method (single minimum)
    y_convex = (x .- 5.0).^2 .+ randn(n) .* 0.1
    fitter_cvx = UnivariateFitter(:convex; nbins=15, left=0.0, right=10.0)
    UnivariateFunctions.fit!(fitter_cvx, x, y_convex)
    @test fitter_cvx.times_through == 1

    # Test :concave method (single maximum)
    y_concave = -(x .- 5.0).^2 .+ randn(n) .* 0.1
    fitter_ccv = UnivariateFitter(:concave; nbins=15, left=0.0, right=10.0)
    UnivariateFunctions.fit!(fitter_ccv, x, y_concave)
    @test fitter_ccv.times_through == 1

    # Test :quasiconvex method
    fitter_qcvx = UnivariateFitter(:quasiconvex; nbins=15, left=0.0, right=10.0)
    UnivariateFunctions.fit!(fitter_qcvx, x, y_convex)
    @test fitter_qcvx.times_through == 1

    # Test :quasiconcave method
    fitter_qccv = UnivariateFitter(:quasiconcave; nbins=15, left=0.0, right=10.0)
    UnivariateFunctions.fit!(fitter_qccv, x, y_concave)
    @test fitter_qccv.times_through == 1

    # Test unknown method
    fitter_bad = UnivariateFitter(:nonexistent; left=0.0, right=10.0)
    @test_throws ErrorException UnivariateFunctions.fit!(fitter_bad, x, y_increasing)
end

@testset "UnivariateFitter Iterative Fitting Tests" begin
    using UnivariateFunctions
    using Random

    Random.seed!(123)
    n = 200
    x = sort(rand(n) .* 10.0)
    y = 2.0 .* x .+ randn(n) .* 0.1

    # Test blending: weight_on_new < 1.0
    fitter = UnivariateFitter(:increasing; nbins=15, weight_on_new=0.7, left=0.0, right=10.0)
    UnivariateFunctions.fit!(fitter, x, y)
    @test fitter.times_through == 1
    val_after_first = evaluate(fitter, 5.0)

    # Second fit with same data - blends with previous
    UnivariateFunctions.fit!(fitter, x, y)
    @test fitter.times_through == 2
    val_after_second = evaluate(fitter, 5.0)
    # With weight_on_new=0.7, the result should still be close
    @test abs(val_after_second - val_after_first) < 2.0

    # Test simplification_frequency
    fitter_simp = UnivariateFitter(:increasing; nbins=15, weight_on_new=1.0,
                                    simplification_frequency=2, left=0.0, right=10.0)
    UnivariateFunctions.fit!(fitter_simp, x, y)
    @test fitter_simp.times_through == 1
    # times_through=1 after first fit, so second fit (times_through becomes 1 before mod check...
    # actually times_through is incremented at end)
    # After first fit: times_through=1
    # The simplification fires when times_through % simplification_frequency == 0
    # times_through is checked before increment... no, it increments at the end.
    # Actually: the code checks (fitter.times_through % fitter.simplification_frequency == 0)
    # BEFORE incrementing. So at first call, times_through=0, 0%2==0, so simplification fires.
    # Let's verify the function is still valid after simplification
    @test fitter_simp.fun isa UnivariateFunction
    @test abs(evaluate(fitter_simp, 5.0)) < 20.0  # sanity check

    UnivariateFunctions.fit!(fitter_simp, x, y)
    @test fitter_simp.times_through == 2
    @test fitter_simp.fun isa UnivariateFunction
end

@testset "UnivariateFitter Callable Tests" begin
    using UnivariateFunctions

    tol = 10 * eps()

    # Test that fitter is callable as a functor and via evaluate
    fitter = UnivariateFitter(:increasing; left=0.0, right=10.0)
    @test abs(fitter(5.0) - evaluate(fitter, 5.0)) < tol

    # After fitting, callable should still work
    x = collect(range(0.0, 10.0, length=50))
    y = 2.0 .* x
    UnivariateFunctions.fit!(fitter, x, y)
    @test abs(fitter(5.0) - evaluate(fitter, 5.0)) < tol
end
