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

# ====================================================================
# UnivariateAdjustedFitter tests
# ====================================================================

@testset "UnivariateAdjustedFitter Constructor Tests" begin
    using UnivariateFunctions

    tol = 10 * eps()

    # Default constructor
    fitter = UnivariateAdjustedFitter(:increasing)
    @test fitter.method == :increasing
    @test fitter.times_through == 0
    @test fitter.nbins == 10
    @test fitter.equally_spaced_bins == true
    @test fitter.weight_on_new == 1.0
    @test fitter.simplification_frequency == 0
    @test fitter.left_ == -Inf
    @test fitter.right_ == Inf
    @test fitter.adjust_for_groups == true
    @test fitter.fit_intercept == true
    @test fitter.coefficient_bounds == ((-1.0, 1.0), (0.1, 2.5))
    @test isempty(fitter.coefficients)
    # Initial function is zero
    @test abs(fitter.fun(0.5)) < tol

    # Custom constructor
    fitter2 = UnivariateAdjustedFitter(:decreasing; nbins=20, weight_on_new=0.5,
                                        simplification_frequency=3, left=0.0, right=10.0,
                                        adjust_for_groups=false, fit_intercept=false,
                                        coefficient_bounds=((-5.0, 5.0), (0.5, 3.0)))
    @test fitter2.method == :decreasing
    @test fitter2.nbins == 20
    @test fitter2.weight_on_new == 0.5
    @test fitter2.simplification_frequency == 3
    @test fitter2.left_ == 0.0
    @test fitter2.right_ == 10.0
    @test fitter2.adjust_for_groups == false
    @test fitter2.fit_intercept == false
    @test fitter2.coefficient_bounds == ((-5.0, 5.0), (0.5, 3.0))

    # Type stability of coefficient_bounds
    @test fitter.coefficient_bounds isa Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}
end

@testset "UnivariateAdjustedFitter Single Group Tests" begin
    using UnivariateFunctions
    using Random

    Random.seed!(99)
    n = 200
    x = sort(rand(n) .* 10.0)
    y = 2.0 .* x .+ 1.0 .+ randn(n) .* 0.1
    groups = fill(:A, n)

    fitter = UnivariateAdjustedFitter(:increasing; nbins=15, left=0.0, right=10.0,
                                       coefficient_bounds=((-5.0, 5.0), (0.1, 5.0)))
    UnivariateFunctions.fit!(fitter, x, y, groups)

    # Basic state checks
    @test fitter.times_through == 1
    @test fitter.fun isa UnivariateFunction
    @test haskey(fitter.coefficients, :A)

    # Coefficients should be finite and within bounds
    a_A, b_A = fitter.coefficients[:A]
    @test isfinite(a_A) && isfinite(b_A)
    @test a_A >= -5.0 && a_A <= 5.0
    @test b_A >= 0.1 && b_A <= 5.0

    # Evaluate via group should return a finite number
    val = evaluate(fitter, 5.0, :A)
    @test isfinite(val)

    # Functor call should match evaluate
    @test abs(fitter(5.0, :A) - evaluate(fitter, 5.0, :A)) < 1e-12

    # The fitted values should roughly track the data
    @test abs(fitter(0.0, :A) - 1.0) < 3.0
    @test abs(fitter(10.0, :A) - 21.0) < 3.0
end

@testset "UnivariateAdjustedFitter Multi-Group Tests" begin
    using UnivariateFunctions
    using Random

    Random.seed!(77)
    n = 150

    # Shared shape: f(x) = x (increasing)
    # Group A: y = 1.0 + 2.0 * x + noise
    # Group B: y = -0.5 + 1.5 * x + noise
    x = sort(rand(n) .* 10.0)
    groups = [isodd(i) ? :A : :B for i in 1:n]

    y = Vector{Float64}(undef, n)
    for i in 1:n
        if groups[i] == :A
            y[i] = 1.0 + 2.0 * x[i] + randn() * 0.1
        else
            y[i] = -0.5 + 1.5 * x[i] + randn() * 0.1
        end
    end

    fitter = UnivariateAdjustedFitter(:increasing; nbins=15, left=0.0, right=10.0,
                                       coefficient_bounds=((-5.0, 5.0), (0.1, 5.0)))
    UnivariateFunctions.fit!(fitter, x, y, groups)

    @test fitter.times_through == 1
    @test haskey(fitter.coefficients, :A)
    @test haskey(fitter.coefficients, :B)

    # Both groups should produce reasonable predictions
    @test isfinite(fitter(5.0, :A))
    @test isfinite(fitter(5.0, :B))

    # Group A should predict higher than Group B at the same x (since a_A > a_B and b_A > b_B)
    @test fitter(5.0, :A) > fitter(5.0, :B)

    # Coefficients should differ between groups
    a_A, b_A = fitter.coefficients[:A]
    a_B, b_B = fitter.coefficients[:B]
    @test !(a_A ≈ a_B && b_A ≈ b_B)
end

@testset "UnivariateAdjustedFitter New Group Onboarding" begin
    using UnivariateFunctions
    using Random

    Random.seed!(55)
    n = 100
    x = sort(rand(n) .* 10.0)
    y = 2.0 .* x .+ randn(n) .* 0.1

    fitter = UnivariateAdjustedFitter(:increasing; nbins=15, left=0.0, right=10.0,
                                       coefficient_bounds=((-5.0, 5.0), (0.1, 5.0)))

    # First fit with group A
    UnivariateFunctions.fit!(fitter, x, y, fill(:A, n))
    @test haskey(fitter.coefficients, :A)
    @test !haskey(fitter.coefficients, :B)

    # Second fit introduces group B
    y2 = 3.0 .* x .+ randn(n) .* 0.1
    UnivariateFunctions.fit!(fitter, x, y2, fill(:B, n))
    @test haskey(fitter.coefficients, :B)
    @test fitter.times_through == 2

    # Group A coefficients should still exist (unchanged since not in second batch)
    @test haskey(fitter.coefficients, :A)
end

@testset "UnivariateAdjustedFitter Coefficient Bounds Clamping" begin
    using UnivariateFunctions
    using Random

    Random.seed!(42)
    n = 100
    x = sort(rand(n) .* 10.0)
    # Very large intercept and slope — should get clamped
    y = 100.0 .+ 50.0 .* x .+ randn(n) .* 0.1
    groups = fill(:A, n)

    # Tight bounds
    fitter = UnivariateAdjustedFitter(:increasing; nbins=15, left=0.0, right=10.0,
                                       coefficient_bounds=((-2.0, 2.0), (0.5, 3.0)))
    UnivariateFunctions.fit!(fitter, x, y, groups)

    a_A, b_A = fitter.coefficients[:A]
    @test a_A >= -2.0
    @test a_A <= 2.0
    @test b_A >= 0.5
    @test b_A <= 3.0
end

@testset "UnivariateAdjustedFitter Iterative Blending" begin
    using UnivariateFunctions
    using Random

    Random.seed!(33)
    n = 200
    x = sort(rand(n) .* 10.0)
    y = 2.0 .* x .+ randn(n) .* 0.1
    groups = fill(:A, n)

    fitter = UnivariateAdjustedFitter(:increasing; nbins=15, weight_on_new=0.7,
                                       left=0.0, right=10.0,
                                       coefficient_bounds=((-5.0, 5.0), (0.1, 5.0)))

    UnivariateFunctions.fit!(fitter, x, y, groups)
    val_after_first = fitter(5.0, :A)
    coeffs_after_first = fitter.coefficients[:A]

    UnivariateFunctions.fit!(fitter, x, y, groups)
    val_after_second = fitter(5.0, :A)
    coeffs_after_second = fitter.coefficients[:A]

    @test fitter.times_through == 2
    # Values should be close since data is the same
    @test abs(val_after_second - val_after_first) < 5.0
    # Coefficients should have been blended
    @test coeffs_after_second != coeffs_after_first || val_after_second != val_after_first
end

@testset "UnivariateAdjustedFitter Simplification" begin
    using UnivariateFunctions
    using Random

    Random.seed!(11)
    n = 200
    x = sort(rand(n) .* 10.0)
    y = 2.0 .* x .+ randn(n) .* 0.1
    groups = fill(:A, n)

    fitter = UnivariateAdjustedFitter(:increasing; nbins=15, weight_on_new=1.0,
                                       simplification_frequency=2, left=0.0, right=10.0,
                                       min_bins_for_simplification=25,
                                       coefficient_bounds=((-5.0, 5.0), (0.1, 5.0)))

    # First fit: times_through=0 before check, 0%2==0 but guard requires times_through > 0
    UnivariateFunctions.fit!(fitter, x, y, groups)
    @test fitter.times_through == 1
    @test fitter.fun isa UnivariateFunction

    # Second fit: times_through=1 before check, 1%2!=0 so no simplification
    UnivariateFunctions.fit!(fitter, x, y, groups)
    @test fitter.times_through == 2

    # Third fit: times_through=2 before check, 2%2==0 and times_through>0 — simplification fires
    UnivariateFunctions.fit!(fitter, x, y, groups)
    @test fitter.times_through == 3
    @test fitter.fun isa UnivariateFunction
    # Function should still produce reasonable results after simplification
    @test isfinite(fitter(5.0, :A))
end

@testset "UnivariateAdjustedFitter adjust_for_groups=false" begin
    using UnivariateFunctions
    using Random

    Random.seed!(22)
    n = 100
    x = sort(rand(n) .* 10.0)
    y = 2.0 .* x .+ 1.0 .+ randn(n) .* 0.1
    groups = fill(:A, n)

    # With adjust_for_groups=false, y is passed to shape fitter directly
    fitter = UnivariateAdjustedFitter(:increasing; nbins=15, left=0.0, right=10.0,
                                       adjust_for_groups=false,
                                       coefficient_bounds=((-5.0, 5.0), (0.1, 5.0)))
    UnivariateFunctions.fit!(fitter, x, y, groups)

    @test fitter.times_through == 1
    @test isfinite(fitter(5.0, :A))

    # Coefficients should still be estimated even when adjust_for_groups=false
    a_A, b_A = fitter.coefficients[:A]
    @test isfinite(a_A) && isfinite(b_A)
end

@testset "UnivariateAdjustedFitter All Methods" begin
    using UnivariateFunctions
    using Random

    Random.seed!(88)
    n = 100
    x = sort(rand(n) .* 10.0)
    groups = fill(:A, n)

    for method in [:increasing, :decreasing, :convex, :concave, :quasiconvex, :quasiconcave, :supersmoother]
        y = if method in [:increasing, :quasiconvex]
            2.0 .* x .+ randn(n) .* 0.1
        elseif method in [:decreasing, :quasiconcave]
            -2.0 .* x .+ randn(n) .* 0.1
        elseif method == :convex
            (x .- 5.0).^2 .+ randn(n) .* 0.1
        elseif method == :concave
            -(x .- 5.0).^2 .+ randn(n) .* 0.1
        else  # supersmoother
            sin.(x) .+ randn(n) .* 0.05
        end

        fitter = UnivariateAdjustedFitter(method; nbins=15, left=0.0, right=10.0,
                                           coefficient_bounds=((-50.0, 50.0), (0.1, 10.0)))
        UnivariateFunctions.fit!(fitter, x, y, groups)
        @test fitter.times_through == 1
        @test fitter.fun isa UnivariateFunction
        @test isfinite(fitter(5.0, :A))
    end

    # Unknown method should error
    fitter_bad = UnivariateAdjustedFitter(:nonexistent; left=0.0, right=10.0)
    @test_throws ErrorException UnivariateFunctions.fit!(fitter_bad, x, 2.0 .* x, groups)
end

@testset "UnivariateAdjustedFitter fit_intercept=false Single Group" begin
    using UnivariateFunctions
    using Random

    Random.seed!(44)
    n = 200
    x = sort(rand(n) .* 10.0)
    # y = 3.0 * f(x) with no intercept, where f is roughly linear increasing
    y = 3.0 .* (2.0 .* x) .+ randn(n) .* 0.1
    groups = fill(:A, n)

    fitter = UnivariateAdjustedFitter(:increasing; nbins=15, left=0.0, right=10.0,
                                       fit_intercept=false,
                                       coefficient_bounds=((-5.0, 5.0), (0.1, 10.0)))
    UnivariateFunctions.fit!(fitter, x, y, groups)

    @test fitter.times_through == 1
    @test fitter.fit_intercept == false

    # Intercept should be forced to zero (within bounds and blending)
    a_A, b_A = fitter.coefficients[:A]
    @test a_A == 0.0

    # b should be finite and positive
    @test isfinite(b_A)
    @test b_A > 0.0

    # Evaluate should work
    @test isfinite(fitter(5.0, :A))
end

@testset "UnivariateAdjustedFitter fit_intercept=false Multi-Group" begin
    using UnivariateFunctions
    using Random

    Random.seed!(66)
    n = 200

    # Shared shape: f(x) ≈ x (increasing)
    # Group A: y = 2.0 * f(x) + noise  (scale=2)
    # Group B: y = 0.5 * f(x) + noise  (scale=0.5)
    x = sort(rand(n) .* 10.0)
    groups = [isodd(i) ? :A : :B for i in 1:n]
    y = Vector{Float64}(undef, n)
    for i in 1:n
        if groups[i] == :A
            y[i] = 2.0 * x[i] + randn() * 0.1
        else
            y[i] = 0.5 * x[i] + randn() * 0.1
        end
    end

    fitter = UnivariateAdjustedFitter(:increasing; nbins=15, left=0.0, right=10.0,
                                       fit_intercept=false,
                                       coefficient_bounds=((-5.0, 5.0), (0.1, 10.0)))
    UnivariateFunctions.fit!(fitter, x, y, groups)

    @test haskey(fitter.coefficients, :A)
    @test haskey(fitter.coefficients, :B)

    # Both intercepts should be zero
    a_A, b_A = fitter.coefficients[:A]
    a_B, b_B = fitter.coefficients[:B]
    @test a_A == 0.0
    @test a_B == 0.0

    # Group A should have a larger scale than Group B
    @test b_A > b_B

    # Both should produce finite predictions
    @test isfinite(fitter(5.0, :A))
    @test isfinite(fitter(5.0, :B))

    # Group A predictions should be larger at positive x
    @test fitter(5.0, :A) > fitter(5.0, :B)
end

@testset "UnivariateAdjustedFitter fit_intercept=false Iterative" begin
    using UnivariateFunctions
    using Random

    Random.seed!(55)
    n = 200
    x = sort(rand(n) .* 10.0)
    y = 2.0 .* x .+ randn(n) .* 0.1
    groups = fill(:A, n)

    fitter = UnivariateAdjustedFitter(:increasing; nbins=15, weight_on_new=0.8,
                                       left=0.0, right=10.0, fit_intercept=false,
                                       coefficient_bounds=((-5.0, 5.0), (0.1, 10.0)))

    # Multiple iterations
    for _ in 1:3
        UnivariateFunctions.fit!(fitter, x, y, groups)
    end
    @test fitter.times_through == 3

    # Intercept should remain zero after blending
    a_A, _ = fitter.coefficients[:A]
    @test a_A == 0.0

    @test isfinite(fitter(5.0, :A))
end

@testset "UnivariateAdjustedFitter fit_intercept=false vs true Comparison" begin
    using UnivariateFunctions
    using Random

    Random.seed!(77)
    n = 200
    x = sort(rand(n) .* 10.0)
    # Data with a clear intercept: y = 5.0 + 2.0*x
    y = 5.0 .+ 2.0 .* x .+ randn(n) .* 0.1
    groups = fill(:A, n)

    # With intercept
    fitter_with = UnivariateAdjustedFitter(:increasing; nbins=15, left=0.0, right=10.0,
                                            fit_intercept=true,
                                            coefficient_bounds=((-10.0, 10.0), (0.1, 10.0)))
    UnivariateFunctions.fit!(fitter_with, x, y, groups)

    # Without intercept
    fitter_without = UnivariateAdjustedFitter(:increasing; nbins=15, left=0.0, right=10.0,
                                               fit_intercept=false,
                                               coefficient_bounds=((-10.0, 10.0), (0.1, 10.0)))
    UnivariateFunctions.fit!(fitter_without, x, y, groups)

    # fit_intercept=true should have a nonzero intercept estimate
    a_with, _ = fitter_with.coefficients[:A]
    a_without, _ = fitter_without.coefficients[:A]
    @test a_without == 0.0
    # The with-intercept version should have a better fit at x=0 (closer to 5.0)
    @test abs(fitter_with(0.0, :A) - 5.0) < abs(fitter_without(0.0, :A) - 5.0)
end

@testset "UnivariateAdjustedFitter fit_intercept=false Bounds Clamping" begin
    using UnivariateFunctions
    using Random

    Random.seed!(88)
    n = 100
    x = sort(rand(n) .* 10.0)
    # Very large scale — b should get clamped
    y = 100.0 .* x .+ randn(n) .* 0.1
    groups = fill(:A, n)

    fitter = UnivariateAdjustedFitter(:increasing; nbins=15, left=0.0, right=10.0,
                                       fit_intercept=false,
                                       coefficient_bounds=((-2.0, 2.0), (0.5, 3.0)))
    UnivariateFunctions.fit!(fitter, x, y, groups)

    a_A, b_A = fitter.coefficients[:A]
    @test a_A == 0.0
    @test b_A >= 0.5
    @test b_A <= 3.0
end

@testset "UnivariateFitter Weighted fit! Tests" begin
    using UnivariateFunctions
    using Random, Distributions

    n = 200
    twister = MersenneTwister(88)
    x = sort(rand(twister, n) .* 5.0)
    y = 2.0 .* x .+ 1.0 .+ rand(twister, Normal(0, 0.3), n)
    w = rand(twister, n) .+ 0.1

    # UnivariateFitter with weights
    fitter_w = UnivariateFitter(:increasing; nbins=10)
    UnivariateFunctions.fit!(fitter_w, x, y; weights=w)
    @test fitter_w.fun isa Piecewise_Function
    @test fitter_w(4.0) > fitter_w(1.0)

    # Equal weights should match unweighted
    fitter_uw = UnivariateFitter(:increasing; nbins=10)
    UnivariateFunctions.fit!(fitter_uw, x, y)
    fitter_ew = UnivariateFitter(:increasing; nbins=10)
    UnivariateFunctions.fit!(fitter_ew, x, y; weights=ones(n))
    @test abs(fitter_uw(2.5) - fitter_ew(2.5)) < 1e-6

    # Iterative fitting with weights
    fitter_iter = UnivariateFitter(:increasing; nbins=10, weight_on_new=0.7)
    UnivariateFunctions.fit!(fitter_iter, x, y; weights=w)
    UnivariateFunctions.fit!(fitter_iter, x, y; weights=w)
    @test fitter_iter.times_through == 2
    @test fitter_iter(3.0) isa Real

    # Supersmoother with weights
    fitter_ss = UnivariateFitter(:supersmoother)
    UnivariateFunctions.fit!(fitter_ss, x, y; weights=w)
    @test fitter_ss.fun isa Piecewise_Function

    # Quasiconcave with weights
    y_peak = -(x .- 2.5).^2 .+ 6.0 .+ rand(twister, Normal(0, 0.3), n)
    fitter_qc = UnivariateFitter(:quasiconcave; nbins=10)
    UnivariateFunctions.fit!(fitter_qc, x, y_peak; weights=w)
    @test fitter_qc.fun isa Piecewise_Function
end

@testset "UnivariateAdjustedFitter Weighted fit! Tests" begin
    using UnivariateFunctions
    using Random, Distributions

    n = 200
    twister = MersenneTwister(99)
    x = sort(rand(twister, n) .* 5.0)
    groups = repeat([:A, :B], n ÷ 2)
    y = [g == :A ? 1.0 + 1.5 * xi : 0.5 + 0.8 * xi for (xi, g) in zip(x, groups)]
    y .+= rand(twister, Normal(0, 0.2), n)
    w = rand(twister, n) .+ 0.1

    # Basic weighted fit
    fitter_w = UnivariateAdjustedFitter(:increasing; nbins=10, left=0.0, right=5.0,
                                         coefficient_bounds=((-5.0, 5.0), (0.1, 5.0)))
    UnivariateFunctions.fit!(fitter_w, x, y, groups; weights=w)
    @test fitter_w.fun isa Piecewise_Function
    @test fitter_w.times_through == 1
    @test haskey(fitter_w.coefficients, :A)
    @test haskey(fitter_w.coefficients, :B)

    # Equal weights should match unweighted
    fitter_uw = UnivariateAdjustedFitter(:increasing; nbins=10, left=0.0, right=5.0,
                                          coefficient_bounds=((-5.0, 5.0), (0.1, 5.0)))
    UnivariateFunctions.fit!(fitter_uw, x, y, groups)
    fitter_ew = UnivariateAdjustedFitter(:increasing; nbins=10, left=0.0, right=5.0,
                                          coefficient_bounds=((-5.0, 5.0), (0.1, 5.0)))
    UnivariateFunctions.fit!(fitter_ew, x, y, groups; weights=ones(n))
    @test abs(fitter_uw(2.5, :A) - fitter_ew(2.5, :A)) < 1e-6
    @test abs(fitter_uw(2.5, :B) - fitter_ew(2.5, :B)) < 1e-6
    @test abs(fitter_uw.coefficients[:A][1] - fitter_ew.coefficients[:A][1]) < 1e-6
    @test abs(fitter_uw.coefficients[:A][2] - fitter_ew.coefficients[:A][2]) < 1e-6

    # Iterative weighted fitting
    fitter_iter = UnivariateAdjustedFitter(:increasing; nbins=10, weight_on_new=0.7,
                                            left=0.0, right=5.0,
                                            coefficient_bounds=((-5.0, 5.0), (0.1, 5.0)))
    UnivariateFunctions.fit!(fitter_iter, x, y, groups; weights=w)
    UnivariateFunctions.fit!(fitter_iter, x, y, groups; weights=w)
    @test fitter_iter.times_through == 2

    # fit_intercept=false with weights
    fitter_noint = UnivariateAdjustedFitter(:increasing; nbins=10, left=0.0, right=5.0,
                                             fit_intercept=false,
                                             coefficient_bounds=((-5.0, 5.0), (0.1, 5.0)))
    UnivariateFunctions.fit!(fitter_noint, x, y, groups; weights=w)
    @test fitter_noint.coefficients[:A][1] == 0.0
    @test fitter_noint.coefficients[:B][1] == 0.0
end
