using Test

@testset "Unimodal Regression Tests" begin
    using UnivariateFunctions
    using Dates, Random, DataFrames, Distributions

    tol = 10*eps()
    obs = 200
    twister = MersenneTwister(123)

    # ========== Test unimodal_regression ==========

    # Test quasiconcave (single peak)
    x = collect(range(-2, 2, length=obs))
    y_concave = -x.^2 .+ 4 .+ rand(twister, Normal(0, 0.3), obs)  # peak at x=0

    fit_qc = unimodal_regression(x, y_concave; quasi=true, convex=false)
    @test fit_qc isa Piecewise_Function
    # Peak should be near the middle
    @test fit_qc(0.0) > fit_qc(-1.5)
    @test fit_qc(0.0) > fit_qc(1.5)

    # Test quasiconvex (single trough)
    y_convex = x.^2 .+ rand(twister, Normal(0, 0.3), obs)  # trough at x=0

    fit_qv = unimodal_regression(x, y_convex; quasi=true, convex=true)
    @test fit_qv isa Piecewise_Function
    # Trough should be near the middle
    @test fit_qv(0.0) < fit_qv(-1.5)
    @test fit_qv(0.0) < fit_qv(1.5)

    # Test true concave (decreasing slopes)
    fit_concave = unimodal_regression(x, y_concave; quasi=false, convex=false)
    @test fit_concave isa Piecewise_Function

    # Test true convex (increasing slopes)
    fit_convex = unimodal_regression(x, y_convex; quasi=false, convex=true)
    @test fit_convex isa Piecewise_Function

    # Test DataFrame interface
    dd = DataFrame(x = x, y = y_concave)
    fit_df = unimodal_regression(dd, :x, :y; quasi=true, convex=false)
    @test fit_df isa Piecewise_Function

    # Test with different nbins
    fit_fine = unimodal_regression(x, y_concave; nbins=20)
    @test fit_fine isa Piecewise_Function

    fit_coarse = unimodal_regression(x, y_concave; nbins=5)
    @test fit_coarse isa Piecewise_Function

    # Test equally_spaced_bins option
    fit_quantile = unimodal_regression(x, y_concave; equally_spaced_bins=false)
    @test fit_quantile isa Piecewise_Function

    # ========== Test cv_monotonic_regression ==========

    # Clearly increasing data
    x_inc = collect(range(0, 5, length=100))
    y_inc = 2.0 .* x_inc .+ rand(twister, Normal(0, 0.5), 100)

    cv_mono_inc = cv_monotonic_regression(x_inc, y_inc; nfolds=5, seed=42)
    @test cv_mono_inc isa CVRegressionResult
    @test cv_mono_inc.selected_shape == :increasing
    @test cv_mono_inc.fitted isa Piecewise_Function
    @test haskey(cv_mono_inc.cv_errors, :increasing)
    @test haskey(cv_mono_inc.cv_errors, :decreasing)
    # Test callable
    @test cv_mono_inc(2.5) isa Real

    # Clearly decreasing data
    y_dec = -2.0 .* x_inc .+ 10 .+ rand(twister, Normal(0, 0.5), 100)

    cv_mono_dec = cv_monotonic_regression(x_inc, y_dec; nfolds=5, seed=42)
    @test cv_mono_dec isa CVRegressionResult
    @test cv_mono_dec.selected_shape == :decreasing

    # ========== Test cv_unimodal_regression ==========

    # Concave data (single peak)
    x_uni = collect(range(-3, 3, length=150))
    y_peak = -x_uni.^2 .+ 9 .+ rand(twister, Normal(0, 0.5), 150)

    cv_uni_peak = cv_unimodal_regression(x_uni, y_peak; nfolds=5, seed=42)
    @test cv_uni_peak isa CVRegressionResult
    @test cv_uni_peak.selected_shape in [:quasiconcave, :concave]
    @test cv_uni_peak.fitted isa Piecewise_Function

    # Convex data (single trough)
    y_trough = x_uni.^2 .+ rand(twister, Normal(0, 0.5), 150)

    cv_uni_trough = cv_unimodal_regression(x_uni, y_trough; nfolds=5, seed=42)
    @test cv_uni_trough isa CVRegressionResult
    @test cv_uni_trough.selected_shape in [:quasiconvex, :convex]

    # ========== Test cv_shape_regression ==========

    # Test with :all shapes
    cv_all = cv_shape_regression(x_uni, y_peak; shapes=:all, nfolds=5, seed=42)
    @test cv_all isa CVRegressionResult
    @test length(cv_all.cv_errors) == 6

    # Test with :monotonic shapes
    cv_mono = cv_shape_regression(x_inc, y_inc; shapes=:monotonic, nfolds=5, seed=42)
    @test cv_mono isa CVRegressionResult
    @test length(cv_mono.cv_errors) == 2
    @test cv_mono.selected_shape in [:increasing, :decreasing]

    # Test with :unimodal shapes
    cv_uni = cv_shape_regression(x_uni, y_peak; shapes=:unimodal, nfolds=5, seed=42)
    @test cv_uni isa CVRegressionResult
    @test length(cv_uni.cv_errors) == 4

    # Test with custom shape selection
    cv_custom = cv_shape_regression(x_uni, y_peak; shapes=[:increasing, :quasiconcave], nfolds=5, seed=42)
    @test cv_custom isa CVRegressionResult
    @test length(cv_custom.cv_errors) == 2

    # Test DataFrame interface
    dd_uni = DataFrame(x = x_uni, y = y_peak)
    cv_df = cv_shape_regression(dd_uni, :x, :y; shapes=:all, nfolds=5, seed=42)
    @test cv_df isa CVRegressionResult

    # Test that derivative works on fitted functions
    deriv = derivative(cv_all.fitted)
    @test deriv isa Piecewise_Function

    # ========== Test integration of fitted functions ==========

    # Test definite integral
    area = evaluate_integral(fit_qc, -2.0, 2.0)
    @test area isa Real

    # Test indefinite integral
    integ = indefinite_integral(fit_qc)
    @test integ isa Piecewise_Function

    # ========== Test DataFrame interfaces for CV functions ==========

    dd_inc = DataFrame(x = x_inc, y = y_inc)

    # cv_monotonic_regression DataFrame interface
    cv_mono_df = cv_monotonic_regression(dd_inc, :x, :y; nfolds=5, seed=42)
    @test cv_mono_df isa CVRegressionResult
    @test cv_mono_df.selected_shape == :increasing

    # cv_unimodal_regression DataFrame interface
    cv_uni_df = cv_unimodal_regression(dd, :x, :y; nfolds=5, seed=42)
    @test cv_uni_df isa CVRegressionResult

    # ========== Test arithmetic operations on fitted functions ==========

    scaled = 2.0 * fit_qc
    @test scaled isa Piecewise_Function
    @test scaled(0.0) ≈ 2.0 * fit_qc(0.0)

    shifted = fit_qc + PE_Function(1.0)
    @test shifted isa Piecewise_Function
    @test shifted(0.0) ≈ fit_qc(0.0) + 1.0

    # ========== Test that concave has decreasing slopes ==========

    # For true concave, slopes should decrease
    deriv_concave = derivative(fit_concave)
    # Check slopes at different points
    slope_left = deriv_concave(-1.5)
    slope_middle = deriv_concave(0.0)
    slope_right = deriv_concave(1.5)
    # Slopes should be decreasing (or at least not increasing)
    @test slope_left >= slope_middle - 0.5  # allow small tolerance
    @test slope_middle >= slope_right - 0.5

    # ========== Test evaluation at multiple points ==========

    x_eval = collect(-2.0:0.5:2.0)
    y_eval = fit_qc.(x_eval)
    @test length(y_eval) == length(x_eval)
    @test all(y -> y isa Real, y_eval)

    # ========== Test with smaller nbins (edge case) ==========

    fit_tiny = unimodal_regression(x, y_concave; nbins=3)
    @test fit_tiny isa Piecewise_Function

    # ========== Test CVRegressionResult callable interface ==========

    # Should be callable like a function
    @test cv_all(0.0) isa Real
    @test cv_all(0.0) == cv_all.fitted(0.0)

    # ========== Weighted unimodal_regression ==========

    twister_w = MersenneTwister(77)
    obs_w = 200
    x_w = collect(range(-2, 2, length=obs_w))
    y_w_concave = -x_w.^2 .+ 4 .+ rand(twister_w, Normal(0, 0.3), obs_w)

    # Equal weights should match unweighted
    fit_unweighted = unimodal_regression(x_w, y_w_concave; quasi=true, convex=false, nbins=10)
    fit_equal_w = unimodal_regression(x_w, y_w_concave; quasi=true, convex=false, nbins=10, weights=ones(obs_w))
    @test abs(fit_unweighted(0.0) - fit_equal_w(0.0)) < 1e-6
    @test abs(fit_unweighted(1.0) - fit_equal_w(1.0)) < 1e-6

    # Non-uniform weights
    w_uni = rand(twister_w, obs_w) .+ 0.1
    fit_weighted_qc = unimodal_regression(x_w, y_w_concave; quasi=true, convex=false, weights=w_uni)
    @test fit_weighted_qc isa Piecewise_Function
    @test fit_weighted_qc(0.0) > fit_weighted_qc(-1.5)
    @test fit_weighted_qc(0.0) > fit_weighted_qc(1.5)

    # Quasiconvex with weights
    y_w_convex = x_w.^2 .+ rand(twister_w, Normal(0, 0.3), obs_w)
    fit_weighted_qv = unimodal_regression(x_w, y_w_convex; quasi=true, convex=true, weights=w_uni)
    @test fit_weighted_qv isa Piecewise_Function
    @test fit_weighted_qv(0.0) < fit_weighted_qv(-1.5)
    @test fit_weighted_qv(0.0) < fit_weighted_qv(1.5)

    # True concave with weights
    fit_weighted_concave = unimodal_regression(x_w, y_w_concave; quasi=false, convex=false, weights=w_uni)
    @test fit_weighted_concave isa Piecewise_Function

    # True convex with weights
    fit_weighted_convex = unimodal_regression(x_w, y_w_convex; quasi=false, convex=true, weights=w_uni)
    @test fit_weighted_convex isa Piecewise_Function

    # DataFrame interface with weights
    dd_w = DataFrame(x = x_w, y = y_w_concave)
    fit_df_w = unimodal_regression(dd_w, :x, :y; quasi=true, convex=false, weights=w_uni)
    @test fit_df_w isa Piecewise_Function

    # ========== Weighted CV functions ==========

    twister_cv = MersenneTwister(55)
    obs_cv = 150
    x_cv = collect(range(0, 5, length=obs_cv))
    w_cv = rand(twister_cv, obs_cv) .+ 0.1

    # cv_monotonic_regression with weights
    y_cv_inc = 2.0 .* x_cv .+ rand(twister_cv, Normal(0, 0.5), obs_cv)
    cv_mono_w = cv_monotonic_regression(x_cv, y_cv_inc; nfolds=5, seed=42, weights=w_cv)
    @test cv_mono_w isa CVRegressionResult
    @test cv_mono_w.selected_shape == :increasing
    @test cv_mono_w.fitted isa Piecewise_Function
    @test cv_mono_w(2.5) isa Real

    # Equal weights should match unweighted
    cv_mono_uw = cv_monotonic_regression(x_cv, y_cv_inc; nfolds=5, seed=42)
    cv_mono_ew = cv_monotonic_regression(x_cv, y_cv_inc; nfolds=5, seed=42, weights=ones(obs_cv))
    @test cv_mono_uw.selected_shape == cv_mono_ew.selected_shape
    @test abs(cv_mono_uw.cv_errors[:increasing] - cv_mono_ew.cv_errors[:increasing]) < 1e-6

    # cv_unimodal_regression with weights
    y_cv_peak = -x_cv.^2 .+ 12.5 .* x_cv .+ rand(twister_cv, Normal(0, 0.5), obs_cv)
    cv_uni_w = cv_unimodal_regression(x_cv, y_cv_peak; nfolds=5, seed=42, weights=w_cv)
    @test cv_uni_w isa CVRegressionResult
    @test cv_uni_w.fitted isa Piecewise_Function
    @test cv_uni_w(2.5) isa Real

    # cv_shape_regression with weights
    cv_shape_w = cv_shape_regression(x_cv, y_cv_inc; shapes=:all, nfolds=5, seed=42, weights=w_cv)
    @test cv_shape_w isa CVRegressionResult
    @test length(cv_shape_w.cv_errors) == 6
    @test cv_shape_w.fitted isa Piecewise_Function

    # cv_shape_regression with weights and subset of shapes
    cv_shape_w2 = cv_shape_regression(x_cv, y_cv_inc; shapes=[:increasing, :quasiconcave], nfolds=5, seed=42, weights=w_cv)
    @test cv_shape_w2 isa CVRegressionResult
    @test length(cv_shape_w2.cv_errors) == 2

    # DataFrame interfaces pass weights through
    dd_cv = DataFrame(x = x_cv, y = y_cv_inc)
    cv_mono_df_w = cv_monotonic_regression(dd_cv, :x, :y; nfolds=5, seed=42, weights=w_cv)
    @test cv_mono_df_w isa CVRegressionResult

    cv_uni_df_w = cv_unimodal_regression(dd_cv, :x, :y; nfolds=5, seed=42, weights=w_cv)
    @test cv_uni_df_w isa CVRegressionResult

    cv_shape_df_w = cv_shape_regression(dd_cv, :x, :y; nfolds=5, seed=42, weights=w_cv)
    @test cv_shape_df_w isa CVRegressionResult

    println("Unimodal regression tests passed.")
end
