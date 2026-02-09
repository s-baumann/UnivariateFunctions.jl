using Test
using Random

# Regression (snapshot) tests for numerical stability.
# These assert exact numerical outputs from each interpolation and approximation
# method to detect unintended changes during refactoring or optimisation.
#
# If a test fails after an intentional algorithmic change, update the reference
# values by running test/regression_snapshot.jl and copying the new outputs here.

@testset "Regression Snapshot Tests" begin
    using UnivariateFunctions

    # --- Fixed deterministic data (MersenneTwister seed 42) ---
    rng = MersenneTwister(42)
    n = 50
    x = collect(range(0.0, 10.0, length=n))
    y_linear = 2.0 .* x .+ 1.0 .+ 0.1 .* randn(rng, n)
    y_quad   = -1.0 .* (x .- 5.0).^2 .+ 25.0 .+ 0.1 .* randn(rng, n)
    y_convex = (x .- 5.0).^2 .+ 0.1 .* randn(rng, n)

    eval_pts = [0.0, 2.5, 5.0, 7.5, 10.0]

    # Helper: assert a function matches expected values at eval_pts
    function assert_values(f, expected; atol=0.0)
        for (p, exp) in zip(eval_pts, expected)
            @test evaluate(f, p) ≈ exp atol=atol
        end
    end

    # ---------------------------------------------------------------
    # Interpolation methods
    # ---------------------------------------------------------------

    @testset "create_linear_interpolation" begin
        f = create_linear_interpolation(x, y_linear)
        assert_values(f, [
            1.1210282022049036,
            6.0000948729052785,
            11.089147421599598,
            16.055947218083553,
            21.025357341664453,
        ])
        @test evaluate_integral(f, 0.0, 10.0) ≈ 110.08327223994316
    end

    @testset "create_quadratic_spline" begin
        f = create_quadratic_spline(x, y_linear)
        assert_values(f, [
            1.1210282022049036,
            6.024382500533199,
            11.095244309933515,
            16.066108596392173,
            21.025357341664453,
        ])
        @test evaluate_integral(f, 0.0, 10.0) ≈ 110.08270657572746
    end

    @testset "create_constant_interpolation_to_right" begin
        f = create_constant_interpolation_to_right(x, y_linear)
        assert_values(f, [
            1.1210282022049036,
            5.941915965322743,
            10.934090842029162,
            15.696011300829854,
            21.025357341664453,
        ])
        @test evaluate_integral(f, 0.0, 10.0) ≈ 108.05221824612073
    end

    @testset "create_constant_interpolation_to_left" begin
        f = create_constant_interpolation_to_left(x, y_linear)
        assert_values(f, [
            1.400263404867264,
            6.174631595652884,
            11.244204001170035,
            16.175925857168117,
            21.025357341664453,
        ])
        @test evaluate_integral(f, 0.0, 10.0) ≈ 112.1143262337656
    end

    # ---------------------------------------------------------------
    # Approximation methods
    # ---------------------------------------------------------------

    @testset "create_ols_approximation (degree=1)" begin
        f = create_ols_approximation(y_linear, x, 0.0, 1)
        assert_values(f, [
            0.995648431126949,
            6.0026364830400345,
            11.009624534953119,
            16.016612586866206,
            21.02360063877929,
        ])
    end

    @testset "create_ols_approximation (degree=2)" begin
        f = create_ols_approximation(y_linear, x, 0.0, 2)
        assert_values(f, [
            1.0177485413129332,
            5.999355997934307,
            10.99788385141682,
            16.013332101760472,
            21.045700748965263,
        ])
    end

    @testset "create_chebyshev_approximation" begin
        f = create_chebyshev_approximation(sin, 12, 8, 0.0, 6.0)
        cheb_pts = [0.0, 1.0, 3.14159, 5.0, 6.0]
        cheb_expected = [
            -0.0002992844276874518,
            0.841293049071965,
            -5.942107383071188e-5,
            -0.9590014607732853,
            -0.2793875808542807,
        ]
        for (p, exp) in zip(cheb_pts, cheb_expected)
            @test evaluate(f, p) ≈ exp
        end
    end

    # ---------------------------------------------------------------
    # Regression methods
    # ---------------------------------------------------------------

    @testset "isotonic_regression" begin
        f = isotonic_regression(x, y_linear; increasing=true)
        assert_values(f, [
            1.1210282022049036,
            6.0000948729052785,
            11.089147421599598,
            16.055947218083553,
            21.025357341664453,
        ])
    end

    @testset "monotonic_regression" begin
        f = monotonic_regression(x, y_linear; nbins=10, increasing=true)
        assert_values(f, [
            1.0551437521045872,
            5.9342444089046875,
            11.101456854136567,
            15.995765286718186,
            21.059834889604364,
        ])
        @test evaluate_integral(f, 0.0, 10.0) ≈ 110.08627351913509
    end

    @testset "supersmoother" begin
        f = supersmoother(x, y_linear)
        assert_values(f, [
            1.0602460860640845,
            5.9594326614265345,
            11.065371114572294,
            16.00511913708144,
            21.05946737357979,
        ])
    end

    @testset "unimodal_regression (quasiconcave)" begin
        f = unimodal_regression(x, y_quad; nbins=10, convex=false, quasi=true)
        assert_values(f, [
            0.09776664977527401,
            18.66871851055128,
            25.21862782415963,
            18.63405971006622,
            0.10956074262233706,
        ])
    end

    @testset "unimodal_regression (quasiconvex)" begin
        f = unimodal_regression(x, y_convex; nbins=10, convex=true, quasi=true)
        assert_values(f, [
            24.976898904683974,
            6.317190969547771,
            -0.10943656117559897,
            6.350972845035017,
            24.75688214498734,
        ])
    end

    # ---------------------------------------------------------------
    # Cross-validated regression methods
    # ---------------------------------------------------------------

    @testset "cv_monotonic_regression" begin
        result = cv_monotonic_regression(x, y_linear; nbins=10, nfolds=5, seed=123)
        @test result.selected_shape == :increasing
        assert_values(result.fitted, [
            1.0551437521045872,
            5.9342444089046875,
            11.101456854136567,
            15.995765286718186,
            21.059834889604364,
        ])
    end

    @testset "cv_unimodal_regression" begin
        result = cv_unimodal_regression(x, y_quad; nbins=10, nfolds=5, seed=123)
        @test result.selected_shape == :quasiconcave
        assert_values(result.fitted, [
            0.09776664977527401,
            18.66871851055128,
            25.21862782415963,
            18.63405971006622,
            0.10956074262233706,
        ])
    end

    @testset "cv_shape_regression" begin
        result = cv_shape_regression(x, y_linear; shapes=:all, nbins=10, nfolds=5, seed=123)
        @test result.selected_shape == :convex
        assert_values(result.fitted, [
            1.0680937913648085,
            5.976164630231834,
            10.995403725144863,
            16.01464282005789,
            21.070994992060847,
        ])
    end

    # ---------------------------------------------------------------
    # Simplification
    # ---------------------------------------------------------------

    @testset "simplify" begin
        f_lin = create_linear_interpolation(x, y_linear)

        f_s_lin = simplify(f_lin, 20, 0.0, 10.0, :linear)
        assert_values(f_s_lin, [
            1.1210282022049036,
            5.904876528813651,
            11.048347587826495,
            16.014934593373425,
            21.025357341664453,
        ])

        f_s_quad = simplify(f_lin, 20, 0.0, 10.0, :quadratic)
        assert_values(f_s_quad, [
            1.1210282022049036,
            5.899033011726793,
            11.041038513948381,
            16.007266197954852,
            21.025357341664453,
        ])
    end

    # ---------------------------------------------------------------
    # UnivariateFitter
    # ---------------------------------------------------------------

    @testset "UnivariateFitter" begin
        fitter = UnivariateFitter(:increasing; nbins=10, weight_on_new=1.0)
        UnivariateFunctions.fit!(fitter, x, y_linear)
        assert_values(fitter, [
            1.0551437521045872,
            5.9342444089046875,
            11.101456854136567,
            15.995765286718186,
            21.059834889604364,
        ])
    end
end
