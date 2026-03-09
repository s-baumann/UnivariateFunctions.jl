using Test

@testset "Regression Tests" begin
    using UnivariateFunctions
    using Dates
    using Random
    using GLM
    using StatsBase

    tol = 10*eps()

    Random.seed!(1)
    obs = 1000
    X = rand(obs)
    y = X .+ rand(Normal(),obs) .+ 7

    # Basic use case with 2 degrees
    lm1 = fit(LinearModel, hcat(ones(obs), X, X .^ 2), y)
    glm_preds = predict(lm1, hcat(ones(obs), X, X .^ 2))
    package_approximation = UnivariateFunctions.create_ols_approximation(y, X, 0.0, 2, true)
    package_predictions = UnivariateFunctions.evaluate.(package_approximation, X)
    @test sum(abs.(glm_preds .- package_predictions)) < 1e-10

    # Degree of 1 with no intercept
    lm1 = fit(LinearModel, hcat(X), y)
    glm_preds = predict(lm1, hcat(X))
    package_approximation = UnivariateFunctions.create_ols_approximation(y, X, 0.0, 1, false)
    package_predictions = UnivariateFunctions.evaluate.(package_approximation, X)
    @test sum(abs.(glm_preds .- package_predictions)) < 1e-10

    # Degree of 1 with no intercept
    lm1 = fit(LinearModel, hcat(ones(obs)), y)
    glm_preds = predict(lm1, hcat(ones(obs)))
    package_approximation = UnivariateFunctions.create_ols_approximation(y, X, 0.0, 0, true)
    package_predictions = UnivariateFunctions.evaluate.(package_approximation, X)
    @test sum(abs.(glm_preds .- package_predictions)) < 1e-10

    # With applying a base
    Xmin = X .- 5
    lm1 = fit(LinearModel, hcat(ones(obs), Xmin, Xmin .^ 2), y)
    glm_preds = predict(lm1, hcat(ones(obs), Xmin, Xmin .^ 2))
    package_approximation = UnivariateFunctions.create_ols_approximation(y, X, 5.0, 2, true)
    package_predictions = UnivariateFunctions.evaluate.(package_approximation, X)
    @test sum(abs.(glm_preds .- package_predictions)) < 1e-10

    # With x as dates and a base.
    X = Array{Date}(undef, obs)
    baseDate =  Date(2016, 7, 21)
    StartDate = Date(2018, 7, 21)
    for i in 1:obs
        X[i] = StartDate +Dates.Day(2* (i-1))
    end
    XConverted = years_between.(X, baseDate)
    lm1 = fit(LinearModel, hcat(ones(obs), XConverted, XConverted .^ 2), y)
    glm_preds = predict(lm1, hcat(ones(obs), XConverted, XConverted .^ 2))
    package_approximation = UnivariateFunctions.create_ols_approximation(y, X, baseDate, 2, true)
    package_predictions = UnivariateFunctions.evaluate.(package_approximation, X)
    @test sum(abs.(glm_preds .- package_predictions)) < 1e-10


    # Chebyshev approximation
    func = sin
    nodes  =  12
    degree =  8
    left   = -2.0
    right  =  5.0
    approxim = UnivariateFunctions.create_chebyshev_approximation(func, nodes, degree, left, right)
    X = convert(Array{Float64,1}, left:0.01:right)
    y = func.(X)
    y_approx = evaluate.(Ref(approxim), X)
    @test maximum(abs.(y .- y_approx)) < 0.01

    func = exp
    nodes  =  12
    degree =  8
    left   =  1.0
    right  =  5.0
    approxim = UnivariateFunctions.create_chebyshev_approximation(func, nodes, degree, left, right)
    X = convert(Array{Float64,1}, left:0.01:right)
    y = func.(X)
    y_approx = evaluate.(Ref(approxim), X)
    @test maximum(abs.(y .- y_approx)) < 0.01

    # Getting Chebyshevs
    first_kinds = get_chebyshevs_up_to(4,true)
    @test length(first_kinds) == 4
    second_kinds = get_chebyshevs_up_to(4,false)
    @test length(second_kinds) == 4

    first_kinds = get_chebyshevs_up_to(25,true)
    @test length(first_kinds) == 25
    second_kinds = get_chebyshevs_up_to(25,false)
    @test length(second_kinds) == 25

    # ========== Weighted OLS approximation ==========

    Random.seed!(42)
    obs_w = 200
    X_w = rand(obs_w) .* 5
    y_w = 2.0 .* X_w .+ 3.0 .+ rand(Normal(), obs_w) .* 0.5

    # Equal weights should match unweighted
    ols_unweighted = UnivariateFunctions.create_ols_approximation(y_w, X_w, 0.0, 1, true)
    ols_equal_w = UnivariateFunctions.create_ols_approximation(y_w, X_w, 0.0, 1, true; weights=ones(obs_w))
    @test abs(evaluate(ols_unweighted, 2.5) - evaluate(ols_equal_w, 2.5)) < 1e-8

    # Non-uniform weights
    w_ols = rand(obs_w) .+ 0.1
    ols_weighted = UnivariateFunctions.create_ols_approximation(y_w, X_w, 0.0, 1, true; weights=w_ols)
    @test ols_weighted isa Sum_Of_Functions
    # Weighted OLS should give a reasonable fit
    pred_w = evaluate(ols_weighted, 2.5)
    @test abs(pred_w - (2.0 * 2.5 + 3.0)) < 2.0

    # Verify weighted matches GLM directly
    X_design = hcat(ones(obs_w), X_w)
    lm_weighted = fit(LinearModel, X_design, y_w; weights=StatsBase.FrequencyWeights(Float64.(w_ols)))
    glm_preds_w = predict(lm_weighted, X_design)
    pkg_preds_w = UnivariateFunctions.evaluate.(ols_weighted, X_w)
    @test sum(abs.(glm_preds_w .- pkg_preds_w)) < 1e-8

    # Higher degree with weights
    ols_weighted_quad = UnivariateFunctions.create_ols_approximation(y_w, X_w, 0.0, 2, true; weights=w_ols)
    @test ols_weighted_quad isa Sum_Of_Functions

    # No intercept with weights
    ols_weighted_noint = UnivariateFunctions.create_ols_approximation(y_w, X_w, 0.0, 1, false; weights=w_ols)
    @test ols_weighted_noint isa Sum_Of_Functions

    # With base_x and weights
    ols_weighted_base = UnivariateFunctions.create_ols_approximation(y_w, X_w, 2.5, 2, true; weights=w_ols)
    @test ols_weighted_base isa Sum_Of_Functions

    # Date interface with weights
    X_dates = [Date(2020, 1, 1) + Dates.Day(i) for i in 1:obs_w]
    ols_date_w = UnivariateFunctions.create_ols_approximation(y_w, X_dates; weights=w_ols)
    @test ols_date_w isa Sum_Of_Functions
end
