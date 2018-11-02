using UnivariateFunctions
using Dates
using Random
using GLM

tol = 10*eps()

Random.seed!(1)
obs = 1000000
X = rand(obs)
y = X .+ rand(Normal(),obs) .+ 7

lm1 = fit(LinearModel, hcat(ones(obs), X, X .^ 2), y)
preds = predict(lm1, hcat(ones(obs), X, X .^ 2))

package_approximation = create_ols_approximation(y, X, 0.0, 2, true)
