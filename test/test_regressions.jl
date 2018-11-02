using UnivariateFunctions
using Dates
using Random
using GLM

tol = 10*eps()

Random.seed!(1)
obs = 1000000
X = rand(obs)
y = X .+ rand(Normal(),obs) .+ 7

# Basic use case with 2 degrees
lm1 = fit(LinearModel, hcat(ones(obs), X, X .^ 2), y)
preds = predict(lm1, hcat(ones(obs), X, X .^ 2))

package_approximation = UnivariateFunctions.create_ols_approximation(y, X, 0.0, 2, true)
package_predictions = UnivariateFunctions.evaluate.(package_approximation, X)
sum(abs.(preds .- package_predictions)) < 1e-10

# Degree of 1 with no intercept
lm1 = fit(LinearModel, hcat(X), y)
preds = predict(lm1, hcat(X))

package_approximation = UnivariateFunctions.create_ols_approximation(y, X, 0.0, 1, false)
package_predictions = UnivariateFunctions.evaluate.(package_approximation, X)
sum(abs.(preds .- package_predictions)) < 1e-10

# Degree of 1 with no intercept
lm1 = fit(LinearModel, hcat(ones(obs)), y)
preds = predict(lm1, hcat(ones(obs)))

package_approximation = UnivariateFunctions.create_ols_approximation(y, X, 0.0, 0, true)
package_predictions = UnivariateFunctions.evaluate.(package_approximation, X)
sum(abs.(preds .- package_predictions)) < 1e-10
