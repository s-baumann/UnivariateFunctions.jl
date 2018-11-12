using DataFrames
using Random
using Distributions
using GLM
using UnivariateFunctions

Random.seed!(1992)
nObs = 1000
dd = DataFrame()
dd[:x] = rand( Normal(),nObs) + 0.1 .* rand( Normal(),nObs)
dd[:z] = rand( Normal(),nObs) + 0.1 .* rand( Normal(),nObs)
dd[:w] = (0.5 .* rand( Normal(),nObs)) .+ 0.7.*(dd[:z] .- dd[:x]) + 0.1 .* rand( Normal(),nObs)
dd[:y] = (dd[:x] .*dd[:w] ) .* (dd[:z] .- dd[:w]) .+ dd[:x] + rand( Normal(),nObs)
dd[7,:y] = 1.0

constant = PE_Function(1.0,0.0,0.0,0)
linear  = PE_Function(1.0,0.0,0.0,1)
quadratic = PE_Function(1.0,0.0,0.0,2)

constant_term = Multivariate_PE_Function(1.0, Dict{Symbol,PE_Function}())
x_lin         = Multivariate_PE_Function(1.0, Dict{Symbol,PE_Function}(:x .=> linear))
z_lin         = Multivariate_PE_Function(1.0, Dict{Symbol,PE_Function}(:z .=> linear))
w_lin         = Multivariate_PE_Function(1.0, Dict{Symbol,PE_Function}(:w .=> linear))
w_quad        = Multivariate_PE_Function(1.0, Dict{Symbol,PE_Function}(:w .=> quadratic))
x_lin_z_quad  = Multivariate_PE_Function(1.0, Dict{Symbol,PE_Function}([:x, :z] .=> [linear, quadratic]))
model = constant_term + x_lin + z_lin + w_lin + w_quad + x_lin_z_quad
mod_1, reg_1 = create_ols_approximation(dd, :y, model)
sum(abs.(evaluate(mod_1, dd) - reg_1.rr.mu)) < 1e-12

#mod_2, reg_2 = create_saturated_ols_approximation(dd, :y, [:x, :w, :z], 3, true)
