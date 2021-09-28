
# Examples

## For basic algebra:

Consider we have a two functions `f` and `g` and want to add them, multiply them by some other function `h`, then square it and finally integrate the result between 2.0 and 2.8. This can be done analytically with `UnivariateFunction`s:
```
f = PE_Function(1.0, 2.0, 4.0, 5)
g = PE_Function(1.3, 2.0, 4.3, 2)
h = PE_Function(5.0, 2.2, 1.0,0)
result_of_operations = (h*(f+g))^2
evaluate_integral(result_of_operations, 2.0, 2.8)
```

## For data interpolation

Suppose we have want to approximate some function with some sampled points. First to generate some points
```
using UnivariateFunctions
const global_base_date = Date(2000,1,1)
StartDate = Date(2018, 7, 21)
x = Array{Date}(undef, 1000)
for i in 1:1000
    x[i] = StartDate +Dates.Day(2* (i-1))
end
function ff(x::Date)
    days_between = years_from_global_base(x)
    return log(days_between) + sqrt(days_between)
end
y = ff.(x)
```
Now we can generate a `UnivariateFunction` that can be used to easily interpolate from the sampled points:
```
func = create_quadratic_spline(x,y)
```
And we can evaluate from this function and integrate it and differentiate it in the normal way:
```
evaluate(func, Date(2020,1,1))
evaluate.(Ref(func), [Date(2020,1,1), Date(2021,1,2)])
evaluate(derivative(func), Date(2021,1,2))
evaluate_integral(func, Date(2020,1,1), Date(2021,1,2))
```
If we had wanted to interpolate instead with a constant method (from left or from right) or by linearly interpolating then we could have just generated func with a different method:
* `create_constant_interpolation_to_left`,
* `create_constant_interpolation_to_right`,
* `create_linear_interpolation`.

If we have lots of data that we want to summarise with OLS
```
# Generating example data
using Random
Random.seed!(1)
obs = 1000
X = rand(obs)
y = X .+ rand(Normal(),obs) .+ 7
# And now making an approximation function
approxFunction = create_ols_approximation(y, X, 0.0, 2, true)
```
And if we want to approximate the sin function in the [2.3, 5.6] bound with 7 polynomial terms and 20 approximation nodes:
```
chebyshevs = create_chebyshev_approximation(sin, 20, 7, 2.3, 5.6)
```
We can integrate the above term in the normal way to achieve Gauss-Chebyshev quadrature:
```
evaluate_integral(chebyshevs, 2.3, 5.6)
```
