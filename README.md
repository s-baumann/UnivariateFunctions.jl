# UnivariateFunctions.jl

The build status is https://travis-ci.com/s-baumann/UnivariateFunctions.jl.svg?branch=master

This implements single algebra and evaluation on simple univariate functions.
There are a few ways in which it can be used.
* As in the StochasticIntegrals.jl package this package can be used to define
    functions that will be the integrands in stochastic integrals. This has the beneft
    that the mean and variances implied by these stochastic integrals can be found analytically.
* UnivariateFunctions can be used in the creation of splines. This has the added
    advantage that a spline implemented as a UnivariateFunction can be manipulated
    easily. It can be differentiated and then added to another function, etc.
* Any continuous interpolation scheme can be accomplished with the added benefit that derivatives/integrals/products etc can be found analytically.

## Structs

There are four main UnivariateFunction structs that are part of this package. These are:
* Undefined_Function - An undefined function behaves similarly to "missing" in Julia. Whenever anything is added/multiplied/etc with an undefined function the result is undefined. The integral and derivative of an undefined function is undefined. If an undefined function is evaluated it will throw an error.
* PE_Function - This is the basic function type. It has a form of $a \exp(b(x-base_)) (x-base)^d$.
* Sum_Of_Functions - This is an array of PE_Functions. Note that by adding PE_Functions we can replicate any given polynomial. Hence from Weierstrass' approximation theorem we can approximate any continuous function on a bounded domain to any desired level of accuracy.[^1]
* Piecewise_Function - This defines a different UnivariateFunction for each part of the x domain.

1: The degree to which this is numerically practical will depend on the function in question. Note however that it is possible to use a Piecewise_Function to define a different polynomial for each part of the x domain. PE_Function also contains an exponential term. With these features it should be possible to approximate the majority of functions fairly well.

It is possible to perform any additions, subtractions, multiplications between any two UnivariateFunctions and between Ints/Floats and any UnivariateFunction. No division is allowed and it is not possible to raise a UnivariateFunction to a negative power. This is to ensure that all univariatefunctions are analytically integrable and differentiable. This may change in future releases.

## Major limitations
* It is not possible to divide by univariate functions or raise them by a negative power.
* When multiplying pe_functions with different base dates there is often an issue of very high or very low numbers that go outside machine precision. If one were trying to change a PE_Function from base 2010 to 50, this would not generally be possible. This is because to change $a exp(x-2020)$ to $q exp(x - 50)$ we need to premultiply the first expression by $exp(-1950)$ which is a tiny number. In these cases it is better to do the algebra on paper and rewriting the code accordingly as often base changes cancel out on paper. It is also good to change bases as rarely as possible. If different univariate functions use different bases then there is a need to base change when multiplying them which can result in errors. Note that if base changes are segment in the x domain by means of a piecewise function then they should never interact meaning it is ok to use different bases here.
* There is no support for finding optima, roots, fixedpoints etc. If anyone has an idea of how to do it efficiently then please let me know.
* There is no support for finding a function representing the upper/lower envelope of multiple functions. If anyone has an idea of how to do it efficiently then please let me know.

## Interpolation and Splines
So far this package support the following interpolation schemes:
* Constant interpolation from the left to the right. Such a Piecewise_Function spline can be constructed by the create_constant_interpolation_to_right method.
* Constant interpolation from the right to the left. Such a Piecewise_Function spline can be constructed by the create_constant_interpolation_to_left method.
* Linear interpolation. Such a Piecewise_Function spline can be constructed by the create_linear_interpolation method.
It also supports the following spline (which can also be used for interpolation)
* Schumaker shape preserving spline - Such a Piecewise_Function spline can be constructed by the create_quadratic_spline method.

## Date Handling

* All base dates are immediately converted to floats and are not otherwise saved. Thus there is no difference between a PE_Function created with a base as a float and one created with the matching date. This is done to simplify the code. All date conversions is done by finding the year fractions between the date and the global base date of Date(2000,1,1). This particular global base date should not affect anything as long as it is consistent. It is relatively trivial to change it (in the helpers.jl file) and recompile however if desired.

## TODO (help wanted if anyone feels keen)

* Implement more useful splines and interpolation schemes.
* Get the constructor of sum_of_functions to amalgamate PE_Functions with the same b_, base_ and d_ values. This may involve some base changes in order to get them to match.

# Examples

## For basic algebra:

Consider we have a two functions f and g and want to add them, multiply them by some other function h, then square it and finally integrate the result between 2.0 and 2.8. This can be done analytically with UnivariateFunctions:
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
Now we can generate a UnivariateFunction that can be used to easily interpolate from the sampled points:
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
If we had wanted to interpolate instead with a constant method(from left or from right) or by linearly
interpolating then we could have just generated func with a different method:
create_constant_interpolation_to_left,
create_constant_interpolation_to_right or
create_linear_interpolation.
