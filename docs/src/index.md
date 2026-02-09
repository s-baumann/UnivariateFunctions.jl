## Structs

There are four main `UnivariateFunction` structs that are part of this package. These are:
* `Undefined_Function` - An undefined function behaves similarly to "missing" in Julia. Whenever anything is added/multiplied/etc with an undefined function the result is undefined. The integral and derivative of an undefined function is undefined. If an undefined function is evaluated it will return a missing.
* `PE_Function` - This is the basic function type. It has a form of

$a \exp(b(x-base_)) (x-base)^d$

* `Sum_Of_Functions` - This is an array of `PE_Function`s. Note that by adding `PE_Function`s we can replicate any given polynomial. Hence from Weierstrass' approximation theorem we can approximate any continuous function on a bounded domain to any desired level of accuracy (whether this is practical in numerical computing depends on the function being approximated).
* `Piecewise_Function` - This defines a different `UnivariateFunction` for each part of the x domain.

It is possible to perform any additions, subtractions, multiplications between any two `UnivariateFunction`s and between Ints/Floats and any `UnivariateFunction`. No division is allowed and it is not possible to raise a `UnivariateFunction` to a negative power. This is to ensure that all `Univariatefunction`s are analytically integrable and differentiable. This may change in future releases.

## Interpolation and Splines
So far this package support the following interpolation schemes:
* Constant interpolation from the left to the right. Such a `Piecewise_Function` spline can be constructed by the `create_constant_interpolation_to_right` method.
* Constant interpolation from the right to the left. Such a `Piecewise_Function` spline can be constructed by the `create_constant_interpolation_to_left` method.
* Linear interpolation. Such a `Piecewise_Function` spline can be constructed by the `create_linear_interpolation` method.
It also supports the following spline (which can also be used for interpolation)
* `Schumaker` shape preserving spline - Such a `Piecewise_Function` spline can be constructed by the `create_quadratic_spline` method.

## Approximation and regression
This package supports the creation of the following approximation and regression schemes:
* OLS regression. The `create_ols_approximation` function can create a `UnivariateFunction` approximating a linear relationship. The degree input to this function can be used to specify the number of higher powers of x to be used in approximating y. For instance if the degree is two then y will be approximated as a linear combination of $x$ and $x^2$ as well as an intercept (if the intercept boolean is true).
* Chebyshev polynomials - This will approximate a function using the Chebyshev basis functions. This approximation function can then be integrated to accomplish Chebyshev–Gauss quadrature.

## Regression and Smoothing
The package provides several shape-constrained regression methods that return `Piecewise_Function` objects:
* `supersmoother` - Friedman's SuperSmoother (1984), an adaptive local linear regression that automatically selects bandwidth at each point.
* `isotonic_regression` - Fits a monotonic step function using the Pool Adjacent Violators algorithm.
* `monotonic_regression` - Fits a piecewise linear monotonic function using nonnegative least squares.
* `unimodal_regression` - Fits functions with a single peak (quasiconcave/concave) or trough (quasiconvex/convex).

Cross-validation functions are also provided for automatic shape selection:
* `cv_monotonic_regression` - Automatically selects between increasing and decreasing.
* `cv_unimodal_regression` - Automatically selects among the four unimodal shapes.
* `cv_shape_regression` - Selects from all six shapes (monotonic + unimodal) or a custom subset.

See the [Regression and Smoothing](Regression.md) page for detailed documentation and examples.

## Date Handling

* All base dates are immediately converted to floats and are not otherwise saved. Thus there is no difference between a `PE_Function` created with a base as a float and one created with the matching date. This is done to simplify the code. All date conversions is done by finding the year fractions between the date and the global base date of `ZonedDateTime(1970,1,1,0,0,0,tz"UTC")`. This particular global base date should not affect anything as long as it is consistent. It is relatively trivial to change it (in the `date_conversions.jl` file) and recompile however if desired.

## Major limitations
* It is not possible to divide by univariate functions or raise them by a negative power.
* When multiplying `PE_Function`s with different base dates there is often an issue of very high or very low numbers that go outside machine precision. If one were trying to change a `PE_Function` from base 2010 to 50, this would not generally be possible. This is because to change $a \exp(x-2020)$ to $q \exp(x - 50)$ we need to premultiply the first expression by $\exp(-1950)$ which is a tiny number. In these cases it is better to do the algebra on paper and rewriting the code accordingly as often base changes cancel out on paper. It is also good to change bases as rarely as possible. If different univariate functions use different bases then there is a need to base change when multiplying them which can result in errors. Note that if base changes are segment in the x domain by means of a piecewise function then they should never interact meaning it is ok to use different bases here.
* There is no support for finding optima, roots, fixedpoints etc. If anyone has an idea of how to do it efficiently then please let me know.
* There is no support for finding a function representing the upper/lower envelope of multiple functions. If anyone has an idea of how to do it efficiently then please let me know.
