# Regression and Smoothing

This package provides several regression and smoothing methods that fit `Piecewise_Function` objects to data. Because the fitted functions are `UnivariateFunction`s, you can differentiate, integrate, and combine them analytically.

## Overview of Methods

| Function | Description | Shape Constraint |
|----------|-------------|------------------|
| `supersmoother` | Adaptive local linear regression | None (flexible) |
| `isotonic_regression` | Pool Adjacent Violators algorithm | Monotonic |
| `monotonic_regression` | Binned monotonic regression | Monotonic (piecewise linear) |
| `unimodal_regression` | Single peak/trough regression | Unimodal or curvature-constrained |
| `cv_monotonic_regression` | Auto-select increasing/decreasing | Monotonic with CV |
| `cv_unimodal_regression` | Auto-select unimodal shape | Unimodal with CV |
| `cv_shape_regression` | Auto-select from all shapes | Any shape with CV |

## SuperSmoother

Friedman's SuperSmoother (1984) is a local linear regression with adaptive bandwidth selection. It automatically chooses the amount of smoothing at each point based on the local structure of the data.

```julia
using UnivariateFunctions
using Random, Distributions

# Generate noisy sinusoidal data
Random.seed!(42)
x = collect(range(0, 4π, length=200))
y = sin.(x) .+ rand(Normal(0, 0.3), 200)

# Fit supersmoother
fit = supersmoother(x, y)

# Evaluate at new points
fit(2.5)
fit.(collect(0:0.1:4π))

# Get derivative and integral
deriv = derivative(fit)
integ = indefinite_integral(fit)
```

### Parameters

- `spans`: Vector of candidate spans (fraction of data). Default `[0.05, 0.2, 0.5]`
- `bass`: Bass enhancement (0-10). Higher values favor smoother fits. Default `0.0`

```julia
# More smoothing with bass enhancement
fit_smooth = supersmoother(x, y; bass=8.0)

# Less smoothing with smaller spans
fit_rough = supersmoother(x, y; spans=[0.02, 0.05, 0.1])
```

## Isotonic Regression

Isotonic regression fits a monotonically increasing (or decreasing) step function using the Pool Adjacent Violators (PAV) algorithm. This is useful when you know the relationship should be monotonic but don't want to assume a parametric form.

```julia
# Generate increasing data with noise
x = collect(range(0, 5, length=100))
y = 2.0 .* x .+ rand(Normal(0, 1.0), 100)

# Fit isotonic increasing
fit_inc = isotonic_regression(x, y; increasing=true)

# Fit isotonic decreasing
fit_dec = isotonic_regression(x, y; increasing=false)
```

## Monotonic Regression

`monotonic_regression` fits a piecewise linear function with monotonicity constraints using nonnegative least squares. This gives smoother results than `isotonic_regression` while still enforcing monotonicity.

```julia
# Fit with 15 bins, equally spaced
fit = monotonic_regression(x, y; nbins=15, increasing=true)

# Use observation-based bins (quantile spacing)
fit_quantile = monotonic_regression(x, y; nbins=15, equally_spaced_bins=false)
```

### Parameters

- `nbins`: Number of bins for piecewise linear fit. Default `10`
- `equally_spaced_bins`: If `true`, bins are equally spaced in x; if `false`, based on observation quantiles. Default `true`
- `increasing`: If `true`, fit increasing function; if `false`, decreasing. Default `true`

## Unimodal Regression

`unimodal_regression` fits functions with a single peak (concave/quasiconcave) or single trough (convex/quasiconvex). This is useful for data with a clear maximum or minimum.

### Shape Options

The function supports four shape constraints controlled by `convex` and `quasi` parameters:

| `convex` | `quasi` | Shape | Description |
|----------|---------|-------|-------------|
| `false` | `true` | Quasiconcave | Single peak: slopes go + then - |
| `true` | `true` | Quasiconvex | Single trough: slopes go - then + |
| `false` | `false` | Concave | Slopes monotonically decrease |
| `true` | `false` | Convex | Slopes monotonically increase |

```julia
# Data with a single peak
x = collect(range(-3, 3, length=200))
y = -2 .* (x .- 0.5).^2 .+ 5 .+ rand(Normal(0, 0.4), 200)

# Fit quasiconcave (single peak, no curvature constraint)
fit_qc = unimodal_regression(x, y; convex=false, quasi=true)

# Fit true concave (slopes must decrease)
fit_c = unimodal_regression(x, y; convex=false, quasi=false)

# Data with a single trough
y_trough = 2 .* (x .- 0.5).^2 .+ 1 .+ rand(Normal(0, 0.4), 200)

# Fit quasiconvex (single trough)
fit_qv = unimodal_regression(x, y_trough; convex=true, quasi=true)

# Fit true convex (slopes must increase)
fit_v = unimodal_regression(x, y_trough; convex=true, quasi=false)
```

### Parameters

- `nbins`: Number of bins for piecewise linear fit. Default `10`
- `equally_spaced_bins`: If `true`, bins are equally spaced; if `false`, quantile-based. Default `true`
- `convex`: If `false`, fit peak (concave); if `true`, fit trough (convex). Default `false`
- `quasi`: If `true`, only enforce unimodality; if `false`, also enforce curvature. Default `true`

## Cross-Validation Model Selection

When you're unsure about the shape of the relationship, use cross-validation to automatically select the best shape.

### cv_monotonic_regression

Automatically selects between increasing and decreasing monotonic functions:

```julia
# Let CV decide if relationship is increasing or decreasing
result = cv_monotonic_regression(x, y; nbins=15, nfolds=10, seed=42)

# Access the fitted function
result.fitted(2.5)

# See which shape was selected
result.selected_shape  # :increasing or :decreasing

# Compare CV errors
result.cv_errors  # Dict(:increasing => ..., :decreasing => ...)
```

### cv_unimodal_regression

Automatically selects among the four unimodal shapes:

```julia
result = cv_unimodal_regression(x, y; nbins=12, nfolds=10, seed=42)

result.selected_shape  # :quasiconcave, :quasiconvex, :concave, or :convex
result.cv_errors       # Dict with all four CV errors
```

### cv_shape_regression

The most flexible option - selects from all six shapes (2 monotonic + 4 unimodal):

```julia
# Choose from all shapes
result = cv_shape_regression(x, y; shapes=:all, nbins=12, nfolds=10, seed=42)

# Or restrict to specific categories
result_mono = cv_shape_regression(x, y; shapes=:monotonic)  # increasing/decreasing only
result_uni = cv_shape_regression(x, y; shapes=:unimodal)    # unimodal shapes only

# Or specify exactly which shapes to consider
result_custom = cv_shape_regression(x, y; shapes=[:increasing, :quasiconcave, :convex])
```

### CVRegressionResult

All CV functions return a `CVRegressionResult` struct:

```julia
struct CVRegressionResult
    fitted::Piecewise_Function  # The fitted function
    selected_shape::Symbol       # Which shape was selected
    cv_errors::Dict{Symbol, Float64}  # CV error for each candidate
    nfolds::Int                  # Number of folds used
end
```

The result is callable - you can use it directly as a function:

```julia
result = cv_shape_regression(x, y)
result(2.5)  # Evaluates the fitted function at 2.5
```

## DataFrame Interface

All regression functions support a DataFrame interface:

```julia
using DataFrames

df = DataFrame(x = x, y = y)

# All functions accept DataFrame + column symbols
fit = supersmoother(df, :x, :y)
fit = monotonic_regression(df, :x, :y; increasing=true)
fit = unimodal_regression(df, :x, :y; convex=false)
fit = cv_shape_regression(df, :x, :y; shapes=:all)
```

## Working with Fitted Functions

Since all fitted functions are `Piecewise_Function`s (which are `UnivariateFunction`s), you get full access to analytical operations:

```julia
fit = supersmoother(x, y)

# Derivative - returns another Piecewise_Function
deriv = derivative(fit)
deriv(2.5)  # Evaluate derivative at a point

# Indefinite integral
integ = indefinite_integral(fit)

# Definite integral
area = evaluate_integral(fit, 0.0, 5.0)

# Combine with other functions
g = PE_Function(2.0, 0.0, 0.0, 1)  # g(x) = 2x
combined = fit + g
combined = fit * g
combined = fit^2

# Plot (requires VegaLite)
plot(fit, 0.0, 5.0)
```

## Complete Example

```julia
using UnivariateFunctions
using DataFrames, Random, Distributions

Random.seed!(123)

# Generate data with unknown shape
n = 200
x = collect(range(0, 6, length=n))
y = 3 .* sin.(x) .+ 0.3 .* x .+ rand(Normal(0, 0.5), n)

df = DataFrame(x = x, y = y)

# Use CV to find the best shape
result = cv_shape_regression(x, y; shapes=:all, nbins=15, nfolds=10, seed=42)

println("Selected shape: ", result.selected_shape)
println("CV errors:")
for (shape, err) in sort(collect(result.cv_errors), by=x->x[2])
    println("  $shape: $(round(err, digits=2))")
end

# Use the fitted function
fit = result.fitted

# Compute derivative to find rate of change
deriv = derivative(fit)

# Find total area under curve
area = evaluate_integral(fit, 0.0, 6.0)
println("Area under curve: ", round(area, digits=4))

# Plot fit over data
plt = plot(fit, df; x_name=:x, y_name=:y)
```
