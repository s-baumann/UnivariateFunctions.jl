# UnivariateFunctions.jl Examples
# ===============================
# This file demonstrates the regression and smoothing functionality,
# including supersmoother, monotonic regression, and unimodal regression.

using UnivariateFunctions
using DataFrames
using Random
using Distributions

# Set random seed for reproducibility
Random.seed!(42)

# =============================================================================
# Example 1: SuperSmoother on Noisy Sinusoidal Data
# =============================================================================

println("\n" * "="^60)
println("Example 1: SuperSmoother on Noisy Sinusoidal Data")
println("="^60)

n1 = 200
x1 = collect(range(0, 4π, length=n1))
y1_true = sin.(x1) .+ 0.5 .* cos.(2 .* x1)
y1 = y1_true .+ rand(Normal(0, 0.3), n1)

dd1 = DataFrame(x = x1, y = y1)

# Fit supersmoother
fit_ss = supersmoother(x1, y1)
println("SuperSmoother fit complete.")
println("Fit type: ", typeof(fit_ss))

# Plot supersmoother fit over data
plt_ss = plot(fit_ss, dd1; x_name=:x, y_name=:y)
# To display: plt_ss |> save("supersmoother_example.png")

# Try different smoothing levels
fit_ss_smooth = supersmoother(x1, y1; bass=8.0)  # More smoothing
fit_ss_rough = supersmoother(x1, y1; spans=[0.02, 0.05, 0.1])  # Less smoothing

# Compare multiple fits
plt_ss_compare = plot(
    Dict(:standard => fit_ss, :smooth => fit_ss_smooth, :rough => fit_ss_rough),
    0.0, 4π
)

println("SuperSmoother plots created.")

# =============================================================================
# Example 2: Monotonic Regression - Increasing Data
# =============================================================================

println("\n" * "="^60)
println("Example 2: Monotonic Regression - Increasing Data")
println("="^60)

n2 = 150
x2 = collect(range(0, 5, length=n2))
# Saturating growth curve (like Michaelis-Menten)
y2_true = 10 .* x2 ./ (2 .+ x2)
y2 = y2_true .+ rand(Normal(0, 0.5), n2)

dd2 = DataFrame(x = x2, y = y2)

# Fit monotonic increasing regression
fit_mono_inc = monotonic_regression(x2, y2; increasing=true, nbins=15)
println("Monotonic (increasing) fit complete.")

# Also fit with different settings
fit_mono_inc_fine = monotonic_regression(x2, y2; increasing=true, nbins=25)
fit_mono_inc_quantile = monotonic_regression(x2, y2; increasing=true, nbins=15, equally_spaced_bins=false)

# Plot
plt_mono_inc = plot(fit_mono_inc, dd2; x_name=:x, y_name=:y)

# Compare bin choices
plt_mono_compare = plot(
    Dict(:standard => fit_mono_inc, :fine => fit_mono_inc_fine, :quantile => fit_mono_inc_quantile),
    0.0, 5.0
)

println("Monotonic increasing plots created.")

# =============================================================================
# Example 3: Monotonic Regression - Decreasing Data
# =============================================================================

println("\n" * "="^60)
println("Example 3: Monotonic Regression - Decreasing Data")
println("="^60)

n3 = 150
x3 = collect(range(0, 10, length=n3))
# Exponential decay
y3_true = 10 .* exp.(-0.3 .* x3)
y3 = y3_true .+ rand(Normal(0, 0.5), n3)

dd3 = DataFrame(x = x3, y = y3)

# Fit monotonic decreasing regression
fit_mono_dec = monotonic_regression(x3, y3; increasing=false, nbins=12)
println("Monotonic (decreasing) fit complete.")

plt_mono_dec = plot(fit_mono_dec, dd3; x_name=:x, y_name=:y)

println("Monotonic decreasing plot created.")

# =============================================================================
# Example 4: CV-Selected Monotonic Regression
# =============================================================================

println("\n" * "="^60)
println("Example 4: CV-Selected Monotonic Regression")
println("="^60)

# Use cross-validation to automatically select increasing vs decreasing
cv_mono = cv_monotonic_regression(x2, y2; nbins=15, nfolds=10, seed=123)
println("CV Monotonic Regression complete.")
println("Selected shape: ", cv_mono.selected_shape)
println("CV errors: ", cv_mono.cv_errors)

plt_cv_mono = plot(cv_mono.fitted, dd2; x_name=:x, y_name=:y)

# =============================================================================
# Example 5: Unimodal Regression - Concave (Single Peak)
# =============================================================================

println("\n" * "="^60)
println("Example 5: Unimodal Regression - Concave (Single Peak)")
println("="^60)

n5 = 200
x5 = collect(range(-3, 3, length=n5))
# Quadratic with peak at x=0.5
y5_true = -2 .* (x5 .- 0.5).^2 .+ 5
y5 = y5_true .+ rand(Normal(0, 0.4), n5)

dd5 = DataFrame(x = x5, y = y5)

# Fit quasiconcave (single peak, slopes go + then -)
fit_qconcave = unimodal_regression(x5, y5; convex=false, quasi=true, nbins=12)
println("Quasiconcave fit complete.")

# Fit true concave (slopes monotonically decrease)
fit_concave = unimodal_regression(x5, y5; convex=false, quasi=false, nbins=12)
println("Concave fit complete.")

# Compare
plt_unimodal_peak = plot(
    Dict(:quasiconcave => fit_qconcave, :concave => fit_concave),
    -3.0, 3.0
)

plt_qconcave = plot(fit_qconcave, dd5; x_name=:x, y_name=:y)

println("Concave/quasiconcave plots created.")

# =============================================================================
# Example 6: Unimodal Regression - Convex (Single Trough)
# =============================================================================

println("\n" * "="^60)
println("Example 6: Unimodal Regression - Convex (Single Trough)")
println("="^60)

n6 = 200
x6 = collect(range(-2, 4, length=n6))
# Quadratic with trough at x=1
y6_true = 1.5 .* (x6 .- 1).^2 .+ 2
y6 = y6_true .+ rand(Normal(0, 0.5), n6)

dd6 = DataFrame(x = x6, y = y6)

# Fit quasiconvex (single trough, slopes go - then +)
fit_qconvex = unimodal_regression(x6, y6; convex=true, quasi=true, nbins=12)
println("Quasiconvex fit complete.")

# Fit true convex (slopes monotonically increase)
fit_convex = unimodal_regression(x6, y6; convex=true, quasi=false, nbins=12)
println("Convex fit complete.")

# Compare
plt_unimodal_trough = plot(
    Dict(:quasiconvex => fit_qconvex, :convex => fit_convex),
    -2.0, 4.0
)

plt_qconvex = plot(fit_qconvex, dd6; x_name=:x, y_name=:y)

println("Convex/quasiconvex plots created.")

# =============================================================================
# Example 7: CV-Selected Unimodal Regression
# =============================================================================

println("\n" * "="^60)
println("Example 7: CV-Selected Unimodal Regression")
println("="^60)

# Use CV to automatically select among quasiconcave/quasiconvex/concave/convex
cv_unimodal = cv_unimodal_regression(x5, y5; nbins=12, nfolds=10, seed=456)
println("CV Unimodal Regression complete.")
println("Selected shape: ", cv_unimodal.selected_shape)
println("CV errors: ", cv_unimodal.cv_errors)

plt_cv_unimodal = plot(cv_unimodal.fitted, dd5; x_name=:x, y_name=:y)

# =============================================================================
# Example 8: Full Shape Selection with cv_shape_regression
# =============================================================================

println("\n" * "="^60)
println("Example 8: Full Shape Selection with cv_shape_regression")
println("="^60)

# Test on clearly monotonic data
cv_shape_mono = cv_shape_regression(x2, y2; shapes=:all, nbins=12, nfolds=10, seed=789)
println("\nOn monotonic data:")
println("Selected shape: ", cv_shape_mono.selected_shape)
println("CV errors: ")
for (k, v) in sort(collect(cv_shape_mono.cv_errors), by=x->x[2])
    println("  $k: $(round(v, digits=2))")
end

# Test on unimodal data
cv_shape_uni = cv_shape_regression(x5, y5; shapes=:all, nbins=12, nfolds=10, seed=789)
println("\nOn unimodal (peak) data:")
println("Selected shape: ", cv_shape_uni.selected_shape)
println("CV errors: ")
for (k, v) in sort(collect(cv_shape_uni.cv_errors), by=x->x[2])
    println("  $k: $(round(v, digits=2))")
end

# =============================================================================
# Example 9: Comparing All Methods on Complex Data
# =============================================================================

println("\n" * "="^60)
println("Example 9: Comparing All Methods on Complex Data")
println("="^60)

n9 = 300
x9 = collect(range(0, 6, length=n9))
# Non-monotonic, non-symmetric data
y9_true = 3 .* sin.(x9) .+ 0.5 .* x9
y9 = y9_true .+ rand(Normal(0, 0.5), n9)

dd9 = DataFrame(x = x9, y = y9)

# Fit all methods
fit_ss9 = supersmoother(x9, y9)
fit_mono_inc9 = monotonic_regression(x9, y9; increasing=true, nbins=15)
fit_mono_dec9 = monotonic_regression(x9, y9; increasing=false, nbins=15)
fit_qconcave9 = unimodal_regression(x9, y9; convex=false, quasi=true, nbins=15)
fit_qconvex9 = unimodal_regression(x9, y9; convex=true, quasi=true, nbins=15)

# Compare all fits
plt_compare_all = plot(
    Dict(
        :supersmoother => fit_ss9,
        :mono_inc => fit_mono_inc9,
        :mono_dec => fit_mono_dec9,
        :quasiconcave => fit_qconcave9,
        :quasiconvex => fit_qconvex9
    ),
    0.0, 6.0
)

# Use CV to pick the best
cv_best = cv_shape_regression(x9, y9; shapes=:all, nbins=15, nfolds=10, seed=999)
println("Best shape via CV: ", cv_best.selected_shape)

plt_best = plot(cv_best.fitted, dd9; x_name=:x, y_name=:y)

println("Comparison plots created.")

# =============================================================================
# Example 10: Working with Derivatives and Integrals
# =============================================================================

println("\n" * "="^60)
println("Example 10: Working with Derivatives and Integrals")
println("="^60)

# Use the supersmoother fit from example 1
deriv_ss = derivative(fit_ss)
integ_ss = indefinite_integral(fit_ss)

println("Derivative type: ", typeof(deriv_ss))
println("Integral type: ", typeof(integ_ss))

# Plot original, derivative, and integral
plt_calculus = plot(
    Dict(:original => fit_ss, :derivative => deriv_ss),
    0.0, 4π
)

# Compute definite integral
area = evaluate_integral(fit_ss, 0.0, 2π)
println("Area under curve from 0 to 2π: ", round(area, digits=4))

# =============================================================================
# Summary of Available Plots
# =============================================================================

println("\n" * "="^60)
println("Summary of Available Plots")
println("="^60)
println("""
The following plot objects have been created:
- plt_ss: SuperSmoother on sinusoidal data
- plt_ss_compare: Comparing different smoothness levels
- plt_mono_inc: Monotonic increasing regression
- plt_mono_compare: Comparing bin choices for monotonic
- plt_mono_dec: Monotonic decreasing regression
- plt_cv_mono: CV-selected monotonic regression
- plt_qconcave: Quasiconcave (peak) regression with data
- plt_unimodal_peak: Comparing quasiconcave vs concave
- plt_qconvex: Quasiconvex (trough) regression with data
- plt_unimodal_trough: Comparing quasiconvex vs convex
- plt_cv_unimodal: CV-selected unimodal regression
- plt_compare_all: All methods on complex data
- plt_best: CV-selected best fit
- plt_calculus: Original function and its derivative

To display a plot in VS Code or a notebook, just evaluate the variable.
To save a plot: plt_ss |> VegaLite.save("filename.png")
""")

println("\nExamples complete!")
