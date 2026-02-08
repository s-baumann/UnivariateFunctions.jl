```@meta
CurrentModule = UnivariateFunctions
```

# Internal Functions

```@index
Pages = ["api.md"]
```

### Main Structs

```@docs
UnivariateFunction
Undefined_Function
PE_Function
Sum_Of_Functions
Piecewise_Function
```

### Function Evaluation and Calculus

Note that in addition to the below functions the following operators:

`+`, `-`, `\`, `*`, `^`

have also been overloaded so that a function will be returned with the analytical
sum, difference, product, quotient, power. The restrictions are that you cannot
divide by a function (although you can divide by a scalar) and only positive
integer powers can be taken.

```@docs
evaluate
derivative
indefinite_integral
evaluate_integral
right_integral
left_integral
```

### Interpolation and Simplification

```@docs
create_quadratic_spline
create_constant_interpolation_to_right
create_constant_interpolation_to_left
create_linear_interpolation
simplify
```

### Iterative Fitting

```@docs
UnivariateFitter
fit!
```

### Approximation

```@docs
create_ols_approximation
create_chebyshev_approximation
```

### Smoothing

```@docs
supersmoother
```

### Monotonic Regression

```@docs
isotonic_regression
monotonic_regression
```

### Unimodal Regression

```@docs
unimodal_regression
```

### Cross-Validation Model Selection

```@docs
cv_monotonic_regression
cv_unimodal_regression
cv_shape_regression
CVRegressionResult
```

### Internal Functions

```@docs
change_base_of_PE_Function
trim_piecewise_function
convert_to_linearly_rescale_inputs
get_chevyshevs_up_to
_kfold_indices
_cv_error
_local_linear_smooth
_local_linear_loo_residuals
_smooth_values
_supersmoother_values
```

### Date Conversions

```@docs
seconds_between
days_between
years_between
period_length
years_from_global_base_date
zdt2unix
unix2zdt
unix2dt
unix2d
```
