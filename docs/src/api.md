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

### Interpolation

```@docs
create_quadratic_spline
create_constant_interpolation_to_right
create_constant_interpolation_to_left
create_linear_interpolation
```

### Approximation

```@docs
create_ols_approximation
create_chebyshev_approximation
```

### Internal Functions

```@docs
change_base_of_PE_Function
trim_piecewise_function
convert_to_linearly_rescale_inputs
get_chevyshevs_up_to
```

### Date Conversions

```@docs
years_between
years_from_global_base
period_length
```
