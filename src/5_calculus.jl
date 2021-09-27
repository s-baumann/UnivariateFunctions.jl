import SchumakerSpline.evaluate
import SchumakerSpline.evaluate_integral

"""
    right_integral(f::UnivariateFunction, left::Real)
    right_integral(f::UnivariateFunction, left::Union{Date,DateTime})

This calculates the integral of a function from a left limit and returns it as a UnivariateFunction.
So if you were to then evaluate this integral function at a point x then you would
get the integral between left and x.

If a `Date`, `DateTime` is input then it is first converted to a scalar with the `years_from_global_base` function.

### Inputs
* `f` - A `UnivariateFunction`.
* `left` - A left limit (scalar)
### Returns
* A `UnivariateFunction`.
"""
function right_integral(f::UnivariateFunction, left::Real)
    indef_int = indefinite_integral(f)
    left_constant = evaluate(indef_int, left)
    return indef_int - left_constant
end

function right_integral(f::UnivariateFunction, left::Union{Date,DateTime})
    left_float = years_from_global_base(left)
    return right_integral(f, left_float)
end

function right_integral(f::Piecewise_Function, left::R) where R<:Real
    whole_number_of_intervals = length(f.starts_)
    which_interval_contains_left = searchsortedlast(f.starts_, left)
    number_of_intervals = whole_number_of_intervals - which_interval_contains_left + 1
    starts_    =  Array{R,1}(undef, number_of_intervals)
    starts_[1] = left
    first_indefinite_integral = indefinite_integral(f.functions_[which_interval_contains_left])
    functions_    =  Array{UnivariateFunction}(undef, number_of_intervals)
    if typeof(first_indefinite_integral) == UnivariateFunctions.Undefined_Function
        return Undefined_Function()
    end
    functions_[1] = first_indefinite_integral - evaluate(first_indefinite_integral, left)
    if number_of_intervals == 1
        return Piecewise_Function(starts_, functions_)
    end
    displacement  = which_interval_contains_left - 1
    for i in 2:number_of_intervals
        starts_[i]    = f.starts_[i + displacement]
        indef_int = indefinite_integral(f.functions_[i + displacement])
        if typeof(indef_int) == UnivariateFunctions.Undefined_Function
            return Piecewise_Function(starts_[1:i], vcat(functions_[1:(i-1)], Undefined_Function()  ) )
        end
        value_at_previous_end = evaluate(functions_[i-1], starts_[i])
        functions_[i] = indef_int - evaluate(indef_int, starts_[i]) + value_at_previous_end
    end
    return Piecewise_Function(starts_, functions_)
end

"""
    left_integral(f::UnivariateFunction, right::Real)
    left_integral(f::UnivariateFunction, right::Union{Date,DateTime})

This calculates the integral of a function from a right limit and returns it as a UnivariateFunction.
So if you were to then evaluate this integral function at a point x then you would
get the integral between right and x.

If a `Date`, `DateTime` is input then it is first converted to a scalar with the `years_from_global_base` function.

### Inputs
* `f` - A `UnivariateFunction`.
* `right` - A right limit (scalar)
### Returns
* A `UnivariateFunction`.
"""
function left_integral(f::UnivariateFunction, right::Real)
    indef_int = indefinite_integral(f)
    right_constant = evaluate(indef_int, right)
    return right_constant - indef_int
end

function left_integral(f::UnivariateFunction, right::Union{Date,DateTime})
    right_float = years_from_global_base(right)
    return left_integral(f, right_float)
end

function left_integral(f::Piecewise_Function, right::R) where R<:Real
    whole_number_of_intervals = length(f.starts_)
    which_interval_contains_right = searchsortedlast(f.starts_, right)
    number_of_intervals = which_interval_contains_right
    starts_    =  Array{R,1}(undef, number_of_intervals)
    starts_[number_of_intervals] = f.starts_[which_interval_contains_right]
    last_indefinite_integral = indefinite_integral(f.functions_[which_interval_contains_right])
    functions_    =  Array{UnivariateFunction}(undef, number_of_intervals)
    if typeof(last_indefinite_integral) == UnivariateFunctions.Undefined_Function
        return Undefined_Function()
    end
    functions_[number_of_intervals] = evaluate(last_indefinite_integral, right) - last_indefinite_integral
    if number_of_intervals == 1
        return Piecewise_Function(starts_, functions_)
    end
    for i in reverse(1:(which_interval_contains_right-1))
        starts_[i]    = f.starts_[i]
        value_at_next_start = evaluate(functions_[i+1], starts_[i+1])
        indef_int = indefinite_integral(f.functions_[i])
        if typeof(indef_int) == UnivariateFunctions.Undefined_Function
            return Piecewise_Function(starts_[(i+1):which_interval_contains_right], functions_[(i+1):which_interval_contains_right]   )
        end
        functions_[i] = evaluate(indef_int, starts_[i+1] ) - indef_int + value_at_next_start
    end
    return Piecewise_Function(starts_, functions_)
end

"""
    evaluate_integral(f::UnivariateFunction,left::Real, right::Real)
    evaluate_integral(f::UnivariateFunction,left::Union{Date,DateTime}, right::Union{Date,DateTime})

This calculates the integral of a function from a left limit to a right limit and returns a scalar.

If a `Date`, `DateTime` is input then it is first converted to a scalar with the `years_from_global_base` function.

### Inputs
* `f` - A `UnivariateFunction`.
* `left` - A left limit (scalar)
* `right` - A right limit (scalar)
### Returns
* A `Real`.
"""
function evaluate_integral(f::UnivariateFunction,left::Real, right::Real)
    indef_int  = indefinite_integral(f)
    left_eval  = evaluate(indef_int, left)
    right_eval = evaluate(indef_int, right)
    return (right_eval - left_eval)
end
function evaluate_integral(f::UnivariateFunction,left::Union{Date,DateTime}, right::Union{Date,DateTime})
    left_as_float  = years_from_global_base(left)
    right_as_float = years_from_global_base(right)
    return evaluate_integral(f, left_as_float, right_as_float)
end

function evaluate_integral(f::Piecewise_Function,left::Real, right::Real)
    return evaluate(right_integral(f, left), right)
end
