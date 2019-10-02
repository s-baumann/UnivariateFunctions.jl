using UnivariateFunctions: create_linear_interpolation, create_constant_interpolation_to_right, create_quadratic_spline
using UnivariateFunctions: create_constant_interpolation_to_left, years_from_global_base, evaluate
using Dates

tol = 10*eps()
const global_base_date = Date(2000,1,1)
StartDate = Date(2018, 7, 21)

x = StartDate .+ Dates.Day.(2 .* (1:1000 .- 1))

function ff(x::Date)
    days_between = years_from_global_base(x)
    return log(days_between) + sqrt(days_between)
end
y = ff.(x)

spline = create_linear_interpolation(x,y)
# Test if interpolating
all(abs.(evaluate.(Ref(spline), x) .- y) .< tol)

x_float = years_from_global_base.(x)
coefficient_in_first_interval = (y[2] - y[1])/(x_float[2] - x_float[1])
(coefficient_in_first_interval - spline.functions_[1].functions_[2].a_) < tol
other_coefficient = (y[11] - y[10])/(x_float[11] - x_float[10])
(other_coefficient - spline.functions_[10].functions_[2].a_) < tol

left_const = create_constant_interpolation_to_left(x,y)
right_const = create_constant_interpolation_to_right(x,y)
all(abs.(evaluate.(Ref(left_const), x_float .- tol) .- evaluate.(Ref(right_const), x_float)) .< tol)

# Making interpolation splies with date periods.
xx = Array{Month,1}([Month(1), Month(2), Month(3)])
xx2 = Array{DatePeriod,1}([Month(1), Month(2), Year(1)])
yy = [1,2,3]
abs(create_constant_interpolation_to_left(xx,yy)(Day(45)) - 2.0) < eps()
abs(evaluate(create_constant_interpolation_to_right(xx,yy), Day(15)) - 1.0) < eps()
abs(evaluate(create_linear_interpolation(xx,yy), Month(1)) - 1.0) < eps()
abs(evaluate(create_quadratic_spline(xx,yy), Month(2)) - 2.0) < eps()

abs(evaluate(create_constant_interpolation_to_left(xx2,yy), Day(45)) - 2.0) < eps()
abs(evaluate(create_constant_interpolation_to_right(xx2,yy), Day(15)) - 1.0) < eps()
abs(evaluate(create_linear_interpolation(xx2,yy), Month(1)) - 1.0) < eps()
abs(evaluate(create_quadratic_spline(xx2,yy), Month(2)) - 2.0) < eps()
