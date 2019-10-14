using UnivariateFunctions: evaluate, years_between, years_from_global_base, PE_Function, Sum_Of_Functions, change_base_of_PE_Function, derivative, indefinite_integral, right_integral, left_integral, evaluate_integral
using Dates
tol = 10*eps()

today = Date(2000,1,1)
pe_func = PE_Function(1.0,2.0,today, 3)
(pe_func.base_ - years_from_global_base(today))   < tol
date_in_2020 = Date(2020,1,1)
pe_func2 = PE_Function(1.0,2.0,date_in_2020, 3)
(pe_func2.base_ - years_from_global_base(date_in_2020))   < tol
abs(evaluate(pe_func, date_in_2020) - evaluate(pe_func, years_from_global_base(date_in_2020)) ) < tol

#Sum of functions
sum_func = Sum_Of_Functions([pe_func, PE_Function(2.0,2.5,today, 3) ])
abs(evaluate(sum_func, date_in_2020) - evaluate(sum_func, years_from_global_base(date_in_2020)) ) < tol

# left and right integrals
l_int = left_integral(pe_func, today)
(evaluate(l_int, date_in_2020) - evaluate_integral(pe_func,today,date_in_2020)) < tol
r_int = right_integral(pe_func, date_in_2020)
(evaluate(r_int, today) - evaluate_integral(pe_func,today,date_in_2020)) < tol

# With DateTimes
today_time = DateTime(2000,1,1, 10, 10, 1)
pe_func = PE_Function(1.0,2.0,today_time, 3)


# left and right integrals
later_today = DateTime(2000,1,1, 10, 14, 1)
l_int = left_integral(pe_func, today)
(evaluate(l_int, later_today) - evaluate_integral(pe_func,today,later_today)) < tol
r_int = right_integral(pe_func, later_today)
(evaluate(r_int, today) - evaluate_integral(pe_func,today,later_today)) < tol
