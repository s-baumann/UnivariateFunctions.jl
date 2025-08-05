using Test
@testset "Date Tests" begin
    using UnivariateFunctions
    using Dates
    tol = 10*eps()

    today = Date(2000,1,1)
    pe_func = PE_Function(1.0,2.0,today, 3)
    @test abs(pe_func.base_ - years_from_global_base_date(today))   < tol
    date_in_2020 = Date(2000,2,3)
    pe_func2 = PE_Function(1.0,0.00000000002,date_in_2020, 2)
    @test abs(pe_func2.base_ - years_from_global_base_date(date_in_2020))   < tol
    @test abs(evaluate(pe_func, date_in_2020) - evaluate(pe_func, years_from_global_base_date(date_in_2020)) ) < tol

    #Sum of functions
    sum_func = Sum_Of_Functions([pe_func, PE_Function(2.0,0.00000000000025,today, 1) ])
    @test abs(evaluate(sum_func, date_in_2020) - evaluate(sum_func, years_from_global_base_date(date_in_2020)) ) < tol

    # left and right integrals
    l_int = left_integral(pe_func, today)
    @test (evaluate(l_int, date_in_2020) - evaluate_integral(pe_func,today,date_in_2020)) < tol
    r_int = right_integral(pe_func, date_in_2020)
    @test (evaluate(r_int, today) - evaluate_integral(pe_func,today,date_in_2020)) < tol

    # With DateTimes
    today_time = DateTime(2000,1,1, 10, 10, 1)
    pe_func = PE_Function(1.0,0.00000000000020,today_time, 1)


    # left and right integrals
    later_today = DateTime(2000,1,1, 10, 14, 1)
    l_int = left_integral(pe_func, today)
    @test (evaluate(l_int, later_today) - evaluate_integral(pe_func,today,later_today)) < tol
    r_int = right_integral(pe_func, later_today)
    @test (evaluate(r_int, today) - evaluate_integral(pe_func,today,later_today)) < tol
end