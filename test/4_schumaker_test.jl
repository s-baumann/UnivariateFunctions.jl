using Test


@testset "Schumaker Tests" begin
    using UnivariateFunctions
    using Dates

    tol = 10*eps()
    global_base_date = Date(2000,1,1)
    StartDate = Date(2018, 7, 21)

    x = Array{Date}(undef, 1000)
    for i in 1:1000
        x[i] = StartDate +Dates.Day(2* (i-1))
    end

    function ff(x::Date)
        days_between = years_from_global_base_date(x)
        return log(days_between) + sqrt(days_between)
    end
    y = ff.(x)

    spline = create_quadratic_spline(x,y)
    # Test if interpolating
    @test all(abs.(spline.(x) .- y) .< 0.001)

    # Testing third derivatives
    third_derivatives = evaluate.(Ref(derivative(derivative(derivative(spline)))), x)
    @test maximum(third_derivatives) < tol

    # Testing Integrals
    function analytic_integral(lhs,rhs)
        lhs_in_days = years_from_global_base_date(lhs)
        rhs_in_days = years_from_global_base_date(rhs)
        return rhs_in_days*log(rhs_in_days) - rhs_in_days + (2/3)*rhs_in_days^(3/2) - (lhs_in_days*log(lhs_in_days) - lhs_in_days)- (2/3)*lhs_in_days^(3/2)
    end

    lhs = StartDate
    rhs = StartDate + Dates.Month(16)
    numerical_integral  = evaluate_integral(spline, lhs,rhs)
    numerical_integral2 = evaluate(right_integral(spline, lhs), rhs)
    numerical_integral3 = evaluate(left_integral(spline, rhs), lhs)

    analytical = analytic_integral(lhs,rhs)
    @test abs(  analytical - numerical_integral  ) < 0.0001
    @test abs(  analytical - numerical_integral2  ) < 0.0001
    @test abs(  analytical - numerical_integral3  ) < 0.0001
end






