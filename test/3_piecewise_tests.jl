using Test

@testset "Piecewise Tests" begin

    using UnivariateFunctions
    using Dates

    tol = 10*eps()

    # Testing constructors.
    und    = Undefined_Function()
    before = PE_Function(0.0,0.0,0.0,1)
    first_ = PE_Function(1.0,0.0,0.0,1)
    second = PE_Function(20.0,0.0,0.0,2)
    last_  = Sum_Of_Functions([first_, second])

    f0 = Piecewise_Function([-1.0,3.0, 10.0], [before, first_, second])
    @test !ismissing(f0(-0.75))
    f0 = trim_piecewise_function(f0, -0.5, 9.0)
    @test ismissing(f0(-0.75))
    f0 = Piecewise_Function([-Inf,3.0, 10.0], [before, first_, second])
    @test !ismissing(f0(-0.75))

    f1 = Piecewise_Function([-Inf, -1.0,3.0, 10.0], [before, first_, second, last_])
    @test abs(evaluate(f1, -10.0) - 0.0 ) < tol
    @test abs(evaluate(f1,   1.0) - 1.0 ) < tol
    @test abs(evaluate(f1,   5.0) - 500.0) < tol
    @test abs(evaluate(f1,  20.0) - 20.0 - 20*20*20) < tol
    f2 = Piecewise_Function([-0.1, 0.0,2.0, 40.0], [before, first_, second, last_])

    f3 =  Piecewise_Function([-0.1, 0.0,2.0, 40.0], [before, f1, second, last_])
    f3.starts_ == [-Inf, -0.100, 0.00, 2.00, 40.0]
    @test isa(f3.functions_[1],UnivariateFunctions.Undefined_Function)
    @test isa(f3.functions_[2],UnivariateFunctions.PE_Function)
    @test isa(f3.functions_[3],UnivariateFunctions.PE_Function)
    @test isa(f3.functions_[4],UnivariateFunctions.PE_Function)
    @test isa(f3.functions_[5],UnivariateFunctions.Sum_Of_Functions)

    f4 =  Piecewise_Function([-2.1, -1.1,4.0, 40.0], [second, f1, second, last_])
    f4.starts_ == [-Inf, -2.10, -1.10, -1.00, 3.00, 4.00, 40.0]
    @test isa(f4.functions_[1], UnivariateFunctions.Undefined_Function)
    @test isa(f4.functions_[2], UnivariateFunctions.PE_Function)
    @test isa(f4.functions_[3], UnivariateFunctions.PE_Function)
    @test isa(f4.functions_[4], UnivariateFunctions.PE_Function)
    @test isa(f4.functions_[5], UnivariateFunctions.PE_Function)
    @test isa(f4.functions_[6], UnivariateFunctions.PE_Function)
    @test isa(f4.functions_[7], UnivariateFunctions.Sum_Of_Functions)

    # Testing the other Piecewise_Function constructors
    f5 = Piecewise_Function(Date(2020,05,01), second)
    @test ismissing(f5(Date(2019,05,01)))
    f6 = Piecewise_Function([Date(2020,05,01)], [second])
    @test ismissing(f6(Date(2019,05,01)))
    f7 = Piecewise_Function([77.0, Inf], [second, second])
    @test ismissing(f7(Date(2019,05,01)))

    function test_result(func, eval0, eval5, len = 1)
        val_test0 = abs(evaluate(func, 0.0) - eval0) < 1e-09
        if (!val_test0)
            print("Failed Val Test at 0.0")
        end
        val_test5 = abs(evaluate(func, 5.0) - eval5) < 1e-09
        if (!val_test5)
            print("Failed Val Test at 5.0")
        end
        if isa(func, UnivariateFunctions.Piecewise_Function)
            len_test = length(func.functions_) == len
            if (!len_test)
                print("Failed Length Test")
            end
            return all([val_test0,val_test5,len_test])
        else
            return all([val_test0,val_test5])
        end
    end

    # Testing with Ints
    @test test_result(f1 + 5, 5.0, 505.0, 4)
    @test test_result(5 + f1, 5.0, 505.0, 4)
    @test test_result(f1 - 5, -5.0, 495.0, 4)
    @test test_result(5 - f1, 5.0, -495.0, 4)
    @test test_result(f1 * 5, 0.0, 2500.0, 4)
    @test test_result(5 * f1, 0.0, 2500.0, 4)
    @test test_result(f1 / 5, 0.0, 100.0, 4)
    # And with Float64s
    @test test_result(f1 + 5.0, 5.0, 505.0, 4)
    @test test_result(5.0 + f1, 5.0, 505.0, 4)
    @test test_result(f1 - 5.0, -5.0, 495.0, 4)
    @test test_result(5.0 - f1, 5.0, -495.0, 4)
    @test test_result(f1 * 5.0, 0.0, 2500.0, 4)
    @test test_result(5.0 * f1, 0.0, 2500.0, 4)
    @test test_result(f1 / 5.0, 0.0, 100.0, 4)

    # Testing with other functions
    f5 = f4 + f1
    f5.starts_[2] == -2.1
    @test typeof(f5.functions_[1]) == UnivariateFunctions.Undefined_Function

    @test typeof(f4 + und) == UnivariateFunctions.Undefined_Function
    @test typeof(und + f4) == UnivariateFunctions.Undefined_Function
    @test typeof(f4 - und) == UnivariateFunctions.Undefined_Function
    @test typeof(und - f4) == UnivariateFunctions.Undefined_Function
    @test typeof(f4 * und) == UnivariateFunctions.Undefined_Function
    @test typeof(und * f4) == UnivariateFunctions.Undefined_Function

    @test typeof(-1*f4) == UnivariateFunctions.Piecewise_Function

    @test test_result(f1 + first_, 0.0, 505.0, 4)
    @test test_result(first_ + f1, 0.0, 505.0, 4)
    @test test_result(f1 - first_, 0.0, 495.0, 4)
    @test test_result(first_ - f1, 0.0, -495.0, 4)
    @test test_result(f1 * first_, 0.0, 2500.0, 4)
    @test test_result(first_ * f1, 0.0, 2500.0, 4)

    @test test_result(f1 + last_, 0.0, 1005.0, 4)
    @test test_result(last_ + f1, 0.0, 1005.0, 4)
    @test test_result(f1 - last_, 0.0, -5.0, 4)
    @test test_result(last_ - f1, 0.0, 5.0, 4)
    @test test_result(f1 * last_, 0.0, 252500.0, 4)
    @test test_result(last_ * f1, 0.0, 252500.0, 4)

    @test test_result(f1 + f4, 0.0, 1000.0, 8)
    @test test_result(f1 - f4, 0.0, 0.0, 8)
    @test test_result(f1 * f4, 0.0, 250000.0, 8)

    # Testing linear rescaling
    #test function will be f4
    alpha = 0.95
    beta = 3.0
    X = [-0.3, 0.0, 1.3, 1.5, 4.0, 7.0, 20.0, 40.0, 70.0, 100.0]
    X_converted = (alpha .* X) .+ beta
    f4_converted = convert_to_linearly_rescale_inputs(f4, alpha, beta)
    y = evaluate.(f4,X)
    converted_y = evaluate.(f4_converted, X_converted)
    @test sum(abs.(y .- converted_y)) < 1e-10

    # Testing with Undefined
    und = Undefined_Function()
    @test typeof(last_ + und) == Undefined_Function
    @test typeof(last_ - und) == Undefined_Function
    @test typeof(last_ / und) == Undefined_Function
    @test typeof(last_ * und) == Undefined_Function
    @test typeof(und + last_) == Undefined_Function
    @test typeof(und - last_) == Undefined_Function
    @test typeof(und / last_) == Undefined_Function
    @test typeof(und * last_) == Undefined_Function
    @test typeof(und / f5) == Undefined_Function
    @test typeof(f5 / und) == Undefined_Function
end
