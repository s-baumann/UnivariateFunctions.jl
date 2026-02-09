using Test

@testset "Basic Tests" begin
    using UnivariateFunctions
    tol = 10*eps()

    # Make undefined function
    und = Undefined_Function()
    @test ismissing(und(5))
    results = und.([5, 9]) # Test broadcasting
    @test ismissing(results[1])
    @test ismissing(results[2])

    # pe functions
    function test_pe(a::Float64,b::Float64,base::Float64,d::Int,x::Float64)
        return a * exp(b*(x-base)) * (x-base)^d
    end
    f1 = PE_Function(1.0, 2.0,3.0, 4)
    f1_test_result = test_pe(1.0, 2.0,3.0, 4,5.0)
    @test abs(evaluate(f1, 5.0) - f1_test_result) < tol
    f2 = PE_Function(0.0, 2.0,3.0, 4)
    f2_test_result = 0.0
    @test abs(evaluate(f2, 5.0) - f2_test_result) < tol
    f3 = PE_Function{Float64,Int64}(1.0, 8.0,7.0, 8)
    f3_test_result = test_pe(1.0, 8.0,7.0, 8,5.0)
    @test abs(f3(5.0) - f3_test_result) < tol
    f4 = PE_Function(1.0)
    f4_test_result = test_pe(1.0, 0.0,0.0, 0,5.0)
    @test abs(f4(5.0) - f4_test_result) < tol

    # Sum of functions
    sum0 =  Sum_Of_Functions([])
    @test typeof(sum0) == UnivariateFunctions.Sum_Of_Functions
    sum1 =  Sum_Of_Functions([f1])
    sum1 = Sum_Of_Functions(sum1)
    @test typeof(sum1) == UnivariateFunctions.Sum_Of_Functions
    sum3 =  Sum_Of_Functions([f2,f3])
    @test typeof(sum3) == UnivariateFunctions.Sum_Of_Functions
    sum4 =  Sum_Of_Functions([f1,f2,f3])
    @test typeof(sum4) == UnivariateFunctions.Sum_Of_Functions
    @test length(sum4.functions_) == 2
    sum5 =  Sum_Of_Functions([sum3, sum4])
    @test typeof(sum5) == UnivariateFunctions.Sum_Of_Functions
    @test length(sum5.functions_) == 2
    @test abs(sum1(5.0) - f1_test_result) < tol
    @test abs(evaluate(sum5, 5.0) - 2*f3_test_result - f1_test_result) < 100* tol
    @test abs(evaluate(sum0, 5.0)) < tol
    sum6 =  Sum_Of_Functions([sum3, sum4, und])
    @test typeof(sum6)  == Undefined_Function

    fl = 8.9
    integ = 7

    function test_result(func, expected_type, eval_to, len = 1)
        val_test = abs(evaluate(func, 5.0) - eval_to) < 1e-09
        if (!val_test)
            print("Failed Val Test")
        end
        type_test = isa(func, expected_type)
        if (!type_test)
            print("Failed Type Test")
        end
        if typeof(func) == UnivariateFunctions.Sum_Of_Functions
            len_test = length(func.functions_) == len
            if (!len_test)
                print("Failed Length Test")
            end
            return all([val_test,type_test,len_test])
        else
            return all([val_test,type_test])
        end
    end

    # pe and Float
    @test test_result(f1 + fl, UnivariateFunctions.Sum_Of_Functions, f1_test_result + fl,2 )
    @test test_result(f1 - fl, UnivariateFunctions.Sum_Of_Functions, f1_test_result - fl,2 )
    @test test_result(f1 * fl, UnivariateFunctions.PE_Function, (f1_test_result) * fl )
    @test test_result(f1 / fl, UnivariateFunctions.PE_Function, (f1_test_result) / fl )
    @test test_result(fl + f1, UnivariateFunctions.Sum_Of_Functions, f1_test_result + fl,2 )
    @test test_result(fl - f1, UnivariateFunctions.Sum_Of_Functions, fl - f1_test_result,2 )
    @test test_result(fl * f1, UnivariateFunctions.PE_Function, (f1_test_result) * fl )
    # pe and Int
    @test test_result(f1 + integ, UnivariateFunctions.Sum_Of_Functions, f1_test_result + integ,2 )
    @test test_result(f1 - integ, UnivariateFunctions.Sum_Of_Functions, f1_test_result - integ,2 )
    @test test_result(f1 * integ, UnivariateFunctions.PE_Function, (f1_test_result) * integ )
    @test test_result(f1 / integ, UnivariateFunctions.PE_Function, (f1_test_result) / integ )
    @test test_result(integ + f1, UnivariateFunctions.Sum_Of_Functions, f1_test_result + integ,2 )
    @test test_result(integ - f1, UnivariateFunctions.Sum_Of_Functions, -f1_test_result + integ,2 )
    @test test_result(integ * f1, UnivariateFunctions.PE_Function, (f1_test_result) * integ )


    # Sum and Float
    @test test_result(sum4 + fl, UnivariateFunctions.Sum_Of_Functions, f1_test_result + f3_test_result + fl, 3 )
    @test test_result(sum4 - fl, UnivariateFunctions.Sum_Of_Functions, f1_test_result + f3_test_result - fl, 3 )
    @test test_result(sum4 * fl, UnivariateFunctions.Sum_Of_Functions, (f1_test_result + f3_test_result) * fl, 2 )
    @test test_result(sum4 / fl, UnivariateFunctions.Sum_Of_Functions, (f1_test_result + f3_test_result) / fl, 2 )
    @test test_result(fl + sum4, UnivariateFunctions.Sum_Of_Functions, f1_test_result + f3_test_result + fl, 3 )
    @test test_result(fl - sum4, UnivariateFunctions.Sum_Of_Functions, -f1_test_result -f3_test_result + fl, 3 )
    @test test_result(fl * sum4, UnivariateFunctions.Sum_Of_Functions, (f1_test_result + f3_test_result) * fl, 2 )
    # Sum and Int
    @test test_result(sum4 + integ, UnivariateFunctions.Sum_Of_Functions, f1_test_result + f3_test_result + integ, 3 )
    @test test_result(sum4 - integ, UnivariateFunctions.Sum_Of_Functions, f1_test_result + f3_test_result - integ, 3 )
    @test test_result(sum4 * integ, UnivariateFunctions.Sum_Of_Functions, (f1_test_result + f3_test_result) * integ, 2 )
    @test test_result(sum4 / integ, UnivariateFunctions.Sum_Of_Functions, (f1_test_result + f3_test_result) / integ, 2 )
    @test test_result(integ + sum4, UnivariateFunctions.Sum_Of_Functions, f1_test_result + f3_test_result + integ, 3 )
    @test test_result(integ - sum4, UnivariateFunctions.Sum_Of_Functions, -f1_test_result -f3_test_result + integ, 3 )
    @test test_result(integ * sum4, UnivariateFunctions.Sum_Of_Functions, (f1_test_result + f3_test_result) * integ, 2 )

    ### Making sums with addition and subtraction.
    @test evaluate(sum4, 5.0) == evaluate(f1 + f2 + f3, 5.0)
    @test abs(evaluate(sum4, 5.0) - 2*evaluate(f3, 5.0) - evaluate(f1 + f2 - f3, 5.0)) < 0.001

    ### Cross-type subtraction tests (using values large enough to distinguish + from -)
    pe_a = PE_Function(10.0, 0.0, 0.0, 0)  # constant 10
    pe_b = PE_Function(3.0, 0.0, 0.0, 0)   # constant 3
    pe_c = PE_Function(2.0, 0.0, 0.0, 1)   # 2x
    sof_ab = Sum_Of_Functions([pe_a, pe_c])  # 10 + 2x
    # Sum_Of_Functions - PE_Function
    @test abs(evaluate(sof_ab - pe_b, 5.0) - 17.0) < tol  # (10 + 2*5) - 3 = 17
    # PE_Function - Sum_Of_Functions
    @test abs(evaluate(pe_b - sof_ab, 5.0) - (-17.0)) < tol  # 3 - (10 + 2*5) = -17
    # Sum_Of_Functions - Sum_Of_Functions
    sof_cd = Sum_Of_Functions([pe_b, pe_c])  # 3 + 2x
    @test abs(evaluate(sof_ab - sof_cd, 5.0) - 7.0) < tol  # (10+2*5) - (3+2*5) = 7

    ### Changing of base
    f1_unchanged = change_base_of_PE_Function(f1, 3.0)
    @test isa(f1_unchanged, UnivariateFunctions.PE_Function)
    @test abs(f1.base_ - f1_unchanged.base_) < tol

    f1_changed = change_base_of_PE_Function(f1, 4.0)
    @test typeof(f1_changed) == UnivariateFunctions.Sum_Of_Functions
    @test abs(f1_changed.functions_[1].base_ - 4.0) < tol
    @test abs(f1_changed.functions_[3].base_ - 4.0) < tol
    @test abs(evaluate(f1_changed, 5.0) - f1_test_result) < 100*tol # Changing bases should not change this.

    f1_changed_again = change_base_of_PE_Function(f1, -1.0)
    @test typeof(f1_changed_again) == UnivariateFunctions.Sum_Of_Functions
    @test abs(f1_changed_again.functions_[1].base_ + 1.0) < tol
    @test abs(f1_changed_again.functions_[3].base_ + 1.0) < tol
    @test abs(evaluate(f1_changed_again, 5.0) - f1_test_result) < 1000000*tol # Changing bases should not change this.

    f3_changed = change_base_of_PE_Function(f3, 12.0)
    @test typeof(f3_changed) == UnivariateFunctions.Sum_Of_Functions
    @test abs(f3_changed.functions_[1].base_ - 12.0) < tol
    @test abs(f3_changed.functions_[3].base_ - 12.0) < tol
    @test abs(evaluate(f3_changed, 5.0) - f3_test_result) < 100*tol # Changing bases should not change this.

    ### multiplication of functions
    @test test_result( f1 * f3 , UnivariateFunctions.Sum_Of_Functions, (f1_test_result * f3_test_result), 9 )
    @test test_result( f1 * PE_Function(1.0,2.0,3.0,4) , UnivariateFunctions.PE_Function, (f1_test_result * f1_test_result), 1 )
    @test test_result( f1 * sum4 , UnivariateFunctions.Sum_Of_Functions, f1_test_result * (f1_test_result + f3_test_result), 10 )
    @test test_result( sum4 * f1 , UnivariateFunctions.Sum_Of_Functions, f1_test_result * (f1_test_result + f3_test_result), 10 )
    @test test_result( sum4 * sum4 , UnivariateFunctions.Sum_Of_Functions, (f1_test_result + f3_test_result)^2, 11 )

    ### Powers pe
    @test test_result( f1 ^ 0 , UnivariateFunctions.PE_Function, 1.0 )
    @test test_result( f1 ^ 1 , UnivariateFunctions.PE_Function, (f1_test_result) )
    @test test_result( f1 ^ 2 , UnivariateFunctions.PE_Function, (f1_test_result * f1_test_result) )
    @test abs(evaluate(f1 ^ 4 ,5.0) - (f1_test_result * f1_test_result * f1_test_result * f1_test_result)) < 1e-03 # This fails normal tol.

    ### Powers sums
    @test test_result( sum4 ^ 0 , UnivariateFunctions.PE_Function, 1.0 )
    @test test_result( sum4 ^ 1 , UnivariateFunctions.Sum_Of_Functions, (f1_test_result + f3_test_result), 2 )
    @test test_result( sum4 ^ 2 , UnivariateFunctions.Sum_Of_Functions, (f1_test_result + f3_test_result)^2, 11 )
    # Higher powers give an underflow problem with base changes.


    ## Testing Calculus.
    pe_const    = PE_Function(2.0,0.0,5.0,0)
    pe_lin      = PE_Function(2.0,0.0,5.0,1)
    pe_quad     = PE_Function(2.0,0.0,2.0,2)
    pe_exp      = PE_Function(2.0,2.0,1.0,0)
    pe_exp_quad = PE_Function(2.0,2.0,2.0,2)

    # Linear gradient constant
    @test abs(evaluate(derivative(pe_lin),5.0) - evaluate(derivative(pe_lin),1.0) ) < tol
    @test typeof(derivative(pe_lin)) <: UnivariateFunctions.PE_Function
    # This is also linear
    @test abs(evaluate(derivative(derivative(pe_quad)),5.0) -evaluate(derivative(derivative(pe_quad)),9.0) ) < tol
    # derivative into a sum of functions
    @test typeof(derivative(pe_exp_quad)) == UnivariateFunctions.Sum_Of_Functions
    @test abs(evaluate(derivative(pe_exp_quad),5.0) - ( 2.0*exp(2.0*(5.0-2.0))*(5.0-2.0)*(2 + 2.0*(5.0-2.0)) ) ) < tol
    # Derivaitve of sum
    @test abs(evaluate(derivative(sum4),5.0) - evaluate(derivative(f1),5.0) - evaluate(derivative(f3),5.0) ) < 1e-09

    # integral of constant
    @test indefinite_integral(pe_const).d_ == 1
    # Integral of quadratic
    @test abs(evaluate_integral(pe_quad,1.0,2.0) - ((2/3)*2^3 - 4*2^2 + 8*2) + ((2/3)*1^3 - 4*1^2 + 8*1)) < tol
    # Integral of exponential
    @test abs(evaluate_integral(pe_exp,0.2,3.8) - (exp(2.0*(3.8-1.0)) - exp(2.0*(0.2-1.0))  )) < tol
    # Integral of combined.
    combined_analytical_integral = PE_Function(1.0,2.0,2.0,2) - PE_Function(1.0,2.0,2.0,1) + PE_Function(0.5,2.0,2.0,0)
    @test abs(evaluate_integral(pe_exp_quad,0.2,3.8) - ( evaluate(combined_analytical_integral, 3.8) - evaluate(combined_analytical_integral, 0.2))) < tol

    # left and right integrals
    l_int = left_integral(pe_exp_quad, 0.2)
    @test (evaluate(l_int, 3.8) - evaluate_integral(pe_exp_quad,0.2,3.8)) < tol
    r_int = right_integral(pe_exp_quad, 3.8)
    @test (evaluate(r_int, 0.2) - evaluate_integral(pe_exp_quad,0.2,3.8)) < tol

    ## Testing integral guards for piecewise functions with out-of-domain limits
    pw = Piecewise_Function([1.0, 5.0], [pe_const, pe_lin])
    # right_integral before domain start returns Undefined_Function
    @test isa(right_integral(pw, 0.0), Undefined_Function)
    # left_integral before domain start returns Undefined_Function
    @test isa(left_integral(pw, 0.0), Undefined_Function)
    # evaluate_integral before domain start returns missing
    @test ismissing(evaluate_integral(pw, 0.0, 0.5))
    @test ismissing(evaluate_integral(pw, -1.0, 2.0))
    # evaluate_integral within domain works normally
    @test !ismissing(evaluate_integral(pw, 1.0, 5.0))
    # right_integral and left_integral within domain work normally
    @test isa(right_integral(pw, 1.0), Piecewise_Function)
    @test isa(left_integral(pw, 5.0), Piecewise_Function)

    ## Testing higher powers on Sum_Of_Functions
    sum_for_pow = PE_Function(2.0, 0.0, 0.0, 0) + PE_Function(1.0, 0.0, 0.0, 1)  # 2 + x
    pow3 = sum_for_pow ^ 3  # (2+x)^3
    @test abs(evaluate(pow3, 1.0) - 27.0) < 1e-06  # (2+1)^3 = 27
    @test abs(evaluate(pow3, 3.0) - 125.0) < 1e-06  # (2+3)^3 = 125
    pow5 = sum_for_pow ^ 5  # (2+x)^5
    @test abs(evaluate(pow5, 1.0) - 243.0) < 1e-03  # (2+1)^5 = 243
    @test abs(evaluate(pow5, 2.0) - 1024.0) < 1e-03  # (2+2)^5 = 1024

    ## Testing linear rescaling
    # For PE_Functions
    test_func = PE_Function(0.5,1.2,0.9,2)
    inp = 4.0
    alpha = 0.435
    beta = -1.52
    rescaled_input = alpha*inp + beta
    converted_test_func = convert_to_linearly_rescale_inputs(test_func, alpha, beta)
    @test abs(evaluate(test_func, inp) - evaluate(converted_test_func, rescaled_input)) < 1e-10
    # For undefined.
    @test typeof(convert_to_linearly_rescale_inputs(Undefined_Function(), alpha, beta)) == UnivariateFunctions.Undefined_Function
    # SumOfFunctions
    sumFunc = test_func + pe_exp_quad + pe_exp
    converted_sumFunc = convert_to_linearly_rescale_inputs(sumFunc, alpha, beta)
    @test abs(evaluate(sumFunc, inp) - evaluate(converted_sumFunc, rescaled_input)) < 1e-10



    ### Test undefined
    @test typeof(und + und)  == Undefined_Function
    @test typeof(und - und)  == Undefined_Function
    @test typeof(und / und)  == Undefined_Function
    @test typeof(und * und)  == Undefined_Function

    # PE_Function
    @test typeof(f1 + und)  == Undefined_Function
    @test typeof(f1 - und)  == Undefined_Function
    @test typeof(f1 / und)  == Undefined_Function
    @test typeof(f1 * und)  == Undefined_Function
    @test typeof(und + f1)  == Undefined_Function
    @test typeof(und - f1)  == Undefined_Function
    @test typeof(und / f1)  == Undefined_Function
    @test typeof(und * f1)  == Undefined_Function

    # Sum_Of_Functions
    @test typeof(sum0 + und)  == Undefined_Function
    @test typeof(sum0 - und)  == Undefined_Function
    @test typeof(sum0 / und)  == Undefined_Function
    @test typeof(sum0 * und)  == Undefined_Function
    @test typeof(und + sum0)  == Undefined_Function
    @test typeof(und - sum0)  == Undefined_Function
    @test typeof(und / sum0)  == Undefined_Function
    @test typeof(und * sum0)  == Undefined_Function
    @test typeof(und + und)  == Undefined_Function
    @test typeof(und - und)  == Undefined_Function
    @test typeof(und / und)  == Undefined_Function
    @test typeof(und * und)  == Undefined_Function

    # Numbers
    @test typeof(und + 4)  == Undefined_Function
    @test typeof(und - 4)  == Undefined_Function
    @test typeof(und / 4)  == Undefined_Function
    @test typeof(und * 4)  == Undefined_Function
    @test typeof(und ^ 4)  == Undefined_Function
    @test typeof(4 + und)  == Undefined_Function
    @test typeof(4 - und)  == Undefined_Function
    @test typeof(4 * und)  == Undefined_Function

    @test typeof(und + 4.0)  == Undefined_Function
    @test typeof(und - 4.0)  == Undefined_Function
    @test typeof(und / 4.0)  == Undefined_Function
    @test typeof(und * 4.0)  == Undefined_Function
    @test typeof(4.0 + und)  == Undefined_Function
    @test typeof(4.0 - und)  == Undefined_Function
    @test typeof(4.0 * und)  == Undefined_Function

end
