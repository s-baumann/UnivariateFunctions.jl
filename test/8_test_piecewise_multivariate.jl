using UnivariateFunctions
const tol = 1e-15
fun1 = PE_Function(1.0,1.0,1.0,1)
fun2 = PE_Function(1.2,1.0,1.0,1)
fun3 = PE_Function(1.5,1.0,2.0,3)
fun4 = PE_Function(3.0,1.0,3.0,2)
fun5 = PE_Function(0.0,2.0,1.0,1)
fun6 = PE_Function(2.0,0.0,1.0,1)
fun7 = PE_Function(3.0,0.0,2.0,4)
fun8 = PE_Function(3.0,0.0,2.0,2)

mfun1 = Multivariate_PE_Function(1.0, Dict(Symbol.(["x", "y", "z"]) .=> [fun1, fun2, fun3]))
abs(mfun1.multiplier_ - 1.2*1.5) < 1e-10
mfun2 = Multivariate_PE_Function(1.0, Dict(Symbol.(["x", "y"]) .=> [fun3, fun4]))
mfun3 = Multivariate_PE_Function(1.0, Dict(Symbol.(["z", "y"]) .=> [fun4, fun1]))
mfun4 = Multivariate_PE_Function(1.0, Dict(Symbol.(["x", "y"]) .=> [fun6, fun7]))
mfun5 = Multivariate_PE_Function(1.0, Dict(Symbol.(["x", "z"]) .=> [fun6, fun8]))

msum = mfun5 + mfun4

func_array = Array{Union{Multivariate_Sum_Of_Functions,Undefined_Function},3}(undef,3,2,2)
func_array[1,1,1] = Multivariate_Sum_Of_Functions(mfun1)
func_array[1,1,2] = Multivariate_Sum_Of_Functions(mfun2)
func_array[1,2,1] = mfun2 + mfun1
func_array[1,2,2] = mfun2 + mfun4
func_array[2,1,1] = Multivariate_Sum_Of_Functions(2*mfun1)
func_array[2,1,2] = mfun1 + mfun5
func_array[2,2,1] = mfun1 + mfun5
func_array[2,2,2] = mfun2 + mfun1 + mfun5
func_array[3,:,:] .= Undefined_Function()

pw_func = Multivariate_Piecewise_Function(func_array, Symbol.(["x", "y", "z"]), Dict(Symbol.(["x", "y", "z"]) .=> [[0.0, 0.5, 1.0], [-5.0, -1.0], [0.0, 5.0]]))

###############################################################################
coordinates = Dict(Symbol.(["w","x", "y", "z"]) .=> [1.4, 0.5, 2.0, 3.0])
function test_result(func, eval_to, len = 1)
    val_test = abs(evaluate(func, coordinates) - eval_to) < 1e-05
    if (!val_test)
        print("Failed Val Test")
    end
    len_test = length(func.functions_) == len
    return all([val_test,len_test])
end
ismissing(evaluate(pw_func, Dict(Symbol.(["w","x", "y", "z"]) .=> [-Inf, -Inf, -Inf, -Inf])))
test_result(pw_func, evaluate( mfun1 + mfun5, coordinates), 12)
deriv = Dict{Symbol,Int}(Symbol.(["x", "z"]) .=> [2,1])
abs(evaluate(derivative(pw_func, deriv), coordinates) - evaluate(derivative(mfun1 + mfun5, deriv), coordinates)) < tol

lower_lims = Dict{Symbol,Float64}(Symbol.(["x", "y", "z"]) .=> [0.25,-4.0,0.0])
upper_lims = Dict{Symbol,Float64}(Symbol.(["x", "y", "z"]) .=> [0.75,-1.0,3.0])
evaluate_integral(pw_func, upper_lims, lower_lims) - ( evaluate_integral(pw_func, upper_lims, lower_lims)   )

# test algebra between Multivariate_Piecewise_Functions and Multivariate_PE_Functions.
test_result(pw_func + mfun4, evaluate( mfun1 + mfun5, coordinates) + evaluate(mfun4, coordinates), 12)
test_result(pw_func - mfun4, evaluate( mfun1 + mfun5, coordinates) - evaluate(mfun4, coordinates), 12)
test_result(pw_func * mfun4, evaluate( mfun1 + mfun5, coordinates) * evaluate(mfun4, coordinates), 12)
test_result(mfun4 + pw_func, evaluate( mfun1 + mfun5, coordinates) + evaluate(mfun4, coordinates), 12)
test_result(mfun4 - pw_func,  evaluate(mfun4, coordinates) - evaluate( mfun1 + mfun5, coordinates), 12)
test_result(mfun4 * pw_func, evaluate( mfun1 + mfun5, coordinates) * evaluate(mfun4, coordinates), 12)
# Multivariate Sum of Functions
test_result(pw_func + msum, evaluate( mfun1 + mfun5, coordinates) + evaluate(msum, coordinates), 12)
test_result(pw_func - msum, evaluate( mfun1 + mfun5, coordinates) - evaluate(msum, coordinates), 12)
test_result(pw_func * msum, evaluate( mfun1 + mfun5, coordinates) * evaluate(msum, coordinates), 12)
test_result(msum + pw_func, evaluate( mfun1 + mfun5, coordinates) + evaluate(msum, coordinates), 12)
test_result(msum - pw_func,  evaluate(msum, coordinates) - evaluate( mfun1 + mfun5, coordinates), 12)
test_result(msum * pw_func, evaluate( mfun1 + mfun5, coordinates) * evaluate(msum, coordinates), 12)


# Algebra between different Multivariate_Piecewise_Functions
func_array = Array{Union{Undefined_Function,Multivariate_Sum_Of_Functions},3}(undef,3,2,2)
func_array[1,1,1] = Multivariate_Sum_Of_Functions(mfun5)
func_array[1,1,2] = Multivariate_Sum_Of_Functions(mfun3)
func_array[1,2,1] = mfun3 + mfun1
func_array[1,2,2] = Multivariate_Sum_Of_Functions(2*mfun1)
func_array[2,1,1] = mfun2 + mfun4
func_array[2,1,2] = mfun1 + mfun5
func_array[2,2,1] = mfun5 + mfun5
func_array[2,2,2] = mfun4 + mfun1 + mfun5
func_array[3,:,:] .= mfun4 + mfun5 + mfun1

pw_func2 = Multivariate_Piecewise_Function(func_array, Symbol.(["x", "y", "w"]), Dict(Symbol.(["x", "y", "w"]) .=> [[0.2, 0.5, 1.0], [-4.0, -1.0], [0.0, 5.0]]))
test_result(pw_func + pw_func2, evaluate( (mfun1 + mfun5) + (mfun5 + mfun5), coordinates), 3*4*2*2)
test_result(pw_func - pw_func2, evaluate( (mfun1 + mfun5) - (mfun5 + mfun5), coordinates), 3*4*2*2)
test_result(pw_func * pw_func2, evaluate( (mfun1 + mfun5) * (mfun5 + mfun5), coordinates), 3*4*2*2)
