using UnivariateFunctions
using DataFrames
using Dates
using Random
using Distributions
const tol = 1e-14
fun1 = PE_Function(1.0,1.0,1.0,1)
fun2 = PE_Function(1.2,1.0,1.0,1)
fun3 = PE_Function(1.5,1.0,2.0,3)
fun4 = PE_Function(3.0,1.0,3.0,2)
fun5 = PE_Function(0.0,2.0,1.0,1)
fun6 = PE_Function(2.0,0.0,1.0,1)
fun7 = PE_Function(3.0,0.0,2.0,4)
fun8 = PE_Function(3.0,0.0,2.0,2)
fun9 = PE_Function(3.0,0.0,Date(2015,1,1),4)
fun10 = PE_Function(3.0,0.0,Date(2016,1,1),2)

mfun1 = Multivariate_PE_Function(1.0, Dict(Symbol.(["x", "y", "z"]) .=> [fun1, fun2, fun3]))
abs(mfun1.multiplier_ - 1.2*1.5) < 1e-10
mfun2 = Multivariate_PE_Function(1.0, Dict(Symbol.(["x", "y"]) .=> [fun3, fun4]))
mfun3 = Multivariate_PE_Function(1.0, Dict(Symbol.(["z", "y"]) .=> [fun4, fun1]))
mfun4 = Multivariate_PE_Function(1.0, Dict(Symbol.(["x", "y"]) .=> [fun6, fun7]))
mfun5 = Multivariate_PE_Function(1.0, Dict(Symbol.(["x", "z"]) .=> [fun6, fun8]))
mfun6 = Multivariate_PE_Function(1.0, Dict(Symbol.(["x", "y", "z"]) .=> [fun9, fun6, fun10]))
# Testing if zero length if multiplier is zero.
length(Multivariate_PE_Function(1.0, Dict(Symbol.(["x", "y"]) .=> [fun4, fun5])).functions_) == 0
mSumFunction1 = Multivariate_Sum_Of_Functions([mfun1, mfun2])
mSumFunction2 = Multivariate_Sum_Of_Functions([mSumFunction1, mfun3])
mSumFunction3 = Multivariate_Sum_Of_Functions([mSumFunction1, mSumFunction2])

coordinates = Dict(Symbol.(["x", "y", "z"]) .=> [0.5, 2.0, 3.0])

# Conversion
typeof(convert(UnivariateFunctions.Multivariate_Sum_Of_Functions, mfun1)) == UnivariateFunctions.Multivariate_Sum_Of_Functions

function test_result(func, eval_to, len = 1)
    val_test = abs(evaluate(func, coordinates) - eval_to) < 1e-05
    if (!val_test)
        print("Failed Val Test")
    end
    len_test = length(func.functions_) == len
    return all([val_test,len_test])
end



# Testing numerical results
mfun1_value = evaluate(fun1, coordinates[:x]) * evaluate(fun2, coordinates[:y]) * evaluate(fun3, coordinates[:z])
test_result(mfun1, mfun1_value ,3)
mfun2_value = evaluate(fun3, coordinates[:x]) * evaluate(fun4, coordinates[:y])
test_result(mfun2, mfun2_value ,2)
test_result(mSumFunction1, mfun1_value + mfun2_value ,2)
mfun3_value = evaluate(fun4, coordinates[:z]) * evaluate(fun1, coordinates[:y])
test_result(mSumFunction3, 2* mfun1_value + 2*mfun2_value + mfun3_value ,5)

# Additions and subtractions
test_result( ((7 - (5.0 - mfun1))  + 5.0) + 7        , 14 + mfun1_value                ,2)
test_result( ((7 - (5.0 - mSumFunction1))  + 5.0) + 7, 14 +  mfun1_value + mfun2_value ,3)
# multiplications
test_result( (5 * mfun1) * 7.0        , 35* mfun1_value                  ,3)
test_result( (5 * mSumFunction1) * 7.0, 35*  (mfun1_value + mfun2_value) ,2)
# Divisions
test_result( (mfun1/2) / 3.0       ,  mfun1_value/6                  ,3)
test_result( (mSumFunction1/2) /3.0, (mfun1_value + mfun2_value)/6 ,2)
# powers
test_result( mfun1^0, 1,0)
test_result( mSumFunction1^0, 1,0)
test_result( mfun1^1, mfun1_value,3)
test_result( mSumFunction1^1, (mfun1_value + mfun2_value),2)
test_result( mfun1^2, mfun1_value^2,3)
test_result( mSumFunction1^2, (mfun1_value + mfun2_value)^2,26)
# multiplying functions
test_result( mfun1 * mfun2, mfun1_value * mfun2_value, 12)
# Using Dates
dateCoordinates = Dict(Symbol.(["x", "y", "z"]) .=> [Date(2017,1,1), 2, 17.0])
dateCoordinates2 = Dict(Symbol.(["x", "y", "z"]) .=> [years_from_global_base(dateCoordinates[:x]), 2.0, 17.0])
abs(evaluate(mfun6, dateCoordinates) - evaluate(mfun6, dateCoordinates2) ) < 1e-14

# Derivatives
test_result(derivative(mfun4, Symbol("x"), 2), 0.0, 0)
test_result(derivative(mfun4, Symbol("y"), 4), 24 * 6 * (-0.5) , 1)
test_result(derivative(mfun5, Dict(Symbol.(["x", "z"]) .=> [1,1])), 12 * (3-2)^1 , 1)
# Integration
lowers = Dict{Symbol,Float64}(Symbol.(["x", "y"]) .=> [2.0, 3.0])
uppers = Dict{Symbol,Float64}(Symbol.(["x", "y"]) .=> [2.5, 3.5])
abs(evaluate_integral(mfun4, uppers, lowers) - 0.6 * ((uppers[:x] - 1)^2 - (lowers[:x]-1)^2)*((uppers[:y]-2)^5 - (lowers[:y]-2)^5)) < tol
abs(evaluate_integral(Multivariate_Sum_Of_Functions([mfun4, mfun4]), uppers, lowers) - 1.2 * ((uppers[:x] - 1)^2 - (lowers[:x]-1)^2)*((uppers[:y]-2)^5 - (lowers[:y]-2)^5)) < tol


## Dataframe evaluations
Random.seed!(1992)
nObs = 100
dd = DataFrame()
dd[:x] = rand( Normal(),nObs) + 0.1 .* rand( Normal(),nObs)
dd[:z] = rand( Normal(),nObs) + 0.1 .* rand( Normal(),nObs)
dd[:y] = (dd[:x] ) .* (dd[:z]) .+ dd[:x] + rand( Normal(),nObs)
sum(abs.(evaluate(mfun1, dd) + evaluate(mfun2, dd) - evaluate(mSumFunction1, dd))) < tol
