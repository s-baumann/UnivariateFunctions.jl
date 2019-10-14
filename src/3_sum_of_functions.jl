import Base.+, Base.-, Base./, Base.*
import SchumakerSpline.evaluate
function evaluate(sf::Sum_Of_Functions, point::Real)
    total = 0.0
    for func in sf.functions_
        total = total + evaluate(func,point)
    end
    return total
end
function (s::Sum_Of_Functions)(x::Union{Real,Date,DateTime,DatePeriod})
    return evaluate(s, x)
end


function derivative(f::Sum_Of_Functions)
    return Sum_Of_Functions( derivative.(f.functions_))
end

function indefinite_integral(f::Sum_Of_Functions)
    return Sum_Of_Functions(indefinite_integral.(f.functions_))
end

function +(f::Sum_Of_Functions,number::Real)
    constant_function = PE_Function(number, 0.0,0.0,0)
    return Sum_Of_Functions(vcat(f.functions_, [constant_function]))
end
function -(f::Sum_Of_Functions, number::Real)
    return +(f, -number)
end
function *(f::Sum_Of_Functions, number::Real)
    funcs = deepcopy(f.functions_)
    for i in 1:length(funcs)
        funcs[i] = funcs[i] * number
    end
    return Sum_Of_Functions(funcs)
end
function /(f::Sum_Of_Functions, number::Real)
    return *(f, 1/number )
end

function +(f1::Sum_Of_Functions, f2::Sum_Of_Functions)
    return Sum_Of_Functions([f1,f2])
end
function +(f1::Sum_Of_Functions, f2::Piecewise_Function)
    added_functions = (f1 .+ f2.functions_)::Array
    return Piecewise_Function(f2.starts_,added_functions)
end

function +(f1::Piecewise_Function, f2::Sum_Of_Functions)
    return +(f2,f1)
end

function -(f1::Sum_Of_Functions, f2::Sum_Of_Functions)
    return Sum_Of_Functions([f1,-1*f2])
end
function -(f1::Sum_Of_Functions, f2::Piecewise_Function)
    return +(f1, -1 * f2)
end

function -(f1::Piecewise_Function, f2::Sum_Of_Functions)
    return +(f1, -1*f2)
end

function *(f1::Sum_Of_Functions,f2::Sum_Of_Functions)
    results = Array{Sum_Of_Functions}(undef, length(f1.functions_))
    for i in 1:length(f1.functions_)
        new_funcs = f1.functions_[i] * f2
        results[i] = new_funcs
    end
    return Sum_Of_Functions(results)
end
function *(f1::Sum_Of_Functions, f2::Piecewise_Function)
    return Piecewise_Function(f2.starts_, f1 .* f2.functions_)
end

function *(f1::Piecewise_Function, f2::Sum_Of_Functions)
    return *(f2,f1)
end
