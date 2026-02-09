import Base.+, Base.-, Base./, Base.*
import SchumakerSpline.evaluate
function evaluate(sf::Sum_Of_Functions, point::Real)
    total = 0.0
    for func in sf.functions_
        total = total + evaluate(func,point)
    end
    return total
end
function (s::Sum_Of_Functions)(x::Q) where Q<:Union{Real,Date,DateTime,ZonedDateTime,DatePeriod}
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
    funcs = PE_Function[func * number for func in f.functions_]
    return Sum_Of_Functions(funcs)
end
function /(f::Sum_Of_Functions, number::Real)
    return *(f, 1/number )
end

function +(f1::Sum_Of_Functions, f2::Sum_Of_Functions)
    return Sum_Of_Functions([f1,f2])
end
function +(f1::Sum_Of_Functions, f2::Piecewise_Function)
    added_functions = (f1 .+ f2.functions_)::Vector
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
    terms = Vector{UnivariateFunction}()
    for fi in f1.functions_
        for fj in f2.functions_
            push!(terms, fi * fj)
        end
    end
    return Sum_Of_Functions(terms)
end
function *(f1::Sum_Of_Functions, f2::Piecewise_Function)
    return Piecewise_Function(f2.starts_, f1 .* f2.functions_)
end

function *(f1::Piecewise_Function, f2::Sum_Of_Functions)
    return *(f2,f1)
end
