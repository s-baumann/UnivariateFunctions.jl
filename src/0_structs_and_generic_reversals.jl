import Base.+, Base.-, Base./, Base.*, Base.^
import Base.sort
import SchumakerSpline.evaluate
abstract type UnivariateFunction end

struct Undefined_Function <: UnivariateFunction
end
Base.broadcastable(e::Undefined_Function) = Ref(e)
struct PE_Function{F<:Real,I<:Integer} <: UnivariateFunction
    # a exp(b(x-base)) (x-base)^d
    a_::F
    b_::F
    base_::F
    d_::I
    function PE_Function(a_::R,b_::S,base_::T,d_::I) where R<:Real where S<:Real where T<:Real where I<:Integer
        promo_type = promote_type(T,S,T)
        if d_ < 0
            error("Negative polynomial powers are not supported by this package")
            # These are banned due to the complications for calculus. Most
            # divisions are not allowed for the same reason.
        elseif (abs(b_) < tol) & (d_ == 0)
            return new{promo_type,I}(promo_type(a_),promo_type(0.0),promo_type(0.0),I(0))
        elseif abs(a_) < tol
            return new{promo_type,I}(promo_type(0.0),promo_type(0.0),promo_type(0.0),I(0))
        else
            return new{promo_type,I}(promo_type(a_),promo_type(b_),promo_type(base_),I(d_))
        end
    end
    function PE_Function(a_::R,b_::S,base_::Date,d_::I) where R<:Real where S<:Real where I<:Real
        new_base_ = years_from_global_base(base_)
        return PE_Function(a_, b_, new_base_, d_)
    end
    function PE_Function{Q,N}(a_::R,b_::S,base_::T,d_::I) where Q<:Real where N<:Integer where R<:Real where S<:Real where T<:Real where I<:Integer
        return new{Q,N}(Q(a_), Q(b_), Q(base_), N(d_))
    end
end
Base.broadcastable(e::PE_Function) = Ref(e)
struct Sum_Of_Functions <: UnivariateFunction
    functions_::Union{Array{PE_Function,1},Array{PE_Function{<:Real,<:Integer},1}}
    function Sum_Of_Functions(funcs::Sum_Of_Functions)
        return funcs
    end
    function Sum_Of_Functions(funcs::AbstractArray)
        undefined_funcs  = funcs[isa.(funcs, Ref(Undefined_Function))]
        if length(undefined_funcs) > 0
            return Undefined_Function()
        end
        functions_ = clean_array_of_functions(funcs)
        return new(functions_)
    end
end
Base.broadcastable(e::Sum_Of_Functions) = Ref(e)
struct Piecewise_Function <: UnivariateFunction
    starts_::Array{<:Real,1}
    functions_::Array{UnivariateFunction,1}
    function Piecewise_Function(starts::Array{<:Real,1}, functions::Array)
        if !issorted(starts)
            error("Piecewise_Function must be constructed with a sorted increasing list of starts and the corresponding functions.")
        end
        if starts[1] != -Inf
            starts    = vcat(-Inf, starts)
            functions = vcat(Undefined_Function(), functions)
        end
        len = length(starts)
        if starts[len] == Inf
            starts    = starts[1:(len-1)]
            functions = functions[1:(len-1)]
        end
        starts_without_pw_parts, functions_without_pw_parts = deal_with_piecewise_inputs(starts, functions)
        new(starts_without_pw_parts, functions_without_pw_parts)
    end
    function Piecewise_Function(starts::Array{Date,1}, functions::Array)
        starts_ = years_from_global_base.(starts)
        return Piecewise_Function(starts_, functions)
    end
    function Piecewise_Function(start::Real, func)
        return Piecewise_Function([start], [func])
    end
    function Piecewise_Function(start::Date, func)
        start_float = years_from_global_base(start)
        return Piecewise_Function(start_float, func)
    end
end
Base.broadcastable(e::Piecewise_Function) = Ref(e)



function sort(funcs::Array{PE_Function,1})
    len = length(funcs)
    attributes = hcat(map(f -> f.d_, funcs), map(f -> f.base_, funcs), map(f -> f.b_, funcs) )
    ordering = convert(Array{Int}, sortslices(hcat(attributes, 1:len) , dims = 1)[:,4])
    return funcs[ordering]
end

function trim_piecewise_function(func::Piecewise_Function, left_limit::Real, right_limit::Real)
    if func.starts_[1] != -Inf
        starts_ = [-Inf, left_limit, right_limit]
        funcs_  = [Undefined_Function(), func, Undefined_Function()]
    else
        starts_ = [left_limit, right_limit]
        funcs_  = [func, Undefined_Function()]
    end
    return Piecewise_Function(starts_, funcs_)

end

function deal_with_piecewise_inputs(starts::Array{<:Real,1}, functions::Array)
    where_are_pw_bits = findall(typeof.(functions) .== UnivariateFunctions.Piecewise_Function)
    if length(where_are_pw_bits) < 1
        return starts, functions
    else
        first_pw_slice_i = where_are_pw_bits[1]
        pw = functions[first_pw_slice_i]
        from = starts[first_pw_slice_i]
        if length(starts) > first_pw_slice_i
            to = starts[first_pw_slice_i+1]
        else
            to = Inf
        end
        pw_starts, pw_funcs = take_piecewise_slice(pw.starts_, pw.functions_,from, to)
        pw_starts[1] = from
        new_starts = vcat(first_entries(starts, first_pw_slice_i-1), pw_starts, last_entries(starts, length(starts) - first_pw_slice_i))
        new_funcs = vcat(first_entries(functions, first_pw_slice_i-1), pw_funcs, last_entries(functions, length(starts) - first_pw_slice_i))
        return deal_with_piecewise_inputs(new_starts, new_funcs)
    end

end

function first_entries(vec::Array, how_many::Integer)
    if how_many < 1
        return Array{Any,1}()
    elseif how_many > length(vec)
        return vec
    else
        return vec[1:how_many]
    end
end

function last_entries(vec::Array, how_many::Integer)
    if how_many < 1
        return Array{Any,1}()
    elseif how_many > length(vec)
        return vec
    else
        leng =  length(vec)
        return vec[(leng-how_many+1):leng]
    end
end

function take_piecewise_slice(starts::Array{<:Real,1}, functions::Array, from::Real, to::Real)
    from_i = searchsortedlast(starts, from)
    to_i = searchsortedlast(starts, to)
    return starts[from_i:to_i], functions[from_i:to_i]
end

function find(a::BitArray)
    len = length(a)
    indices = Array{Int,1}()
    for i in 1:len
        if a[i]
            append!(indices, i)
        end
    end
    return indices
end

function rationalise_array_of_functions(funcs::Array{PE_Function,1})
    len = length(funcs)
    if len < 2
        return funcs
    end
    attributes = hcat(map(f -> f.b_, funcs), map(f -> f.base_, funcs), map(f -> f.d_, funcs) )
    multipliers = map(f -> f.a_, funcs)
    unique_attributes = unique(attributes, dims = 1)
    new_len = size(unique_attributes)[1]
    new_funcs = Array{PE_Function,1}(undef,new_len)
    for i in 1:new_len
        atts =  transpose(unique_attributes[i,:])
        which_multipliers = all(attributes .== atts, dims=2)
        multiplier = sum(multipliers[find(which_multipliers)])
        new_funcs[i] = PE_Function(multiplier, atts[1], atts[2], convert(Int, atts[3]))
    end
    return new_funcs
end

function clean_array_of_functions(funcs::Array)
    undefined_funcs  = funcs[isa.(funcs, Ref(Undefined_Function))]
    piecewise_funcs  = funcs[isa.(funcs, Ref(Piecewise_Function))]
    if length(undefined_funcs) > 0
        return Undefined_Function()
    elseif length(piecewise_funcs) > 0
        error("You cannot directly construct a Sum_Of_Functions with a Piecewise_Function in the input array. Instead add up the piecewise functions directly, for instance typing 'f1 + f2' for the two piecewise functions.")
    end
    pe_funcs  = funcs[isa.(funcs, Ref(PE_Function))]
    pe_funcs = pe_funcs[abs.(map( x -> x.a_, pe_funcs)) .>= eps()]
    sum_funcs = funcs[isa.(funcs, Ref(Sum_Of_Functions))]
    if length(sum_funcs) > 0
        for sf in sum_funcs
            pe_funcs_in_sf = clean_array_of_functions(sf.functions_)
            pe_funcs       = vcat(pe_funcs, pe_funcs_in_sf)
        end
    end
    if length(pe_funcs) == 0
        return [PE_Function(0.0,0.0,0.0,0)]
    end
    simplified_functions = rationalise_array_of_functions(convert(Array{PE_Function,1}, pe_funcs))
    return simplified_functions
end




function +(number::Real, f::UnivariateFunction)
    return +(f,number)
end
function -(number::Real, f::UnivariateFunction)
    return +(number, -1*f)
end
function *(number::Real, f::UnivariateFunction)
    return *(f,number)
end
function /(number::Real, f::UnivariateFunction)
    error("It is not possible yet to divide scalars by UnivariateFunctions")
end
function ^(number::Real, f::UnivariateFunction)
    error("It is not possible yet to raise to the power of a UnivariateFunctions")
end

function evaluate(f::UnivariateFunction, d::Date)
    date_in_relation_to_global_base = years_from_global_base(d)
    return evaluate(f, date_in_relation_to_global_base)
end
function evaluate(f::UnivariateFunction, x::DatePeriod)
    return evaluate(f, period_length(x))
end

function ^(f1::UnivariateFunction,num::Integer) # This will get overridden for undefined and zeros.
    if num < 0
        error("Cannot raise any univariate function to a negative power")
    elseif num == 0
        return PE_Function(1.0,0.0,0.0,0)
    elseif num == 1
        return f1
    elseif num == 2
        return f1 * f1
    else
        product = f1 * f1
        for i in 1:(num-2)
            product = product * f1
        end
        return product
    end
end


function change_base_of_PE_Function(f::PE_Function, new_base::Real)
    old_base = f.base_
    diff = new_base - old_base
    if abs(diff) < tol
        return f
    end
    # First the exponential part.
    new_a = f.a_ * exp(f.b_*diff)
    if new_a < tol
        error("Underflow problem. Changing to this base cannot be done")
    end
    # Now the polynomial part.
    if f.d_ == 0
        return PE_Function(new_a, f.b_, new_base, 0)
    else
        n = f.d_
        funcs = Array{UnivariateFunction,1}(undef,n+1)
        for r in 0:n
            binom_coeff = factorial(n) / (factorial(r) * factorial(n-r))
            new_multiplier = binom_coeff * new_a * diff^r
            new_func = PE_Function(new_multiplier, f.b_, new_base, n-r )
            funcs[r+1] = new_func
        end
        return Sum_Of_Functions(funcs)
    end
end

# Conversions for linearly rescaling inputs.
function convert_to_linearly_rescale_inputs(f::Undefined_Function, alpha::Real, beta::Real)
    return f
end
function convert_to_linearly_rescale_inputs(f::PE_Function, alpha::Real, beta::Real)
    # We want the change the function so that whenever we put in x it is like we put in alpha x + beta.
    beta = beta / alpha
    alpha = 1.0/alpha
    new_base_ = (f.base_ + beta)/alpha
    new_multiplier = f.a_ * alpha^(f.d_)
    new_power_ = f.b_ * alpha
    return PE_Function(new_multiplier, new_power_, new_base_, f.d_)
end
function convert_to_linearly_rescale_inputs(f::Sum_Of_Functions, alpha::Real, beta::Real)
    funcs = convert_to_linearly_rescale_inputs.(f.functions_, alpha,beta)
    return Sum_Of_Functions(funcs)
end
function convert_to_linearly_rescale_inputs(f::Piecewise_Function, alpha::Real, beta::Real)
    new_starts_ = (f.starts_ .* alpha) .+ beta
    funcs_ = convert_to_linearly_rescale_inputs.(f.functions_, alpha,beta)
    return Piecewise_Function(new_starts_, funcs_)
end
