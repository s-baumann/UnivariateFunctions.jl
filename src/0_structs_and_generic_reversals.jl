import Base.+, Base.-, Base./, Base.*, Base.^
import SchumakerSpline.evaluate

"""
    UnivariateFunction
An abstract type. The concrete structs that have been implemented are Undefined_Function,
    PE_Function, Sum_Of_Functions, Piecewise_Function.
"""
abstract type UnivariateFunction end

"""
    Undefined_Function <: UnivariateFunction

This function throws an error if you ever try to evaluate it.
Think of it as doing the role of missing but for UnivariateFunctions
"""
struct Undefined_Function <: UnivariateFunction
end
Base.broadcastable(e::Undefined_Function) = Ref(e)

"""
    PE_Function{F<:Real,I<:Integer} <: UnivariateFunction

This function has the functional form:
    a exp(b(x-base)) (x-base)^d
Where a,b,base are floats and d is a positive integer. These four are the
 members of the struct.
"""
struct PE_Function{F<:Real,I<:Integer} <: UnivariateFunction
    # a exp(b(x-base)) (x-base)^d
    a_::F
    b_::F
    base_::F
    d_::I
    function PE_Function(a_::Real)
        return PE_Function(a_, 0.0, 0.0, 0)
    end
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
    function PE_Function(a_::R,b_::S,base_::Q,d_::I) where R<:Real where S<:Real where I<:Real where Q<:Union{Date,DateTime,ZonedDateTime}
        new_base_ = years_from_global_base_date(base_)
        return PE_Function(a_, b_, new_base_, d_)
    end
    function PE_Function{Q,N}(a_::R,b_::S,base_::T,d_::I) where Q<:Real where N<:Integer where R<:Real where S<:Real where T<:Real where I<:Integer
        return new{Q,N}(Q(a_), Q(b_), Q(base_), N(d_))
    end
end
Base.broadcastable(e::PE_Function) = Ref(e)

"""
    Sum_Of_Functions <: UnivariateFunction

This function contants a vector of UnivariateFunctions. When evaluted it
adds the evaluations of these functions and returns the sum.
"""
struct Sum_Of_Functions <: UnivariateFunction
    functions_::Union{Vector{PE_Function},Vector{PE_Function{<:Real,<:Integer}}}
    function Sum_Of_Functions(funcs::Sum_Of_Functions)
        return funcs
    end
    function Sum_Of_Functions(funcs::Vector)
        undefined_funcs  = funcs[isa.(funcs, Ref(Undefined_Function))]
        if length(undefined_funcs) > 0
            return Undefined_Function()
        end
        functions_ = clean_array_of_functions(funcs)
        return new(functions_)
    end
end
Base.broadcastable(e::Sum_Of_Functions) = Ref(e)

"""
    Piecewise_Function <: UnivariateFunction

This function contants a vector of locations in the x space and a vector of UnivariateFunctions.
When evaludated it uses these vectors as a lookup. It chooses the correct UnivariateFunction
and evaluates it.
"""
struct Piecewise_Function <: UnivariateFunction
    starts_::Vector{<:Real}
    functions_::Vector{UnivariateFunction}
    function Piecewise_Function(starts::Vector{<:Real}, functions::Vector)
        if !issorted(starts)
            error("Piecewise_Function must be constructed with a sorted increasing list of starts and the corresponding functions.")
        end
        len = length(starts)
        if starts[len] == Inf
            starts    = starts[1:(len-1)]
            functions = functions[1:(len-1)]
        end
        starts_without_pw_parts, functions_without_pw_parts = deal_with_piecewise_inputs(starts, functions)
        new(starts_without_pw_parts, functions_without_pw_parts)
    end
    function Piecewise_Function(starts::Vector{Q}, functions::Vector) where Q<:Union{Date,DateTime,ZonedDateTime}
        starts_ = years_from_global_base_date.(starts)
        return Piecewise_Function(starts_, functions)
    end
    function Piecewise_Function(start::R, func) where R<:Real
        return Piecewise_Function(R[start], [func])
    end
    function Piecewise_Function(start::Q, func) where Q<:Union{Date,DateTime,ZonedDateTime}
        start_float = years_from_global_base_date(start)
        return Piecewise_Function(start_float, func)
    end
end
Base.broadcastable(e::Piecewise_Function) = Ref(e)

"""
    trim_piecewise_function(func::Piecewise_Function, left_limit::Real, right_limit::Real)
This trims the end of a piecewise function. So if there is a piecewise function
with support between -10,10 then you can trim it to only have support between
-5 and 5. Then if it is evaluated outside -5 to 5 it will be undefined.

### Inputs
* `func` - A Piecewise_Function.
* `left_limit` - The left limit.
* `right_limit` - The right limit.
### Returns
* A `Piecewise_Function`.
"""
function trim_piecewise_function(func::Piecewise_Function, left_limit::R, right_limit::Real) where R<:AbstractFloat
    starts_ = R[left_limit, right_limit]
    funcs_  = UnivariateFunction[func, Undefined_Function()]
    return Piecewise_Function(starts_, funcs_)
end

function deal_with_piecewise_inputs(starts::Vector{R}, functions::Vector) where R<:Real
    where_are_pw_bits = findall(typeof.(functions) .== UnivariateFunctions.Piecewise_Function)
    if length(where_are_pw_bits) < 1
        return starts, functions
    else
        stt = Vector{R}()
        fns = Vector{UnivariateFunction}()
        for i in 1:length(starts)
            s = starts[i]
            terminal = i < length(starts) ? starts[i+1] : Inf
            f = functions[i]
            if isa(f, UnivariateFunctions.Piecewise_Function)
                # If we are before the start of the piecewise insert an Undefined_Function.
                if s < f.starts_[1]
                    stt = push!(stt, s)
                    fns = push!(fns, Undefined_Function())
                end
                # Adding on the starts and functions of the Piecewise_Function.
                for j in 1:length(f.starts_)
                    putative_start = max(f.starts_[j], s)
                    putative_endd = j < length(f.starts_) ? min(f.starts_[j+1], terminal) : terminal
                    if putative_endd - putative_start > 10 * eps()
                        stt = push!(stt, putative_start)
                        fns = push!(fns, f.functions_[j])
                    end
                end
            else
                stt = push!(stt, s)
                fns = push!(fns, f)
            end
        end
        return stt, fns
    end
end

function first_entries(vec::Vector, how_many::Integer)
    if how_many < 11
        return Vector{Any}()
    elseif how_many > length(vec)
        return vec
    else
        return vec[1:how_many]
    end
end

function last_entries(vec::Vector, how_many::Integer)
    if how_many < 1
        return Vector{Any}()
    elseif how_many > length(vec)
        return vec
    else
        leng =  length(vec)
        return vec[(leng-how_many+1):leng]
    end
end

function take_piecewise_slice(starts::Vector{<:Real}, functions::Vector, from::Real, to::Real)
    from_i = searchsortedlast(starts, from)
    to_i = searchsortedlast(starts, to)
    return starts[from_i:to_i], functions[from_i:to_i]
end

function find(a::BitArray)
    len = length(a)
    indices = Vector{Int}()
    for i in 1:len
        if a[i]
            append!(indices, i)
        end
    end
    return indices
end

function rationalise_array_of_functions(funcs::Vector{PE_Function})
    len = length(funcs)
    if len < 2
        return funcs
    end
    attributes = hcat(map(f -> f.b_, funcs), map(f -> f.base_, funcs), map(f -> f.d_, funcs) )
    multipliers = map(f -> f.a_, funcs)
    unique_attributes = unique(attributes, dims = 1)
    new_len = size(unique_attributes)[1]
    new_funcs = Vector{PE_Function}(undef,new_len)
    for i in 1:new_len
        atts =  transpose(unique_attributes[i,:])
        which_multipliers = all(attributes .== atts, dims=2)
        multiplier = sum(multipliers[find(which_multipliers)])
        new_funcs[i] = PE_Function(multiplier, atts[1], atts[2], convert(Int, atts[3]))
    end
    return new_funcs
end

function clean_array_of_functions(funcs::Vector)
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
    simplified_functions = rationalise_array_of_functions(convert(Vector{PE_Function}, pe_funcs))
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

"""
    evaluate(f::UnivariateFunction, r::Real)
    evaluate(f::UnivariateFunction, d::Q) where Q<:Union{Date,DateTime,ZonedDateTime}
    evaluate(f::UnivariateFunction, x::DatePeriod)

This evaluates the function at the requested point. If a `Date`, `DateTime` is input
then it is first converted to a scalar with the `years_from_global_base_date` function.
`DatePeriod`s are converted with the `period_length` function.
"""
function evaluate(f::UnivariateFunction, d::Q) where Q<:Union{Date,DateTime,ZonedDateTime}
    date_in_relation_to_global_base = years_from_global_base_date(d)
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

"""
    change_base_of_PE_Function(f::PE_Function, new_base::Real)
This changes the base of a `PE_Function`. So if the base was 2 then it can be
converted to 3 with an additional constant term.

### Inputs
* `f` - A `PE_Function`.
* `new_base` - The new base.
### Returns
* A `PE_Function` or a `Sum_Of_Functions`.
"""
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
        funcs = Vector{UnivariateFunction}(undef,n+1)
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
"""
    convert_to_linearly_rescale_inputs(f::UnivariateFunction, alpha::Real, beta::Real)
This alters a function so that whenever we put in x it is like we put in `alpha * x + beta`.

### Inputs
* `f` - A `UnivariateFunction`.
* `alpha` - The slope of the rescaling.
* `beta` - The level of the rescaling.
### Returns
* A `UnivariateFunction` of the type that you input to the function.
"""
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
