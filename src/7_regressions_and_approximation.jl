
"""
    create_ols_approximation(y::Vector{<:Real}, x::Vector{<:Real}, base_x::Real = 0.0, degree::Integer = 1, intercept::Bool = true)
    create_ols_approximation(y::Vector{<:Real}, x::Vector{Q}, base_x::Union{Date,DateTime} = global_base_date, degree::Integer = 1, intercept::Bool = true) where Q<:Union{Date,DateTime,ZonedDateTime}

An approximation function calculated via OLS.

### Inputs
* `y` - A `Vector` with the y coordinates
* `x` - A `Vector` with the x coordinates
* `base_x` - A real that offsets the x. So a coordinate with x value of 2.0 will be converted to 1.8 if base_x is 0.2.
* `degree` - What the highest power of x should be. So if this is 3 then the equation will have x, x^2, x^3 as predictors.
* `intercept` - Should there be an x intercept.
### Returns
* A `Sum_Of_Functions` containing the approximation function.
"""
function create_ols_approximation(y::Vector{<:Real}, x::Vector{<:Real}, base_x::Real = 0.0, degree::Integer = 1, intercept::Bool = true)
    obs = length(y)
    if degree < 0
        error("Cannot approximate with OLS with a degree that is negative")
    end
    x = x .- base_x
    if intercept
        X = ones(obs)
        for i in 1:degree
            X = hcat(X, (x .^ i))
        end
    else
        X = x
        for i in 2:degree
            X = hcat(X, (x .^ i))
        end
    end

    lm1 = fit(LinearModel,  hcat(X), y)
    beta = lm1.pp.beta0
    func_array = Vector{PE_Function}(undef,convert(Int, intercept) + degree)
    if intercept
        func_array[1] = PE_Function(beta[1], 0.0, base_x, 0)
    end
    for d in 1:degree
        func_array[d+Integer(intercept)] = PE_Function(beta[d+Integer(intercept)], 0.0, base_x, d)
    end
    return Sum_Of_Functions(func_array)
end

function create_ols_approximation(y::Vector{<:Real}, x::Vector{Q}, base_x::WW = minimum(x), degree::Integer = 1, intercept::Bool = true) where Q<:Union{Date,DateTime,ZonedDateTime} where WW<:Union{Date,DateTime,ZonedDateTime}
    base = years_from_global_base_date.(base_x)
    xx   = years_from_global_base_date.(x)
    return create_ols_approximation(y, xx, base, degree, intercept)
end

function get_chebyshev_coefficients(chebyshev::Sum_Of_Functions, y::Vector{<:Real}, unnormalised_nodes::Vector{<:Real})
    chebyshev_on_nodes = evaluate.(Ref(chebyshev), unnormalised_nodes)
    a = sum(y .* chebyshev_on_nodes) / sum(chebyshev_on_nodes .^ 2)
    return a
end

"""
    create_chebyshev_approximation(func::Function, nodes::Integer, degree::Integer, left::Real, right::Real)

An function that will approximate another function via Chebyshev polynomials.

### Inputs
* `func` - A function that you want to approximation
* `nodes` - The number of approximation nodes
* `degree` - The degree of the Chebyshev polynomials.
* `left` - The left limit of the approximation
* `right` - The right limit of the approximation.
### Returns
* A `Sum_Of_Functions` containing the approximation function.
"""
function create_chebyshev_approximation(func::Function, nodes::Integer, degree::Integer, left::Real, right::Real)
    # This is all after Algorithm 6.2 from Judd (1998) Numerical Methods in Economics.
    if nodes <= degree
        error("Need to have more nodes than degree to use a chebyshev approximation")
    end
    k = 1:nodes
    unnormalised_nodes = -cos.( (((2 .* k) .- 1) ./ (2 * nodes)) .* pi    )
    normalised_nodes = ((unnormalised_nodes .+ 1) .* ((right-left)/2)) .+ left
    y = func.(normalised_nodes)
    chebyshevs = get_chebyshevs_up_to(degree, true)
    a = get_chebyshev_coefficients.(chebyshevs, Ref(y), Ref(unnormalised_nodes))
    transformed_chebyshevs = convert_to_linearly_rescale_inputs.(chebyshevs, (right-left)/2, +(right+left)/2)
    # Note that these alpha and beta parameters in the convert_to_linearly_rescale_inputs function differ from those in Judd because those did not work in this context.
    all_terms = a .* transformed_chebyshevs
    final_func = Sum_Of_Functions(all_terms)
    return final_func
end
