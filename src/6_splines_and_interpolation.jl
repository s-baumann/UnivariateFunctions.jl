"""
    create_quadratic_spline(x::Vector{Q},y::Vector{<:Real} ; gradients::Union{Missing,Vector{<:Real}} = missing, extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve),
                                 left_gradient::Union{Missing,Real} = missing, right_gradient::Union{Missing,Real} = missing) where Q<:Union{Date,DateTime,ZonedDateTime}
    create_quadratic_spline(x::Vector{<:Real},y::Vector{<:Real} ; gradients::Union{Missing,Vector{<:Real}} = missing, extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve),
                                 left_gradient::Union{Missing,Real} = missing, right_gradient::Union{Missing,Real} = missing)
    create_quadratic_spline(x::Union{Vector{D},Vector{<:DatePeriod}},y::Vector{<:Real}; gradients::Union{Missing,Vector{<:Real}} = missing, extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve),
                                 left_gradient::Union{Missing,Real} = missing, right_gradient::Union{Missing,Real} = missing) where D<:DatePeriod
Makes a quadratic shape-preserving interpolation spline using the SchumakerSpline.jl package. This is returned as a `Piecewise_Function` rather than as a `Schumaker` struct.

### Inputs
* `x` - A `Vector` with the x coordinates
* `y` - A `Vector` with the y coordinates
* `gradients` - A `Vector` with the gradiants at each x location. This is calculated if not provided.
* `extrapolation` - A tuple of enum value describing how to extrapolate (on the left and right sides).
* `left_gradient` - The gradiant to impose on the left edge (ie the first x coordinate).
* `right_gradient` - The gradiant to impose on the right edge (ie the last x coordinate).
### Returns
* A `Piecewise_Function` containing the spline.


    create_quadratic_spline(schum::Schumaker)
This converts a spline represented by a `SchumakerSpline.Schumaker` struct into the same
spline but represented by a `Piecewise_Function`.
### Inputs
* `schum` - A `Schumaker` struct.
### Returns
* A `Piecewise_Function` containing the spline.
"""
function create_quadratic_spline(x::Vector{Q},y::Vector{<:Real} ; gradients::Union{Missing,Vector{<:Real}} = missing, extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve), left_gradient::Union{Missing,Real} = missing, right_gradient::Union{Missing,Real} = missing) where Q<:Union{Date,DateTime,ZonedDateTime}
    x_as_Floats = years_from_global_base_date.(x)
    return create_quadratic_spline(x_as_Floats, y; gradients = gradients, extrapolation = extrapolation, left_gradient = left_gradient, right_gradient = right_gradient)
end

function create_quadratic_spline(x::Vector{<:Real},y::Vector{<:Real} ; gradients::Union{Missing,Vector{<:Real}} = missing, extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve),
                                 left_gradient::Union{Missing,Real} = missing, right_gradient::Union{Missing,Real} = missing)
    schum = Schumaker(x, y; gradients = gradients, extrapolation = extrapolation, left_gradient = left_gradient, right_gradient = right_gradient)
    return create_quadratic_spline(schum)
end
function create_quadratic_spline(schum::Schumaker)
    starts_ = schum.IntStarts_
    coefficients = schum.coefficient_matrix_
    number_of_intervals = size(coefficients)[1]
    funcs_ = Vector{Sum_Of_Functions}(undef, number_of_intervals)
    for i in 1:number_of_intervals
        quadratic = PE_Function(coefficients[i,1], 0.0, starts_[i], 2)
        linear    = PE_Function(coefficients[i,2], 0.0, starts_[i], 1)
        constant  = PE_Function(coefficients[i,3], 0.0, 0.0       , 0)
        polynomial = Sum_Of_Functions([quadratic, linear, constant])
        funcs_[i] = polynomial
    end
    return Piecewise_Function(starts_, funcs_)
end
function create_quadratic_spline(x::Union{Vector{D},Vector{<:DatePeriod}},y::Vector{<:Real}; gradients::Union{Missing,Vector{<:Real}} = missing, extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve),
                                 left_gradient::Union{Missing,Real} = missing, right_gradient::Union{Missing,Real} = missing) where D<:DatePeriod
    x_as_Floats = period_length.(x)
    return create_quadratic_spline(x_as_Floats, y; gradients = gradients, extrapolation = extrapolation, left_gradient = left_gradient, right_gradient = right_gradient)
end


"""
    create_constant_interpolation_to_right(x::Vector{Date},y::Vector{<:Real})
    create_constant_interpolation_to_right(x::Vector{<:Real},y::Vector{<:Real})
    create_constant_interpolation_to_right(x::Union{Vector{D},Vector{<:DatePeriod}},y::Vector{<:Real}) where D<:DatePeriod

Makes a piecewise constant interpolation function. values from the left are copied to the right.

### Inputs
* `x` - A `Vector` with the x coordinates
* `y` - A `Vector` with the y coordinates
### Returns
* A `Piecewise_Function` containing the interpolation function.
"""
function create_constant_interpolation_to_right(x::Vector{Date},y::Vector{<:Real})
    x_Float = years_from_global_base_date.(x)
    return create_constant_interpolation_to_right(x_Float,y)
end

function create_constant_interpolation_to_right(x::Vector{<:Real},y::Vector{<:Real})
    x_ = vcat(-Inf,x)
    y = vcat(y[1], y)
    funcs_ = PE_Function.(y,0.0,0.0,0)
    return Piecewise_Function(x_, funcs_)
end
function create_constant_interpolation_to_right(x::Union{Vector{D},Vector{<:DatePeriod}},y::Vector{<:Real}) where D<:DatePeriod
    x_as_Floats = period_length.(x)
    return create_constant_interpolation_to_right(x_as_Floats,y)
end

"""
    create_constant_interpolation_to_left(x::Vector{Date},y::Vector{<:Real})
    create_constant_interpolation_to_left(x::Vector{<:Real},y::Vector{<:Real})
    create_constant_interpolation_to_left(x::Union{Vector{D},Vector{<:DatePeriod}},y::Vector{<:Real}) where D<:DatePeriod

Makes a piecewise constant interpolation function. values from the right are copied to the left.

### Inputs
* `x` - A `Vector` with the x coordinates
* `y` - A `Vector` with the y coordinates
### Returns
* A `Piecewise_Function` containing the interpolation function.
"""
function create_constant_interpolation_to_left(x::Vector{Date},y::Vector{<:Real})
    x_Float = years_from_global_base_date.(x)
    return create_constant_interpolation_to_left(x_Float,y)
end

function create_constant_interpolation_to_left(x::Vector{<:Real},y::Vector{<:Real})
    x_ = vcat(-Inf,x[1:(length(x)-1)])
    funcs_ = PE_Function.(y,0.0,0.0,0)
    return Piecewise_Function(x_, funcs_)
end

function create_constant_interpolation_to_left(x::Union{Vector{D},Vector{<:DatePeriod}},y::Vector{<:Real}) where D<:DatePeriod
    x_as_Floats = period_length.(x)
    return create_constant_interpolation_to_left(x_as_Floats,y)
end

"""
    create_linear_interpolation(x::Vector{Date},y::Vector{<:Real})
    create_linear_interpolation(x::Union{Vector{D},Vector{<:DatePeriod}},y::Vector{<:Real}) where D<:DatePeriod
    create_linear_interpolation(x::Vector{R},y::Vector{<:Real}) where R<:Real
Makes a piecewise linear interpolation function. This is continuous.

### Inputs
* `x` - A `Vector` with the x coordinates
* `y` - A `Vector` with the y coordinates
### Returns
* A `Piecewise_Function` containing the interpolation function.
"""
function create_linear_interpolation(x::Vector{Q},y::Vector{<:Real}) where Q<:Union{Date,DateTime,ZonedDateTime}
    x_Float = years_from_global_base_date.(x)
    return create_linear_interpolation(x_Float,y)
end

function create_linear_interpolation(x::Union{Vector{D},Vector{<:DatePeriod}},y::Vector{<:Real}) where D<:DatePeriod
    x_as_Floats = period_length.(x)
    return create_linear_interpolation(x_as_Floats,y)
end

function create_linear_interpolation(x::Vector{R},y::Vector{<:Real}) where R<:Real
    len = length(x)
    if len < 2
        error("Insufficient data to linearly interpolate")
    end
    starts_ = Vector{R}(undef, len-1)
    funcs_  = Vector{UnivariateFunction}(undef, len-1)
    coefficient = (y[2] - y[1])/(x[2] - x[1])
    starts_[1] = -Inf
    funcs_[1]  = Sum_Of_Functions([PE_Function(y[1],0.0,0.0,0), PE_Function(coefficient,0.0,x[1],1)])
    if len > 2
        for i in 2:(len-1)
            starts_[i] = x[i]
            coefficient = (y[i+1] - y[i])/(x[i+1] - x[i])
            funcs_[i]  = Sum_Of_Functions([PE_Function(y[i],0.0,0.0,0), PE_Function(coefficient,0.0,x[i],1)])
        end
    end
    return Piecewise_Function(starts_, funcs_)
end

"""
    simplify(f::Piecewise_Function, n_points::Integer, left::Real, right::Real, method::Symbol = :linear)
    simplify(f::Piecewise_Function, n_points::Integer, left::Q, right::Q, method::Symbol = :linear) where Q<:Union{Date,DateTime,ZonedDateTime}

Approximates a `Piecewise_Function` with a simpler one by evaluating at `n_points` evenly
spaced points over `[left, right]` and re-interpolating.

### Inputs
* `f` - A `Piecewise_Function` to simplify.
* `n_points` - The number of evenly spaced sample points.
* `left` - The left boundary of the domain.
* `right` - The right boundary of the domain.
* `method` - The interpolation method. One of `:constant_right`, `:constant_left`, `:linear`, `:quadratic`.

### Returns
* A `Piecewise_Function`.
"""
function simplify(f::Piecewise_Function, n_points::Integer, left::Real, right::Real, method::Symbol = :linear)
    if n_points < 2
        error("n_points must be at least 2")
    end
    x = collect(range(left, right, length = n_points))
    y = evaluate.(f, x)
    if method == :constant_right
        return create_constant_interpolation_to_right(x, y)
    elseif method == :constant_left
        return create_constant_interpolation_to_left(x, y)
    elseif method == :linear
        return create_linear_interpolation(x, y)
    elseif method == :quadratic
        return create_quadratic_spline(x, y)
    else
        error("Unknown method :$method. Use :constant_right, :constant_left, :linear, or :quadratic.")
    end
end

function simplify(f::Piecewise_Function, n_points::Integer, left::Q, right::Q, method::Symbol = :linear) where Q<:Union{Date,DateTime,ZonedDateTime}
    return simplify(f, n_points, years_from_global_base_date(left), years_from_global_base_date(right), method)
end


"""
    UnivariateFitter

A mutable struct for iteratively fitting univariate functions with shape constraints.
Each call to `fit!` fits new data and optionally blends with the previous fit via
`weight_on_new`. Periodic simplification reduces the complexity of the accumulated
`Piecewise_Function`.

### Fields
* `fun` - The current fitted `UnivariateFunction`.
* `method` - Fitting method. One of `:increasing`, `:decreasing`, `:convex`, `:concave`, `:quasiconvex`, `:quasiconcave`, `:supersmoother`.
* `times_through` - Number of times `fit!` has been called.
* `simplification_frequency` - Simplify the function every this many calls to `fit!`. `0` disables simplification.
* `nbins` - Number of bins for the regression and for simplification.
* `equally_spaced_bins` - If `true`, bins are equally spaced in x; if `false`, based on observation quantiles.
* `weight_on_new` - Blending weight in `[0,1]`. `1.0` means only the new fit is used; lower values blend with the previous fit.
* `left_` - Left boundary for simplification domain.
* `right_` - Right boundary for simplification domain.

### Constructor
    UnivariateFitter(method::Symbol; nbins=10, equally_spaced_bins=true,
                     weight_on_new=1.0, simplification_frequency=0,
                     left=-Inf, right=Inf)

Creates a `UnivariateFitter` initialised with a zero function.
"""
mutable struct UnivariateFitter
    fun::UnivariateFunctions.UnivariateFunction
    method::Symbol
    times_through::Int
    simplification_frequency::Int
    nbins::Int
    equally_spaced_bins::Bool
    weight_on_new::Float64
    left_::Float64
    right_::Float64
end

function UnivariateFitter(method::Symbol; nbins::Int=10, equally_spaced_bins::Bool=true,
                          weight_on_new::Float64=1.0, simplification_frequency::Int=0,
                          left::Real=-Inf, right::Real=Inf)
    return UnivariateFitter(PE_Function(0.0, 0.0, 0.0, 0), method, 0, simplification_frequency,
                            nbins, equally_spaced_bins, weight_on_new, Float64(left), Float64(right))
end

function evaluate(fitter::UnivariateFitter, x)
    return fitter.fun(x)
end
function (fitter::UnivariateFitter)(x)
    return evaluate(fitter, x)
end

"""
    fit!(fitter::UnivariateFitter, x_new, y_new)

Fit the `UnivariateFitter` to new data `x_new`, `y_new`. The fitted function is blended
with the previous fit according to `fitter.weight_on_new`. If `simplification_frequency > 0`
and the current iteration is a multiple, the function is simplified via resampling.
"""
function fit!(fitter::UnivariateFitter, x_new, y_new)
    newfun = begin
        if fitter.method == :increasing
            UnivariateFunctions.monotonic_regression(x_new, y_new; nbins=fitter.nbins, equally_spaced_bins=fitter.equally_spaced_bins, increasing=true)
        elseif fitter.method == :decreasing
            UnivariateFunctions.monotonic_regression(x_new, y_new; nbins=fitter.nbins, equally_spaced_bins=fitter.equally_spaced_bins, increasing=false)
        elseif fitter.method == :convex
            UnivariateFunctions.unimodal_regression(x_new, y_new; nbins=fitter.nbins, equally_spaced_bins=fitter.equally_spaced_bins, convex=true, quasi=false)
        elseif fitter.method == :concave
            UnivariateFunctions.unimodal_regression(x_new, y_new; nbins=fitter.nbins, equally_spaced_bins=fitter.equally_spaced_bins, convex=false, quasi=false)
        elseif fitter.method == :quasiconvex
            UnivariateFunctions.unimodal_regression(x_new, y_new; nbins=fitter.nbins, equally_spaced_bins=fitter.equally_spaced_bins, convex=true, quasi=true)
        elseif fitter.method == :quasiconcave
            UnivariateFunctions.unimodal_regression(x_new, y_new; nbins=fitter.nbins, equally_spaced_bins=fitter.equally_spaced_bins, convex=false, quasi=true)
        elseif fitter.method == :supersmoother
            UnivariateFunctions.supersmoother(x_new, y_new)
        else
            error("Unknown maximal correlation method: $(fitter.method)")
        end
    end
    if fitter.times_through > 0
        newfun = (fitter.weight_on_new * newfun) + ((1.0 - fitter.weight_on_new) * fitter.fun)
    end
    if fitter.simplification_frequency > 0 && (fitter.times_through % fitter.simplification_frequency == 0)
        newfun = UnivariateFunctions.simplify(newfun, fitter.nbins, fitter.left_, fitter.right_)
    end
    fitter.fun = newfun
    fitter.times_through += 1
end
