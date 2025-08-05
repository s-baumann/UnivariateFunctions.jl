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
