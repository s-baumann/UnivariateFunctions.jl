var documenterSearchIndex = {"docs":
[{"location":"api/","page":"API","title":"API","text":"CurrentModule = UnivariateFunctions","category":"page"},{"location":"api/#Internal-Functions","page":"API","title":"Internal Functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Pages = [\"api.md\"]","category":"page"},{"location":"api/#Main-Structs","page":"API","title":"Main Structs","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"UnivariateFunction\nUndefined_Function\nPE_Function\nSum_Of_Functions\nPiecewise_Function","category":"page"},{"location":"api/#UnivariateFunctions.UnivariateFunction","page":"API","title":"UnivariateFunctions.UnivariateFunction","text":"UnivariateFunction\n\nAn abstract type. The concrete structs that have been implemented are UndefinedFunction,     PEFunction, SumOfFunctions, Piecewise_Function.\n\n\n\n\n\n","category":"type"},{"location":"api/#UnivariateFunctions.Undefined_Function","page":"API","title":"UnivariateFunctions.Undefined_Function","text":"Undefined_Function <: UnivariateFunction\n\nThis function throws an error if you ever try to evaluate it. Think of it as doing the role of missing but for UnivariateFunctions\n\n\n\n\n\n","category":"type"},{"location":"api/#UnivariateFunctions.PE_Function","page":"API","title":"UnivariateFunctions.PE_Function","text":"PE_Function{F<:Real,I<:Integer} <: UnivariateFunction\n\nThis function has the functional form:     a exp(b(x-base)) (x-base)^d Where a,b,base are floats and d is a positive integer. These four are the  members of the struct.\n\n\n\n\n\n","category":"type"},{"location":"api/#UnivariateFunctions.Sum_Of_Functions","page":"API","title":"UnivariateFunctions.Sum_Of_Functions","text":"Sum_Of_Functions <: UnivariateFunction\n\nThis function contants a vector of UnivariateFunctions. When evaluted it adds the evaluations of these functions and returns the sum.\n\n\n\n\n\n","category":"type"},{"location":"api/#UnivariateFunctions.Piecewise_Function","page":"API","title":"UnivariateFunctions.Piecewise_Function","text":"Piecewise_Function <: UnivariateFunction\n\nThis function contants a vector of locations in the x space and a vector of UnivariateFunctions. When evaludated it uses these vectors as a lookup. It chooses the correct UnivariateFunction and evaluates it.\n\n\n\n\n\n","category":"type"},{"location":"api/#Function-Evaluation-and-Calculus","page":"API","title":"Function Evaluation and Calculus","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Note that in addition to the below functions the following operators:","category":"page"},{"location":"api/","page":"API","title":"API","text":"+, -, \\, *, ^","category":"page"},{"location":"api/","page":"API","title":"API","text":"have also been overloaded so that a function will be returned with the analytical sum, difference, product, quotient, power. The restrictions are that you cannot divide by a function (although you can divide by a scalar) and only positive integer powers can be taken.","category":"page"},{"location":"api/","page":"API","title":"API","text":"evaluate\nderivative\nindefinite_integral\nevaluate_integral\nright_integral\nleft_integral","category":"page"},{"location":"api/#SchumakerSpline.evaluate","page":"API","title":"SchumakerSpline.evaluate","text":"evaluate(f::UnivariateFunction, r::Real)\nevaluate(f::UnivariateFunction, d::Union{Date,DateTime})\nevaluate(f::UnivariateFunction, x::DatePeriod)\n\nThis evaluates the function at the requested point. If a Date, DateTime is input then it is first converted to a scalar with the yearsfromglobalbase function. DatePeriods are converted with the periodlength function.\n\n\n\n\n\n","category":"function"},{"location":"api/#UnivariateFunctions.derivative","page":"API","title":"UnivariateFunctions.derivative","text":"derivative(f::UnivariateFunction)\n\nThis calculates the derivative of the function and returns it as a UnivariateFunction.\n\nInputs\n\nf - A UnivariateFunction.\n\nReturns\n\nA UnivariateFunction.\n\n\n\n\n\n","category":"function"},{"location":"api/#UnivariateFunctions.indefinite_integral","page":"API","title":"UnivariateFunctions.indefinite_integral","text":"indefinite_integral(f::UnivariateFunction)\n\nThis calculates the indefinite integral of a UnivariateFunction.\n\nInputs\n\nf - A UnivariateFunction.\n\nReturns\n\nA UnivariateFunction.\n\n\n\n\n\n","category":"function"},{"location":"api/#SchumakerSpline.evaluate_integral","page":"API","title":"SchumakerSpline.evaluate_integral","text":"evaluate_integral(f::UnivariateFunction,left::Real, right::Real)\nevaluate_integral(f::UnivariateFunction,left::Union{Date,DateTime}, right::Union{Date,DateTime})\n\nThis calculates the integral of a function from a left limit to a right limit and returns a scalar.\n\nIf a Date, DateTime is input then it is first converted to a scalar with the yearsfromglobal_base function.\n\nInputs\n\nf - A UnivariateFunction.\nleft - A left limit (scalar)\nright - A right limit (scalar)\n\nReturns\n\nA Real.\n\n\n\n\n\n","category":"function"},{"location":"api/#UnivariateFunctions.right_integral","page":"API","title":"UnivariateFunctions.right_integral","text":"right_integral(f::UnivariateFunction, left::Real)\nright_integral(f::UnivariateFunction, left::Union{Date,DateTime})\n\nThis calculates the integral of a function from a left limit and returns it as a UnivariateFunction. So if you were to then evaluate this integral function at a point x then you would get the integral between left and x.\n\nIf a Date, DateTime is input then it is first converted to a scalar with the yearsfromglobal_base function.\n\nInputs\n\nf - A UnivariateFunction.\nleft - A left limit (scalar)\n\nReturns\n\nA UnivariateFunction.\n\n\n\n\n\n","category":"function"},{"location":"api/#UnivariateFunctions.left_integral","page":"API","title":"UnivariateFunctions.left_integral","text":"left_integral(f::UnivariateFunction, right::Real)\nleft_integral(f::UnivariateFunction, right::Union{Date,DateTime})\n\nThis calculates the integral of a function from a right limit and returns it as a UnivariateFunction. So if you were to then evaluate this integral function at a point x then you would get the integral between right and x.\n\nIf a Date, DateTime is input then it is first converted to a scalar with the yearsfromglobal_base function.\n\nInputs\n\nf - A UnivariateFunction.\nright - A right limit (scalar)\n\nReturns\n\nA UnivariateFunction.\n\n\n\n\n\n","category":"function"},{"location":"api/#Interpolation","page":"API","title":"Interpolation","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"create_quadratic_spline\ncreate_constant_interpolation_to_right\ncreate_constant_interpolation_to_left\ncreate_linear_interpolation","category":"page"},{"location":"api/#UnivariateFunctions.create_quadratic_spline","page":"API","title":"UnivariateFunctions.create_quadratic_spline","text":"create_quadratic_spline(x::Union{Vector{DateTime},Vector{Date},Vector{Union{Date,DateTime}}},y::Vector{<:Real} ; gradients::Union{Missing,Vector{<:Real}} = missing, extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve),\n                             left_gradient::Union{Missing,Real} = missing, right_gradient::Union{Missing,Real} = missing)\ncreate_quadratic_spline(x::Vector{<:Real},y::Vector{<:Real} ; gradients::Union{Missing,Array{<:Real}} = missing, extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve),\n                             left_gradient::Union{Missing,Real} = missing, right_gradient::Union{Missing,Real} = missing)\ncreate_quadratic_spline(x::Union{Vector{D},Vector{<:DatePeriod}},y::Vector{<:Real}; gradients::Union{Missing,Vector{<:Real}} = missing, extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve),\n                             left_gradient::Union{Missing,Real} = missing, right_gradient::Union{Missing,Real} = missing) where D<:DatePeriod\n\nMakes a quadratic shape-preserving interpolation spline using the SchumakerSpline.jl package. This is returned as a Piecewise_Function rather than as a Schumaker struct.\n\nInputs\n\nx - A Vector with the x coordinates\ny - A Vector with the y coordinates\ngradients - A Vector with the gradiants at each x location. This is calculated if not provided.\nextrapolation - A tuple of enum value describing how to extrapolate (on the left and right sides).\nleft_gradient - The gradiant to impose on the left edge (ie the first x coordinate).\nright_gradient - The gradiant to impose on the right edge (ie the last x coordinate).\n\nReturns\n\nA Piecewise_Function containing the spline.\n\ncreate_quadratic_spline(schum::Schumaker)\n\nThis converts a spline represented by a SchumakerSpline.Schumaker struct into the same spline but represented by a Piecewise_Function.\n\nInputs\n\nschum - A Schumaker struct.\n\nReturns\n\nA Piecewise_Function containing the spline.\n\n\n\n\n\n","category":"function"},{"location":"api/#UnivariateFunctions.create_constant_interpolation_to_right","page":"API","title":"UnivariateFunctions.create_constant_interpolation_to_right","text":"create_constant_interpolation_to_right(x::Vector{Date},y::Vector{<:Real})\ncreate_constant_interpolation_to_right(x::Vector{<:Real},y::Vector{<:Real})\ncreate_constant_interpolation_to_right(x::Union{Vector{D},Vector{<:DatePeriod}},y::Vector{<:Real}) where D<:DatePeriod\n\nMakes a piecewise constant interpolation function. values from the left are copied to the right.\n\nInputs\n\nx - A Vector with the x coordinates\ny - A Vector with the y coordinates\n\nReturns\n\nA Piecewise_Function containing the interpolation function.\n\n\n\n\n\n","category":"function"},{"location":"api/#UnivariateFunctions.create_constant_interpolation_to_left","page":"API","title":"UnivariateFunctions.create_constant_interpolation_to_left","text":"create_constant_interpolation_to_left(x::Vector{Date},y::Vector{<:Real})\ncreate_constant_interpolation_to_left(x::Vector{<:Real},y::Vector{<:Real})\ncreate_constant_interpolation_to_left(x::Union{Vector{D},Vector{<:DatePeriod}},y::Vector{<:Real}) where D<:DatePeriod\n\nMakes a piecewise constant interpolation function. values from the right are copied to the left.\n\nInputs\n\nx - A Vector with the x coordinates\ny - A Vector with the y coordinates\n\nReturns\n\nA Piecewise_Function containing the interpolation function.\n\n\n\n\n\n","category":"function"},{"location":"api/#UnivariateFunctions.create_linear_interpolation","page":"API","title":"UnivariateFunctions.create_linear_interpolation","text":"create_linear_interpolation(x::Vector{Date},y::Vector{<:Real})\ncreate_linear_interpolation(x::Union{Vector{D},Vector{<:DatePeriod}},y::Vector{<:Real}) where D<:DatePeriod\ncreate_linear_interpolation(x::Vector{R},y::Vector{<:Real}) where R<:Real\n\nMakes a piecewise linear interpolation function. This is continuous.\n\nInputs\n\nx - A Vector with the x coordinates\ny - A Vector with the y coordinates\n\nReturns\n\nA Piecewise_Function containing the interpolation function.\n\n\n\n\n\n","category":"function"},{"location":"api/#Approximation","page":"API","title":"Approximation","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"create_ols_approximation\ncreate_chebyshev_approximation","category":"page"},{"location":"api/#UnivariateFunctions.create_ols_approximation","page":"API","title":"UnivariateFunctions.create_ols_approximation","text":"create_ols_approximation(y::Vector{<:Real}, x::Vector{<:Real}, base_x::Real = 0.0, degree::Integer = 1, intercept::Bool = true)\ncreate_ols_approximation(y::Vector{<:Real}, x::Union{Vector{DateTime},Vector{Date},Vector{Union{Date,DateTime}}}, base_x::Union{Date,DateTime} = global_base_date, degree::Integer = 1, intercept::Bool = true)\n\nAn approximation function calculated via OLS.\n\nInputs\n\ny - A Vector with the y coordinates\nx - A Vector with the x coordinates\nbase_x - A real that offsets the x. So a coordinate with x value of 2.0 will be converted to 1.8 if base_x is 0.2.\ndegree - What the highest power of x should be. So if this is 3 then the equation will have x, x^2, x^3 as predictors.\nintercept - Should there be an x intercept.\n\nReturns\n\nA Sum_Of_Functions containing the approximation function.\n\n\n\n\n\n","category":"function"},{"location":"api/#UnivariateFunctions.create_chebyshev_approximation","page":"API","title":"UnivariateFunctions.create_chebyshev_approximation","text":"create_chebyshev_approximation(func::Function, nodes::Integer, degree::Integer, left::Real, right::Real)\n\nAn function that will approximate another function via Chebyshev polynomials.\n\nInputs\n\nfunc - A function that you want to approximation\nnodes - The number of approximation nodes\ndegree - The degree of the Chebyshev polynomials.\nleft - The left limit of the approximation\nright - The right limit of the approximation.\n\nReturns\n\nA Sum_Of_Functions containing the approximation function.\n\n\n\n\n\n","category":"function"},{"location":"api/#Internal-Functions-2","page":"API","title":"Internal Functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"change_base_of_PE_Function\ntrim_piecewise_function\nsort\nconvert_to_linearly_rescale_inputs\nget_chevyshevs_up_to","category":"page"},{"location":"api/#UnivariateFunctions.change_base_of_PE_Function","page":"API","title":"UnivariateFunctions.change_base_of_PE_Function","text":"change_base_of_PE_Function(f::PE_Function, new_base::Real)\n\nThis changes the base of a PE_Function. So if the base was 2 then it can be converted to 3 with an additional constant term.\n\nInputs\n\nf - A PE_Function.\nnew_base - The new base.\n\nReturns\n\nA PE_Function or a Sum_Of_Functions.\n\n\n\n\n\n","category":"function"},{"location":"api/#UnivariateFunctions.trim_piecewise_function","page":"API","title":"UnivariateFunctions.trim_piecewise_function","text":"trim_piecewise_function(func::Piecewise_Function, left_limit::Real, right_limit::Real)\n\nThis trims the end of a piecewise function. So if there is a piecewise function with support between -10,10 then you can trim it to only have support between -5 and 5. Then if it is evaluated outside -5 to 5 it will be undefined.\n\nInputs\n\nfunc - A Piecewise_Function.\nleft_limit - The left limit.\nright_limit - The right limit.\n\nReturns\n\nA Piecewise_Function.\n\n\n\n\n\n","category":"function"},{"location":"api/#UnivariateFunctions.convert_to_linearly_rescale_inputs","page":"API","title":"UnivariateFunctions.convert_to_linearly_rescale_inputs","text":"convert_to_linearly_rescale_inputs(f::UnivariateFunction, alpha::Real, beta::Real)\n\nThis alters a function so that whenever we put in x it is like we put in alpha x + beta.\n\nInputs\n\nf - A UnivariateFunction.\nalpha - The slope of the rescaling.\nbeta - The level of the rescaling.\n\nReturns\n\nA UnivariateFunction of the type that you input to the function.\n\n\n\n\n\n","category":"function"},{"location":"api/#UnivariateFunctions.get_chevyshevs_up_to","page":"API","title":"UnivariateFunctions.get_chevyshevs_up_to","text":"get_chevyshevs_up_to(N::Integer, first_kind::Bool = true)\n\nGet the first N chebyshev polynomials returned as a vector of UnivariateFunctions. The first 20 polynomials of each are precompiled into the binaries for speed. If you need more than that they will be calculated at runtime.\n\nThese can be from either the first kind or second kind polynomial sequence.\n\nInputs\n\nN - How many chebyshev polynomials do you want.\nfirst_kind - A Bool. If true you get first kind polynomials. If false you get second kind.\n\nReturns\n\nA Vector of UnivariateFunctions for each polynomial.\n\n\n\n\n\n","category":"function"},{"location":"api/#Date-Conversions","page":"API","title":"Date Conversions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"years_between\nyears_from_global_base\nperiod_length","category":"page"},{"location":"api/#UnivariateFunctions.years_between","page":"API","title":"UnivariateFunctions.years_between","text":"years_between(a::Union{DateTime,Date}, b::Union{DateTime,Date})\n\nThe number of years between two dates. This is returned as a scalar and assumes 365.2422 days per year.\n\nInputs\n\na - The end date\nb - The start date.\n\nReturns\n\nA Real.\n\n\n\n\n\n","category":"function"},{"location":"api/#UnivariateFunctions.years_from_global_base","page":"API","title":"UnivariateFunctions.years_from_global_base","text":"years_from_global_base(a::Union{DateTime,Date})\n\nThe number of years (calculated by the years_between function) between a given date and the 1st of January 2000.\n\nInputs\n\ndate - A date.\n\nReturns\n\nA Real.\n\n\n\n\n\n","category":"function"},{"location":"api/#UnivariateFunctions.period_length","page":"API","title":"UnivariateFunctions.period_length","text":"period_length(a::Dates.DatePeriod, base::Date = global_base_date)\n\nPeriod length is designed to convert TimePeriod objects to a float in a consistent way to years_from_global_base. So effectively the years_between method is calculated with start and end dates being those at the start and end of a Dates.DatePeriod. This is slightly complicated because a period like Month(3) might have slightly different numbers of total days depending on when in the year it is. So a base date has to be input. The period is then measured starting from this base date.\n\nInputs\n\nperiod - A period.\nbase - A date from which the period will be measured from.\n\nReturns\n\nA Real.\n\n\n\n\n\n","category":"function"},{"location":"#Structs","page":"Index","title":"Structs","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"There are four main UnivariateFunction structs that are part of this package. These are:","category":"page"},{"location":"","page":"Index","title":"Index","text":"Undefined_Function - An undefined function behaves similarly to \"missing\" in Julia. Whenever anything is added/multiplied/etc with an undefined function the result is undefined. The integral and derivative of an undefined function is undefined. If an undefined function is evaluated it will return a missing.\nPEFunction - This is the basic function type. It has a form of a \\exp(b(x-base)) (x-base)^d$.\nSumOfFunctions - This is an array of PEFunctions. Note that by adding PEFunctions we can replicate any given polynomial. Hence from Weierstrass' approximation theorem we can approximate any continuous function on a bounded domain to any desired level of accuracy (whether this is practical in numerical computing depends on the function being approximated).\nPiecewise_Function - This defines a different UnivariateFunction for each part of the x domain.","category":"page"},{"location":"","page":"Index","title":"Index","text":"It is possible to perform any additions, subtractions, multiplications between any two UnivariateFunctions and between Ints/Floats and any UnivariateFunction. No division is allowed and it is not possible to raise a UnivariateFunction to a negative power. This is to ensure that all univariatefunctions are analytically integrable and differentiable. This may change in future releases.","category":"page"},{"location":"#Major-limitations","page":"Index","title":"Major limitations","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"It is not possible to divide by univariate functions or raise them by a negative power.\nWhen multiplying pefunctions with different base dates there is often an issue of very high or very low numbers that go outside machine precision. If one were trying to change a PEFunction from base 2010 to 50, this would not generally be possible. This is because to change a exp(x-2020) to q exp(x - 50) we need to premultiply the first expression by exp(-1950) which is a tiny number. In these cases it is better to do the algebra on paper and rewriting the code accordingly as often base changes cancel out on paper. It is also good to change bases as rarely as possible. If different univariate functions use different bases then there is a need to base change when multiplying them which can result in errors. Note that if base changes are segment in the x domain by means of a piecewise function then they should never interact meaning it is ok to use different bases here.\nThere is no support for finding optima, roots, fixedpoints etc. If anyone has an idea of how to do it efficiently then please let me know.\nThere is no support for finding a function representing the upper/lower envelope of multiple functions. If anyone has an idea of how to do it efficiently then please let me know.","category":"page"},{"location":"#Interpolation-and-Splines","page":"Index","title":"Interpolation and Splines","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"So far this package support the following interpolation schemes:","category":"page"},{"location":"","page":"Index","title":"Index","text":"Constant interpolation from the left to the right. Such a PiecewiseFunction spline can be constructed by the createconstantinterpolationto_right method.\nConstant interpolation from the right to the left. Such a PiecewiseFunction spline can be constructed by the createconstantinterpolationto_left method.\nLinear interpolation. Such a PiecewiseFunction spline can be constructed by the createlinear_interpolation method.","category":"page"},{"location":"","page":"Index","title":"Index","text":"It also supports the following spline (which can also be used for interpolation)","category":"page"},{"location":"","page":"Index","title":"Index","text":"Schumaker shape preserving spline - Such a PiecewiseFunction spline can be constructed by the createquadratic_spline method.","category":"page"},{"location":"#Approximation-and-regression","page":"Index","title":"Approximation and regression","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"So for this package supports the creation of the following approximation schemes:","category":"page"},{"location":"","page":"Index","title":"Index","text":"OLS regression. The createolsapproximation function can create a UnivariateFunction approximating a linear relationship. The degree input to this function can be used to specify the number of higher powers of x to be used in approximating y. For instance if the degree is two then y will be approximated as a linear combination of x and x^2 as well as an intercept (if the intercept boolean is true).\nChebyshev polynomials - This will approximate a function using the Chebyshev basis functions. This approximation function can then be integrated to accomplish Chebyshev–Gauss quadrature.","category":"page"},{"location":"#Date-Handling","page":"Index","title":"Date Handling","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"All base dates are immediately converted to floats and are not otherwise saved. Thus there is no difference between a PEFunction created with a base as a float and one created with the matching date. This is done to simplify the code. All date conversions is done by finding the year fractions between the date and the global base date of Date(2000,1,1). This particular global base date should not affect anything as long as it is consistent. It is relatively trivial to change it (in the dateconversions.jl file) and recompile however if desired.","category":"page"},{"location":"#Examples","page":"Index","title":"Examples","text":"","category":"section"},{"location":"#For-basic-algebra:","page":"Index","title":"For basic algebra:","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"Consider we have a two functions f and g and want to add them, multiply them by some other function h, then square it and finally integrate the result between 2.0 and 2.8. This can be done analytically with UnivariateFunctions:","category":"page"},{"location":"","page":"Index","title":"Index","text":"f = PE_Function(1.0, 2.0, 4.0, 5)\ng = PE_Function(1.3, 2.0, 4.3, 2)\nh = PE_Function(5.0, 2.2, 1.0,0)\nresult_of_operations = (h*(f+g))^2\nevaluate_integral(result_of_operations, 2.0, 2.8)","category":"page"},{"location":"#For-data-interpolation","page":"Index","title":"For data interpolation","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"Suppose we have want to approximate some function with some sampled points. First to generate some points","category":"page"},{"location":"","page":"Index","title":"Index","text":"using UnivariateFunctions\nconst global_base_date = Date(2000,1,1)\nStartDate = Date(2018, 7, 21)\nx = Array{Date}(undef, 1000)\nfor i in 1:1000\n    x[i] = StartDate +Dates.Day(2* (i-1))\nend\nfunction ff(x::Date)\n    days_between = years_from_global_base(x)\n    return log(days_between) + sqrt(days_between)\nend\ny = ff.(x)","category":"page"},{"location":"","page":"Index","title":"Index","text":"Now we can generate a UnivariateFunction that can be used to easily interpolate from the sampled points:","category":"page"},{"location":"","page":"Index","title":"Index","text":"func = create_quadratic_spline(x,y)","category":"page"},{"location":"","page":"Index","title":"Index","text":"And we can evaluate from this function and integrate it and differentiate it in the normal way:","category":"page"},{"location":"","page":"Index","title":"Index","text":"evaluate(func, Date(2020,1,1))\nevaluate.(Ref(func), [Date(2020,1,1), Date(2021,1,2)])\nevaluate(derivative(func), Date(2021,1,2))\nevaluate_integral(func, Date(2020,1,1), Date(2021,1,2))","category":"page"},{"location":"","page":"Index","title":"Index","text":"If we had wanted to interpolate instead with a constant method(from left or from right) or by linearly interpolating then we could have just generated func with a different method: createconstantinterpolationtoleft, createconstantinterpolationtoright or createlinearinterpolation.","category":"page"},{"location":"","page":"Index","title":"Index","text":"If we have lots of data that we want to summarise with OLS","category":"page"},{"location":"","page":"Index","title":"Index","text":"# Generating example data\nusing Random\nRandom.seed!(1)\nobs = 1000\nX = rand(obs)\ny = X .+ rand(Normal(),obs) .+ 7\n# And now making an approximation function\napproxFunction = create_ols_approximation(y, X, 0.0, 2, true)","category":"page"},{"location":"","page":"Index","title":"Index","text":"And if we want to approximate the sin function in the [2.3, 5.6] bound with 7 polynomial terms and 20 approximation nodes:","category":"page"},{"location":"","page":"Index","title":"Index","text":"chebyshevs = create_chebyshev_approximation(sin, 20, 7, 2.3, 5.6)","category":"page"},{"location":"","page":"Index","title":"Index","text":"We can integrate the above term in the normal way to achieve Gauss-Chebyshev quadrature:","category":"page"},{"location":"","page":"Index","title":"Index","text":"evaluate_integral(chebyshevs, 2.3, 5.6)","category":"page"}]
}
