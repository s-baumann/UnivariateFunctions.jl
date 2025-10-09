
function isotonic_regression(x::Vector{R}, y::Vector{R}; increasing::Bool = true) where R<:Real
    if increasing == false
        y = -1.0 .* y
    end
    idx = sortperm(x)
    x = x[idx]
    y = y[idx]
    n = length(x)
    yhat = MultivariateStats.isotonic(x, y)

    starts_ = Vector{Float64}([-Inf])
    funcs_ = Vector{UnivariateFunction}([PE_Function(yhat[1])])

    last_x = x[1]
    last_y = yhat[1]
    now_x  = x[2]
    now_y  = yhat[2]
    constant_term_has_been_set = false

    gradient = (now_y - last_y) / (now_x - last_x)
    if (abs(gradient) < 1e-10)
        # The first two points are the same.
        push!(starts_, last_x)
        push!(funcs_, PE_Function(last_y))
        constant_term_has_been_set = true
    else
        # The first two points are different.
        push!(starts_, last_x)
        funn = PE_Function(last_y) + PE_Function(gradient, 0.0, last_x, 1)
        push!(funcs_, funn)
        constant_term_has_been_set = false
    end

    last_x = now_x
    last_y = now_y

    for i in 3:length(x)
        now_x = x[i]
        now_y = yhat[i]
        #
        gradient = (now_y - last_y) / (now_x - last_x)
        if (abs(gradient) < 1e-10) & (constant_term_has_been_set == false)
            # The first two points are the same.
            push!(starts_, last_x)
            push!(funcs_, PE_Function(last_y))
            constant_term_has_been_set = true
        elseif (gradient <= 1e-10)
            #
        else
            # The first two points are different.
            push!(starts_, last_x)
            funn = PE_Function(last_y) + PE_Function(gradient, 0.0, last_x, 1)
            push!(funcs_, funn)
            constant_term_has_been_set = false
        end
        # Update last_x and last_y
        last_x = now_x
        last_y = now_y
    end
    push!(starts_, x[end])
    push!(funcs_, PE_Function(yhat[end]))
    aa = Piecewise_Function(starts_, funcs_)
    if increasing == false
        aa = -1.0 .* aa
    end
    return aa
end

function isotonic_regression(dd::DataFrame, xvar::Symbol, yvar::Symbol; increasing::Bool = true)
    return isotonic_regression(dd[!, xvar], dd[!, yvar]; increasing=increasing)
end

function make_obs_grid(x::Vector{R}; nbins::Int=10) where R<:Real
    n = length(x)
    edges = Vector{R}()
    push!(edges, x[1])
    obs_spacing = (n/(nbins))
    for i in 0:(nbins-2)
        next_bit = Int(round((i+1)*obs_spacing))
        en = x[next_bit]
        push!(edges, en)
    end
    push!(edges, x[n])
    return edges
end


function monotonic_regression(x::Vector{R}, y::Vector{R}; nbins::Integer = 10, equally_spaced_bins::Bool=true, increasing::Bool = true) where R<:Real
    idx = sortperm(x)
    x = x[idx]
    y = y[idx]
    n = length(x)

    if increasing == false
        y = -1.0 .* y
    end

    edges = (equally_spaced_bins ? 
             collect(range(first(x), last(x), length=nbins+1)) :
             make_obs_grid(x; nbins=nbins)
    )

    m = length(edges) - 1
    Δ = diff(edges)

    A = zeros(n, m)
    for i in 1:n
        for j in 1:m
            if x[i] > edges[j+1]
                A[i,j] = Δ[j]
            elseif x[i] > edges[j]
                A[i,j] = x[i] - edges[j]
            else
                A[i,j] = 0.0
            end
        end
    end
    # Now we need to fit an intercept. The algorithm here does nonneg-least-squares.
    # So intercept coefficient needs to be positive. So we will do two constants one all 1 and one all -1 and it can choose.
    A = hcat(ones(n), -ones(n), A)
    
    # Solve nonnegative least squares for slopes
    ĝ = NonNegLeastSquares.nonneg_lsq(A, y)

    â = sum(ĝ[1] - ĝ[2])  # intercept is difference of two nonnegative coefficients

    starts_ = Vector{Float64}()
    funcs_ = Vector{UnivariateFunction}()
    start_y_accum = â
    for i in 1:(length(edges)-1)
        st = i == 1 ? -Inf : edges[i][1]
        gradient = ĝ[i+2,1]
        funn = PE_Function(start_y_accum) + PE_Function(gradient, 0.0, edges[i][1], 1)
        push!(starts_, st)
        push!(funcs_, funn)
        start_y_accum = start_y_accum + gradient * (i < length(edges) ? Δ[i] : 0.0)
    end
    aa = Piecewise_Function(starts_, funcs_)

    if increasing == false
        aa = -1.0 .* aa
    end

    return aa
end

function monotonic_regression(dd::DataFrame, xvar::Symbol, yvar::Symbol; nbins::Integer = 10, equally_spaced_bins::Bool=true, increasing::Bool = true)
    return monotonic_regression(dd[!, xvar], dd[!, yvar]; nbins=nbins, equally_spaced_bins=equally_spaced_bins, increasing=increasing)
end
