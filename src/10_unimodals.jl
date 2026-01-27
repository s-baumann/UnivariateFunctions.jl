"""
    _kfold_indices(n::Int, nfolds::Int; seed::Union{Int,Nothing}=nothing)

Generate random fold assignments for k-fold CV.
Returns a vector of length n with values 1:nfolds.
"""
function _kfold_indices(n::Int, nfolds::Int; seed::Union{Int,Nothing}=nothing)
    if !isnothing(seed)
        Random.seed!(seed)
    end
    folds = repeat(1:nfolds, ceil(Int, n / nfolds))[1:n]
    return shuffle(folds)
end

"""
    _cv_error(x, y, fold_indices, fit_func; kwargs...)

Compute total cross-validation SSE for a given fitting function.
"""
function _cv_error(x::Vector{R}, y::Vector{R}, fold_indices::Vector{Int},
                   fit_func::Function; kwargs...) where R<:Real
    nfolds = maximum(fold_indices)
    total_sse = 0.0

    for k in 1:nfolds
        test_mask = fold_indices .== k
        train_mask = .!test_mask

        x_train = x[train_mask]
        y_train = y[train_mask]
        x_test = x[test_mask]
        y_test = y[test_mask]

        # Skip if train or test is empty
        if length(x_train) < 2 || length(x_test) < 1
            continue
        end

        # Fit on training data
        fitted = fit_func(x_train, y_train; kwargs...)

        # Evaluate on test data
        ŷ_test = fitted.(x_test)
        total_sse += sum((y_test .- ŷ_test).^2)
    end

    return total_sse
end

"""
    unimodal_regression(x, y; nbins=10, equally_spaced_bins=true, convex=false, quasi=true)

Fit a unimodal regression function to the data.

# Arguments
- `x`: Independent variable values
- `y`: Dependent variable values
- `nbins`: Number of bins for piecewise linear fit
- `equally_spaced_bins`: If true, bins are equally spaced; if false, based on observation quantiles
- `convex`: If false, fits concave/quasiconcave (single maximum); if true, fits convex/quasiconvex (single minimum)
- `quasi`: If true, only enforces unimodality (slopes change sign once); if false, also enforces curvature

# Shape constraints:
- `convex=false, quasi=true`:  Quasiconcave - slopes go from + to - (single peak)
- `convex=true,  quasi=true`:  Quasiconvex  - slopes go from - to + (single trough)
- `convex=false, quasi=false`: Concave - slopes monotonically decrease
- `convex=true,  quasi=false`: Convex  - slopes monotonically increase

Returns a `Piecewise_Function`.
"""
function unimodal_regression(x::Vector{R}, y::Vector{R};
                              nbins::Integer = 10,
                              equally_spaced_bins::Bool = true,
                              convex::Bool = false,
                              quasi::Bool = true) where R<:Real
    idx = sortperm(x)
    x = x[idx]
    y = y[idx]
    n = length(x)

    # For convex (single minimum), we flip y, fit quasiconcave, then flip back
    if convex
        y = -1.0 .* y
    end

    edges = (equally_spaced_bins ?
             collect(range(first(x), last(x), length=nbins+1)) :
             make_obs_grid(x; nbins=nbins))

    m = length(edges) - 1
    Δ = diff(edges)

    if quasi
        # Quasiconcave: Find optimal peak location, fit increasing then decreasing
        best_sse = Inf
        best_result = nothing

        # Try each possible peak location (bin boundary)
        for peak_bin in 1:m
            # Build design matrix:
            # - Bins 1:peak_bin have positive slopes (increasing part)
            # - Bins (peak_bin+1):m have negative slopes (decreasing part)

            A = zeros(n, m + 2)  # +2 for positive and negative intercept
            A[:, 1] .= 1.0       # positive intercept
            A[:, 2] .= -1.0      # negative intercept

            for i in 1:n
                for j in 1:m
                    if x[i] > edges[j+1]
                        contribution = Δ[j]
                    elseif x[i] > edges[j]
                        contribution = x[i] - edges[j]
                    else
                        contribution = 0.0
                    end

                    if j <= peak_bin
                        # Increasing part: positive slope
                        A[i, j + 2] = contribution
                    else
                        # Decreasing part: negative slope (we negate so NNLS gives positive coef)
                        A[i, j + 2] = -contribution
                    end
                end
            end

            # Solve nonnegative least squares
            ĝ = NonNegLeastSquares.nonneg_lsq(A, y)
            ŷ = A * ĝ
            sse = sum((y .- ŷ).^2)

            if sse < best_sse
                best_sse = sse
                best_result = (ĝ = ĝ, peak_bin = peak_bin, edges = edges, Δ = Δ)
            end
        end

        # Build piecewise function from best result
        ĝ = best_result.ĝ
        peak_bin = best_result.peak_bin

        â = ĝ[1] - ĝ[2]  # intercept

        starts_ = Vector{Float64}()
        funcs_ = Vector{UnivariateFunction}()
        start_y_accum = â

        for i in 1:m
            st = i == 1 ? -Inf : edges[i]

            if i <= peak_bin
                gradient = ĝ[i + 2]  # positive slope
            else
                gradient = -ĝ[i + 2]  # negative slope
            end

            funn = PE_Function(start_y_accum) + PE_Function(gradient, 0.0, edges[i], 1)
            push!(starts_, st)
            push!(funcs_, funn)
            start_y_accum = start_y_accum + gradient * Δ[i]
        end

        aa = Piecewise_Function(starts_, funcs_)

    else
        # True concave: slopes must be monotonically decreasing
        # We parameterize as: g_1 >= g_2 >= ... >= g_m
        # Let g_j = h_1 + h_2 + ... + h_j where h_k >= 0 for the "excess" slope
        # Actually easier: g_j = s - (d_1 + d_2 + ... + d_{j-1}) where d_k >= 0
        # So g_1 = s, g_2 = s - d_1, g_3 = s - d_1 - d_2, etc.
        # We need s to be free (can be negative), so use s = s+ - s-

        # Design matrix columns: [1, -1, x_contrib_for_s, -x_contrib_for_d1, -x_contrib_for_d2, ...]
        # where x_contrib_for_s affects all bins, and d_j affects bins j+1 onwards

        A = zeros(n, m + 3)  # intercept+, intercept-, s+, s-, d_1, d_2, ..., d_{m-1}
        # Actually let's restructure: [int+, int-, s+, s-, d_1, ..., d_{m-1}]
        # That's 2 + 2 + (m-1) = m + 3 columns

        A = zeros(n, 2 + 2 + (m - 1))
        A[:, 1] .= 1.0   # positive intercept
        A[:, 2] .= -1.0  # negative intercept

        for i in 1:n
            for j in 1:m
                if x[i] > edges[j+1]
                    contribution = Δ[j]
                elseif x[i] > edges[j]
                    contribution = x[i] - edges[j]
                else
                    contribution = 0.0
                end

                # s+ contribution (positive initial slope)
                A[i, 3] += contribution
                # s- contribution (negative initial slope)
                A[i, 4] -= contribution

                # d_k contributions: d_k decreases slope for bins k+1 onwards
                for k in 1:(j-1)
                    if k <= m - 1
                        A[i, 4 + k] -= contribution
                    end
                end
            end
        end

        ĝ = NonNegLeastSquares.nonneg_lsq(A, y)

        â = ĝ[1] - ĝ[2]
        s = ĝ[3] - ĝ[4]

        # Compute slopes: g_j = s - sum(d_1:d_{j-1})
        slopes = zeros(m)
        slopes[1] = s
        for j in 2:m
            slopes[j] = slopes[j-1] - ĝ[4 + j - 1]
        end

        starts_ = Vector{Float64}()
        funcs_ = Vector{UnivariateFunction}()
        start_y_accum = â

        for i in 1:m
            st = i == 1 ? -Inf : edges[i]
            gradient = slopes[i]
            funn = PE_Function(start_y_accum) + PE_Function(gradient, 0.0, edges[i], 1)
            push!(starts_, st)
            push!(funcs_, funn)
            start_y_accum = start_y_accum + gradient * Δ[i]
        end

        aa = Piecewise_Function(starts_, funcs_)
    end

    if convex
        aa = -1.0 .* aa
    end

    return aa
end

function unimodal_regression(dd::DataFrame, xvar::Symbol, yvar::Symbol; kwargs...)
    return unimodal_regression(dd[!, xvar], dd[!, yvar]; kwargs...)
end

"""
    CVRegressionResult{F<:UnivariateFunction}

Result type for cross-validated regression model selection.

# Fields
- `fitted::F`: The fitted `UnivariateFunction` (typically a `Piecewise_Function`)
- `selected_shape::Symbol`: The shape that was selected (e.g., `:increasing`, `:quasiconcave`)
- `cv_errors::Dict{Symbol, Float64}`: Cross-validation errors for each candidate shape
- `nfolds::Int`: Number of folds used in cross-validation

The result is callable - you can use it directly as a function:
```julia
result = cv_shape_regression(x, y)
result(2.5)  # equivalent to result.fitted(2.5)
```
"""
struct CVRegressionResult{F<:UnivariateFunction}
    fitted::F
    selected_shape::Symbol
    cv_errors::Dict{Symbol, Float64}
    nfolds::Int
end

# Make it callable like the function itself
(r::CVRegressionResult)(x) = r.fitted(x)

"""
    cv_monotonic_regression(x, y; nbins=10, equally_spaced_bins=true, nfolds=10, seed=nothing)

Fit monotonic regression, automatically selecting increasing vs decreasing
based on k-fold cross-validation error.

Returns a `CVRegressionResult` containing:
- `fitted`: The fitted function (using full dataset)
- `selected_shape`: Either `:increasing` or `:decreasing`
- `cv_errors`: Dict mapping each shape to its CV error
- `nfolds`: Number of folds used
"""
function cv_monotonic_regression(x::Vector{R}, y::Vector{R};
                                  nbins::Integer = 10,
                                  equally_spaced_bins::Bool = true,
                                  nfolds::Integer = 10,
                                  seed::Union{Int,Nothing} = nothing) where R<:Real

    fold_indices = _kfold_indices(length(x), nfolds; seed=seed)

    candidates = [
        (increasing = true,  key = :increasing),
        (increasing = false, key = :decreasing),
    ]

    cv_errors = Dict{Symbol, Float64}()
    best_error = Inf
    best_candidate = candidates[1]

    for candidate in candidates
        cv_err = _cv_error(x, y, fold_indices, monotonic_regression;
                           nbins=nbins,
                           equally_spaced_bins=equally_spaced_bins,
                           increasing=candidate.increasing)

        cv_errors[candidate.key] = cv_err

        if cv_err < best_error
            best_error = cv_err
            best_candidate = candidate
        end
    end

    # Fit final model on full data
    fitted = monotonic_regression(x, y;
                                   nbins=nbins,
                                   equally_spaced_bins=equally_spaced_bins,
                                   increasing=best_candidate.increasing)

    return CVRegressionResult(fitted, best_candidate.key, cv_errors, nfolds)
end

function cv_monotonic_regression(dd::DataFrame, xvar::Symbol, yvar::Symbol; kwargs...)
    return cv_monotonic_regression(dd[!, xvar], dd[!, yvar]; kwargs...)
end

"""
    cv_unimodal_regression(x, y; nbins=10, equally_spaced_bins=true, nfolds=10, seed=nothing)

Fit unimodal regression, automatically selecting among:
- `:quasiconcave` (convex=false, quasi=true)
- `:quasiconvex`  (convex=true,  quasi=true)
- `:concave`      (convex=false, quasi=false)
- `:convex`       (convex=true,  quasi=false)

based on k-fold cross-validation error.

Returns a `CVRegressionResult` containing:
- `fitted`: The fitted function (using full dataset)
- `selected_shape`: One of the four shape symbols
- `cv_errors`: Dict mapping each shape to its CV error
- `nfolds`: Number of folds used
"""
function cv_unimodal_regression(x::Vector{R}, y::Vector{R};
                                 nbins::Integer = 10,
                                 equally_spaced_bins::Bool = true,
                                 nfolds::Integer = 10,
                                 seed::Union{Int,Nothing} = nothing) where R<:Real

    fold_indices = _kfold_indices(length(x), nfolds; seed=seed)

    candidates = [
        (convex = false, quasi = true,  key = :quasiconcave),
        (convex = true,  quasi = true,  key = :quasiconvex),
        (convex = false, quasi = false, key = :concave),
        (convex = true,  quasi = false, key = :convex),
    ]

    cv_errors = Dict{Symbol, Float64}()
    best_error = Inf
    best_candidate = candidates[1]

    for candidate in candidates
        cv_err = _cv_error(x, y, fold_indices, unimodal_regression;
                           nbins=nbins,
                           equally_spaced_bins=equally_spaced_bins,
                           convex=candidate.convex,
                           quasi=candidate.quasi)

        cv_errors[candidate.key] = cv_err

        if cv_err < best_error
            best_error = cv_err
            best_candidate = candidate
        end
    end

    # Fit final model on full data
    fitted = unimodal_regression(x, y;
                                  nbins=nbins,
                                  equally_spaced_bins=equally_spaced_bins,
                                  convex=best_candidate.convex,
                                  quasi=best_candidate.quasi)

    return CVRegressionResult(fitted, best_candidate.key, cv_errors, nfolds)
end

function cv_unimodal_regression(dd::DataFrame, xvar::Symbol, yvar::Symbol; kwargs...)
    return cv_unimodal_regression(dd[!, xvar], dd[!, yvar]; kwargs...)
end

"""
    cv_shape_regression(x, y; shapes=:all, nbins=10, equally_spaced_bins=true, nfolds=10, seed=nothing)

Fit regression with automatic shape selection via cross-validation.

`shapes` can be:
- `:monotonic` - choose between increasing/decreasing
- `:unimodal`  - choose between quasiconcave/quasiconvex/concave/convex
- `:all`       - choose from all 6 shapes
- A vector of symbols, e.g. `[:increasing, :quasiconcave, :convex]`

Available shapes: `:increasing`, `:decreasing`, `:quasiconcave`, `:quasiconvex`, `:concave`, `:convex`

Returns a `CVRegressionResult` with the fitted function and selection metadata.
"""
function cv_shape_regression(x::Vector{R}, y::Vector{R};
                              shapes::Union{Symbol, Vector{Symbol}} = :all,
                              nbins::Integer = 10,
                              equally_spaced_bins::Bool = true,
                              nfolds::Integer = 10,
                              seed::Union{Int,Nothing} = nothing) where R<:Real

    all_candidates = Dict(
        :increasing   => (func = monotonic_regression, kwargs = (increasing = true,)),
        :decreasing   => (func = monotonic_regression, kwargs = (increasing = false,)),
        :quasiconcave => (func = unimodal_regression,  kwargs = (convex = false, quasi = true)),
        :quasiconvex  => (func = unimodal_regression,  kwargs = (convex = true,  quasi = true)),
        :concave      => (func = unimodal_regression,  kwargs = (convex = false, quasi = false)),
        :convex       => (func = unimodal_regression,  kwargs = (convex = true,  quasi = false)),
    )

    # Determine which candidates to try
    if shapes == :all
        candidate_keys = [:increasing, :decreasing, :quasiconcave, :quasiconvex, :concave, :convex]
    elseif shapes == :monotonic
        candidate_keys = [:increasing, :decreasing]
    elseif shapes == :unimodal
        candidate_keys = [:quasiconcave, :quasiconvex, :concave, :convex]
    else
        candidate_keys = shapes
    end

    fold_indices = _kfold_indices(length(x), nfolds; seed=seed)

    cv_errors = Dict{Symbol, Float64}()
    best_error = Inf
    best_key = candidate_keys[1]

    for key in candidate_keys
        candidate = all_candidates[key]
        cv_err = _cv_error(x, y, fold_indices, candidate.func;
                           nbins=nbins,
                           equally_spaced_bins=equally_spaced_bins,
                           candidate.kwargs...)

        cv_errors[key] = cv_err

        if cv_err < best_error
            best_error = cv_err
            best_key = key
        end
    end

    # Fit final model on full data
    best = all_candidates[best_key]
    fitted = best.func(x, y;
                        nbins=nbins,
                        equally_spaced_bins=equally_spaced_bins,
                        best.kwargs...)

    return CVRegressionResult(fitted, best_key, cv_errors, nfolds)
end

function cv_shape_regression(dd::DataFrame, xvar::Symbol, yvar::Symbol; kwargs...)
    return cv_shape_regression(dd[!, xvar], dd[!, yvar]; kwargs...)
end
