"""
    _local_linear_smooth(x::Vector{R}, y::Vector{R}, span::R) where R<:Real

Compute local linear regression smoothed values for all points.
`span` is the fraction of data points to use in each local fit (0 < span ≤ 1).
"""
function _local_linear_smooth(x::Vector{R}, y::Vector{R}, span::R) where R<:Real
    n = length(x)
    k = max(2, Int(ceil(span * n)))  # number of neighbors to use
    ŷ = zeros(n)
    
    for i in 1:n
        # Find k nearest neighbors by x distance
        distances = abs.(x .- x[i])
        neighbor_idx = sortperm(distances)[1:k]
        
        x_local = x[neighbor_idx]
        y_local = y[neighbor_idx]
        
        # Tricube weights based on distance
        max_dist = maximum(abs.(x_local .- x[i])) + 1e-10
        u = abs.(x_local .- x[i]) ./ max_dist
        w = (1.0 .- u.^3).^3  # tricube kernel
        
        # Weighted linear regression: ŷ = a + b*x
        sum_w = sum(w)
        sum_wx = sum(w .* x_local)
        sum_wy = sum(w .* y_local)
        sum_wxx = sum(w .* x_local.^2)
        sum_wxy = sum(w .* x_local .* y_local)
        
        denom = sum_w * sum_wxx - sum_wx^2
        if abs(denom) < 1e-12
            # Fallback to weighted mean
            ŷ[i] = sum_wy / sum_w
        else
            b = (sum_w * sum_wxy - sum_wx * sum_wy) / denom
            a = (sum_wy - b * sum_wx) / sum_w
            ŷ[i] = a + b * x[i]
        end
    end
    
    return ŷ
end

"""
    _local_linear_loo_residuals(x::Vector{R}, y::Vector{R}, span::R) where R<:Real

Compute leave-one-out residuals for local linear regression at given span.
Returns |yᵢ - ŷ₋ᵢ| for each point (prediction without using point i).
"""
function _local_linear_loo_residuals(x::Vector{R}, y::Vector{R}, span::R) where R<:Real
    n = length(x)
    k = max(3, Int(ceil(span * n)))  # need at least 3 so LOO still has 2
    residuals = zeros(n)
    
    for i in 1:n
        # Find k nearest neighbors, excluding self for LOO
        distances = abs.(x .- x[i])
        distances[i] = Inf  # exclude self
        neighbor_idx = sortperm(distances)[1:(k-1)]
        
        x_local = x[neighbor_idx]
        y_local = y[neighbor_idx]
        
        # Tricube weights
        max_dist = maximum(abs.(x_local .- x[i])) + 1e-10
        u = abs.(x_local .- x[i]) ./ max_dist
        w = (1.0 .- u.^3).^3
        
        # Weighted linear regression
        sum_w = sum(w)
        sum_wx = sum(w .* x_local)
        sum_wy = sum(w .* y_local)
        sum_wxx = sum(w .* x_local.^2)
        sum_wxy = sum(w .* x_local .* y_local)
        
        denom = sum_w * sum_wxx - sum_wx^2
        if abs(denom) < 1e-12
            ŷ_i = sum_wy / sum_w
        else
            b = (sum_w * sum_wxy - sum_wx * sum_wy) / denom
            a = (sum_wy - b * sum_wx) / sum_w
            ŷ_i = a + b * x[i]
        end
        
        residuals[i] = abs(y[i] - ŷ_i)
    end
    
    return residuals
end

"""
    _smooth_values(x::Vector{R}, v::Vector{R}, span::R) where R<:Real

Smooth a vector v using local mean with given span. Used to smooth the span selection.
"""
function _smooth_values(x::Vector{R}, v::Vector{R}, span::R) where R<:Real
    n = length(x)
    k = max(1, Int(ceil(span * n)))
    v_smooth = zeros(n)
    
    for i in 1:n
        distances = abs.(x .- x[i])
        neighbor_idx = sortperm(distances)[1:k]
        v_smooth[i] = mean(v[neighbor_idx])
    end
    
    return v_smooth
end

"""
    _supersmoother_values(x::Vector{R}, y::Vector{R}; 
                          spans::Vector{R} = [0.05, 0.2, 0.5],
                          bass::R = 0.0) where R<:Real

Core SuperSmoother algorithm. Returns smoothed y values.

- `spans`: The three candidate spans to try (fraction of data)
- `bass`: Bass enhancement parameter (0-10). Higher values produce smoother results.
"""
function _supersmoother_values(x::Vector{R}, y::Vector{R}; 
                                spans::Vector{R} = R[0.05, 0.2, 0.5],
                                bass::R = R(0.0)) where R<:Real
    n = length(x)
    n_spans = length(spans)
    
    # Step 1: Compute smoothed values at each candidate span
    smooths = [_local_linear_smooth(x, y, s) for s in spans]
    
    # Step 2: Compute LOO residuals at each span
    residuals = [_local_linear_loo_residuals(x, y, s) for s in spans]
    
    # Step 3: Smooth the residuals (to stabilize span selection)
    mid_span = spans[div(n_spans + 1, 2)]  # use middle span for smoothing
    smoothed_residuals = [_smooth_values(x, r, mid_span) for r in residuals]
    
    # Step 4: Apply bass enhancement (increases preference for larger spans)
    if bass > 0
        for j in 1:n_spans
            smoothed_residuals[j] .= smoothed_residuals[j] .* (spans[j] / spans[1])^(-bass / 10)
        end
    end
    
    # Step 5: At each point, pick the span with lowest smoothed residual
    best_span_idx = zeros(Int, n)
    for i in 1:n
        best_j = 1
        best_res = smoothed_residuals[1][i]
        for j in 2:n_spans
            if smoothed_residuals[j][i] < best_res
                best_res = smoothed_residuals[j][i]
                best_j = j
            end
        end
        best_span_idx[i] = best_j
    end
    
    # Step 6: Smooth the span selection to avoid erratic switching
    # Convert to continuous span values, smooth, then use for interpolation
    selected_spans = [spans[best_span_idx[i]] for i in 1:n]
    smoothed_spans = _smooth_values(x, selected_spans, mid_span)
    
    # Step 7: Compute final smoothed values by interpolating between span results
    ŷ = zeros(n)
    for i in 1:n
        s = smoothed_spans[i]
        
        # Find which two spans we're between and interpolate
        if s <= spans[1]
            ŷ[i] = smooths[1][i]
        elseif s >= spans[end]
            ŷ[i] = smooths[end][i]
        else
            # Find bracketing spans
            for j in 1:(n_spans-1)
                if spans[j] <= s <= spans[j+1]
                    t = (s - spans[j]) / (spans[j+1] - spans[j])
                    ŷ[i] = (1 - t) * smooths[j][i] + t * smooths[j+1][i]
                    break
                end
            end
        end
    end
    
    return ŷ
end

"""
    supersmoother(x::Vector{R}, y::Vector{R}; 
                  spans::Vector{R} = [0.05, 0.2, 0.5],
                  bass::R = 0.0) where R<:Real

Friedman's SuperSmoother (1984) - a local linear regression with adaptive bandwidth.

Returns a `Piecewise_Function` (linear interpolation of smoothed values).

# Arguments
- `x`: Independent variable values
- `y`: Dependent variable values  
- `spans`: Candidate spans to consider (fraction of data points). Default `[0.05, 0.2, 0.5]`
- `bass`: Bass enhancement (0-10). Higher values favor smoother fits. Default `0.0`

# Example
```julia
f = supersmoother(x, y)
f(2.5)  # evaluate at new point
derivative(f)  # get derivative function
```
"""
function supersmoother(x::Vector{R}, y::Vector{R}; 
                        spans::Vector{R} = R[0.05, 0.2, 0.5],
                        bass::R = R(0.0)) where R<:Real
    # Sort by x
    idx = sortperm(x)
    x_sorted = x[idx]
    y_sorted = y[idx]
    
    # Get smoothed values
    ŷ = _supersmoother_values(x_sorted, y_sorted; spans=spans, bass=bass)
    
    # Return as piecewise linear interpolation
    return create_linear_interpolation(x_sorted, ŷ)
end

function supersmoother(dd::DataFrame, xvar::Symbol, yvar::Symbol; kwargs...)
    return supersmoother(dd[!, xvar], dd[!, yvar]; kwargs...)
end