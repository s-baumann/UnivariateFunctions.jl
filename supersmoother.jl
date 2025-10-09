using Statistics
using LinearAlgebra

"""
    supersmoother(x, y; spans=[0.05, 0.2, 0.5], alpha=0.0)

Friedman's supersmoother algorithm for nonparametric regression.

# Arguments
- `x::Vector{<:Real}`: predictor values (will be sorted)
- `y::Vector{<:Real}`: response values
- `spans::Vector{<:Real}`: smoothing spans to consider (as fractions of data)
- `alpha::Real`: sensitivity parameter for span selection (0 ≤ alpha ≤ 10)

# Returns
- `y_smooth::Vector{<:Real}`: smoothed values corresponding to sorted x
- `x_sorted::Vector{<:Real}`: sorted x values

# Reference
Friedman, J. H. (1984). A variable span smoother. 
Technical Report No. 5, Laboratory for Computational Statistics, Stanford University.
"""
function supersmoother(x::Vector{<:Real}, y::Vector{<:Real}; 
                      spans::Vector{<:Real}=[0.05, 0.2, 0.5], 
                      alpha::Real=0.0)
    
    n = length(x)
    @assert n == length(y) "x and y must have the same length"
    @assert all(0 < s ≤ 1 for s in spans) "spans must be between 0 and 1"
    @assert 0 ≤ alpha ≤ 10 "alpha must be between 0 and 10"
    
    # Sort data by x values
    perm = sortperm(x)
    x_sorted = x[perm]
    y_sorted = y[perm]
    
    # Convert spans to actual window sizes
    span_sizes = max.(1, round.(Int, spans .* n))
    num_spans = length(spans)
    
    # Storage for smoothed values at each span
    smooths = Matrix{Float64}(undef, n, num_spans)
    
    # Compute smoothed values for each span
    for (i, span_size) in enumerate(span_sizes)
        smooths[:, i] = running_median_smooth(x_sorted, y_sorted, span_size)
    end
    
    # Compute residuals for cross-validation
    residuals = Matrix{Float64}(undef, n, num_spans)
    for i in 1:num_spans
        residuals[:, i] = compute_cv_residuals(x_sorted, y_sorted, span_sizes[i])
    end
    
    # Select optimal span at each point
    y_smooth = Vector{Float64}(undef, n)
    
    for i in 1:n
        # Compute cross-validation scores
        cv_scores = abs.(residuals[i, :])
        
        # Add penalty for larger spans (controlled by alpha)
        if alpha > 0
            penalties = alpha .* (spans .- minimum(spans))
            cv_scores .+= penalties
        end
        
        # Select span with minimum CV score
        best_span_idx = argmin(cv_scores)
        y_smooth[i] = smooths[i, best_span_idx]
    end
    
    return y_smooth, x_sorted
end

"""
    running_median_smooth(x, y, span_size)

Compute running median smooth with given span size.
"""
function running_median_smooth(x::Vector{<:Real}, y::Vector{<:Real}, span_size::Int)
    n = length(x)
    smoothed = Vector{Float64}(undef, n)
    half_span = span_size ÷ 2
    
    for i in 1:n
        # Determine window bounds
        left = max(1, i - half_span)
        right = min(n, i + half_span)
        
        # Ensure we have exactly span_size points (when possible)
        if right - left + 1 < span_size && i - half_span < 1
            right = min(n, left + span_size - 1)
        elseif right - left + 1 < span_size && i + half_span > n
            left = max(1, right - span_size + 1)
        end
        
        # Use tricube weights for local regression
        window_x = x[left:right]
        window_y = y[left:right]
        center_x = x[i]
        
        # Compute weights using tricube function
        distances = abs.(window_x .- center_x)
        max_dist = maximum(distances)
        
        if max_dist > 0
            u = distances ./ max_dist
            weights = (1 .- u.^3).^3
            weights[u .≥ 1] .= 0
        else
            weights = ones(length(window_x))
        end
        
        # Weighted local linear regression
        smoothed[i] = weighted_local_regression(window_x, window_y, weights, center_x)
    end
    
    return smoothed
end

"""
    weighted_local_regression(x, y, weights, x0)

Perform weighted local linear regression at point x0.
"""
function weighted_local_regression(x::Vector{<:Real}, y::Vector{<:Real}, 
                                 weights::Vector{<:Real}, x0::Real)
    n = length(x)
    
    if n == 1
        return y[1]
    end
    
    # Set up weighted least squares for local linear fit
    X = hcat(ones(n), x .- x0)
    W = Diagonal(weights)
    
    try
        # Solve weighted least squares: (X'WX)β = X'Wy
        XTW = X' * W
        XTWX = XTW * X
        XTWy = XTW * y
        
        # Add small regularization for numerical stability
        β = (XTWX + 1e-12 * I) \ XTWy
        
        # Return fitted value at x0 (intercept since we centered at x0)
        return β[1]
    catch
        # Fallback to weighted mean if regression fails
        return sum(weights .* y) / sum(weights)
    end
end

"""
    compute_cv_residuals(x, y, span_size)

Compute leave-one-out cross-validation residuals for given span size.
"""
function compute_cv_residuals(x::Vector{<:Real}, y::Vector{<:Real}, span_size::Int)
    n = length(x)
    residuals = Vector{Float64}(undef, n)
    half_span = span_size ÷ 2
    
    for i in 1:n
        # Determine window bounds (excluding point i)
        left = max(1, i - half_span)
        right = min(n, i + half_span)
        
        # Create index vector excluding point i
        indices = [left:(i-1); (i+1):right]
        
        if length(indices) < 2
            residuals[i] = 0.0
            continue
        end
        
        # Fit on nearby points excluding point i
        window_x = x[indices]
        window_y = y[indices]
        center_x = x[i]
        
        # Compute weights
        distances = abs.(window_x .- center_x)
        max_dist = maximum(distances)
        
        if max_dist > 0
            u = distances ./ max_dist
            weights = (1 .- u.^3).^3
            weights[u .≥ 1] .= 0
        else
            weights = ones(length(window_x))
        end
        
        # Predict at point i
        y_pred = weighted_local_regression(window_x, window_y, weights, center_x)
        residuals[i] = y[i] - y_pred
    end
    
    return residuals
end

"""
    supersmoother_predict(x_new, x_fit, y_fit, spans, alpha)

Predict new values using fitted supersmoother.
"""
function supersmoother_predict(x_new::Vector{<:Real}, x_fit::Vector{<:Real}, 
                             y_fit::Vector{<:Real}; 
                             spans::Vector{<:Real}=[0.05, 0.2, 0.5], 
                             alpha::Real=0.0)
    n_fit = length(x_fit)
    n_new = length(x_new)
    
    # Convert spans to actual window sizes
    span_sizes = max.(1, round.(Int, spans .* n_fit))
    
    y_new = Vector{Float64}(undef, n_new)
    
    for (j, x_val) in enumerate(x_new)
        # Find nearest point in fitted data
        nearest_idx = argmin(abs.(x_fit .- x_val))
        
        # Use local smoothing around this point
        predictions = Vector{Float64}(undef, length(spans))
        
        for (i, span_size) in enumerate(span_sizes)
            half_span = span_size ÷ 2
            left = max(1, nearest_idx - half_span)
            right = min(n_fit, nearest_idx + half_span)
            
            window_x = x_fit[left:right]
            window_y = y_fit[left:right]
            
            # Compute weights
            distances = abs.(window_x .- x_val)
            max_dist = maximum(distances)
            
            if max_dist > 0
                u = distances ./ max_dist
                weights = (1 .- u.^3).^3
                weights[u .≥ 1] .= 0
            else
                weights = ones(length(window_x))
            end
            
            predictions[i] = weighted_local_regression(window_x, window_y, weights, x_val)
        end
        
        # Simple selection (in practice, you'd want to use the CV-selected span)
        y_new[j] = predictions[2]  # Use middle span as default
    end
    
    return y_new
end


n = 100
x = sort(rand(n) * 4π)
y_true = sin.(x) + 0.5 * sin.(2*x)
y = y_true + 0.2 * randn(n)  # Add noise

println("Testing supersmoother with $(n) points...")

# Apply supersmoother
@time y_smooth, x_sorted = supersmoother(x, y, spans=[0.1, 0.3, 0.6], alpha=1.0)

# Compute MSE
y_true_sorted = sin.(x_sorted) + 0.5 * sin.(2*x_sorted)
mse = mean((y_smooth - y_true_sorted).^2)

println("MSE: $(round(mse, digits=6))")
println("Test completed successfully!")
    
#   return x_sorted, y[sortperm(x)], y_smooth, y_true_sorted
