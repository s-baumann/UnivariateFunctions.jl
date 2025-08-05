function isotonic_regression(dd::DataFrame, xvar::Symbol, yvar::Symbol)
    sort!(dd, [xvar])
    dd[!, :yhat] = isotonic(dd[!, xvar], dd[!, yvar])

    starts_ = Vector{Float64}([-Inf])
    funcs_ = Vector{UnivariateFunction}([PE_Function(dd[1, :yhat])])

    last_x = dd[1, xvar]
    last_y = dd[1, :yhat]
    now_x  = dd[2, xvar]
    now_y  = dd[2, :yhat]
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

    for i in 3:nrow(dd)
        now_x = dd[i, :x]
        now_y = dd[i, :yhat]
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
        #
        last_x = now_x
        last_y = now_y
    end
    push!(starts_, dd[nrow(dd), xvar])
    push!(funcs_, PE_Function(dd[nrow(dd), :yhat]))
    aa = Piecewise_Function(starts_, funcs_)
    return aa
end
