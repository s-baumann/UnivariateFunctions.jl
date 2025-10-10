function plot(fun::Union{UnivariateFunction,Undefined_Function, PE_Function, Sum_Of_Functions, Piecewise_Function}, left::Real = 0.0, right::Real = 1.0)
    X = convert(Array{Float64,1}, left:0.01:right)
    y = UnivariateFunctions.evaluate.(Ref(fun), X)
    df = DataFrame(x = X, y = y)
    plt = df |> @vlplot(:line, x = :x, y = :y, colour = "blue")
    return plt
end
function plot(fun::Union{UnivariateFunction,Undefined_Function, PE_Function, Sum_Of_Functions, Piecewise_Function}, left::Real, right::Real, dd; width = 1000, height = 1000, x_name = :x, y_name = :y)
    X = convert(Array{Float64,1}, left:0.01:right)
    y = UnivariateFunctions.evaluate.(Ref(fun), X)
    df = DataFrame(x = X, y = y)
    plt_points = @vlplot(mark = {:point, color = "red", opacity = 0.5}, data = dd, encoding = {x = {x_name, type = "quantitative"}, y = {y_name, type = "quantitative"}})
    plt_line = @vlplot(mark = {:line, color = "blue"}, data = df, encoding = {x = {x_name, type = "quantitative"}, y = {y_name, type = "quantitative"}})
    return @vlplot(width=width,height=height) + plt_points + plt_line
end
function plot(fun::Union{UnivariateFunction,Undefined_Function, PE_Function, Sum_Of_Functions, Piecewise_Function}, dd; width = 1000, height = 1000, x_name = :x, y_name = :y)
    return plot(fun, minimum(dd[!, x_name]), maximum(dd[!, x_name]), dd; width=width, height=height, x_name=x_name, y_name=y_name)
end 