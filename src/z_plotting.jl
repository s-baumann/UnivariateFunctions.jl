function expand_grid(vector_of_vectors, nombres::Vector{Symbol})
      if length(vector_of_vectors) != length(nombres) error("The two input vectors have different lengths.") end
      v = vec(collect(Base.product(  vector_of_vectors...  )))
      df = DataFrame()
      n = length(v)
      if n < 1 return df end
      for j in 1:length(v[1])
          df[!, nombres[j]] = Vector{typeof(v[1][j])}(undef,n)
      end
      for i in 1:n
            for j in 1:length(v[1])
                df[i, nombres[j]] = v[i][j]
            end
      end
      return df
end

function plot(fun::UnivariateFunction, left::Real = 0.0, right::Real = 1.0)
    X = convert(Array{Float64,1}, left:0.01:right)
    y = UnivariateFunctions.evaluate.(Ref(fun), X)
    df = DataFrame(x = X, y = y)
    plt = df |> @vlplot(mark = {:line, color = "blue"}, x = :x, y = :y)
    return plt
end
function plot(fun::UnivariateFunction, left::Real, right::Real, dd; width = 1000, height = 1000, x_name = :x, y_name = :y)
    X = convert(Array{Float64,1}, left:0.01:right)
    y = UnivariateFunctions.evaluate.(Ref(fun), X)
    df = DataFrame(x = X, y = y)
    plt_points = @vlplot(mark = {:point, color = "red", opacity = 0.5}, data = dd, encoding = {x = {x_name, type = "quantitative"}, y = {y_name, type = "quantitative"}})
    plt_line = @vlplot(mark = {:line, color = "blue"}, data = df, encoding = {x = {:x, type = "quantitative"}, y = {:y, type = "quantitative"}})
    return @vlplot(width=width,height=height) + plt_points + plt_line
end
function plot(fun::UnivariateFunction, dd; width = 1000, height = 1000, x_name = :x, y_name = :y)
    return plot(fun, minimum(dd[!, x_name]), maximum(dd[!, x_name]), dd; width=width, height=height, x_name=x_name, y_name=y_name)
end 

function plot(fun::Vector{<:UnivariateFunction}, left::Real = 0.0, right::Real = 1.0)
    ff = Dict{Int,UnivariateFunction}()
    for i in 1:length(fun)
        ff[i] = fun[i]
    end
    return plot(ff, left, right)
end

function plot(fun::Dict{ss,tt}, left::Real = 0.0, right::Real = 1.0) where ss<:Union{Symbol,Integer,String} where tt<:UnivariateFunction
    X = convert(Array{Float64,1}, left:0.01:right)
    labs = collect(keys(fun))
    dd = expand_grid(Vector[X, labs], [:x, :lab])
    dd[!, :y] = Vector{Float64}(undef, nrow(dd))
    for i in 1:nrow(dd)
        lab = dd[i, :lab]
        dd[i, :y] = UnivariateFunctions.evaluate(fun[lab], dd[i, :x])
    end
    plt = dd |> @vlplot(
        mark = :line,
        data = dd,
        encoding = {
            x = {"x", type = "quantitative"},
            y = {"y", type = "quantitative"},
            color = {"lab:n", title = "Label"}
        }
    )
    return plt
end

