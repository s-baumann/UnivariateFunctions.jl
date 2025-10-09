function DataFrames.DataFrame(uf::Undefined_Function)
    return DataFrame(:type => "undefined_function", :a => NaN, :b => NaN, :base => NaN, :d => NaN)
end
function DataFrames.DataFrame(uf::PE_Function)
    return DataFrame(:type => "pe_function", :a => uf.a_, :b => uf.b_, :base => uf.base_, :d => uf.d_)
end

function DataFrames.DataFrame(uf::Sum_Of_Functions)
    dfs = [DataFrames.DataFrame(f) for f in uf.functions_]
    df = vcat(dfs...)
    df[!, :group] .= UUIDs.uuid4()
    return df
end

function DataFrames.DataFrame(uf::Piecewise_Function)
    ll = []
    for i in 1:length(uf.starts_)
        s = uf.starts_[i]
        subdd = DataFrames.DataFrame(uf.functions_[i])
        subdd[!, :segment_start] .= s
        push!(ll, subdd)
    end
    return vcat(ll...)
end


function PE_Function(dd::Union{DataFrame,DataFrameRow})
    if isnan(dd.a[1]) | isnan(dd.b[1]) | isnan(dd.base[1]) | isnan(dd.d[1])
        return Undefined_Function()
    end
    return PE_Function(dd.a[1], dd.b[1], dd.base[1], dd.d[1])
end

function Sum_Of_Functions(dd::DataFrame)
    funs = [PE_Function(subdd) for subdd in eachrow(dd)]
    return Sum_Of_Functions(funs)
end

function Piecewise_Function(dd::DataFrame)
    grps = unique(dd.group)
    starts = Vector{Float64}()
    funs = []
    for g in grps
        subdd = dd[dd.group .== g, :]
        fun = Sum_Of_Functions(subdd)
        push!(starts, subdd.segment_start[1])
        push!(funs, fun)
    end
    idx = sortperm(starts)
    starts = starts[idx]
    funs = funs[idx]
    return Piecewise_Function(starts, funs)
end

function UnivariateFunction(df::DataFrame)
    noms = Symbol.(names(df))
    if :segment_start in noms
        return Piecewise_Function(df)
    elseif :group in noms
        return Sum_Of_Functions(df)
    elseif (:a in noms) & (:b in noms) & (:base in noms) & (:d in noms) & (nrow(df) == 1)
        return PE_Function(df)
    else
        error("DataFrame does not have the right columns to convert to a UnivariateFunction. You dont have :segment_start so its not a Piecewise_Function, you dont have :group so its not a Sum_Of_Functions, and you dont have the right columns to be a PE_Function or you have more than one row.")
    end
end
