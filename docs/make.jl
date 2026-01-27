using Documenter, UnivariateFunctions

makedocs(
    format = Documenter.HTML(),
    sitename = "UnivariateFunctions",
    modules = [UnivariateFunctions],
    pages = Any["Overview" => "index.md",
                "Examples" => "Examples.md",
                "Regression and Smoothing" => "Regression.md",
                "API" => "api.md"]
)

deploydocs(
    repo   = "github.com/s-baumann/UnivariateFunctions.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
