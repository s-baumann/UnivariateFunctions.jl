using Documenter, UnivariateFunctions

makedocs(
    format = Documenter.HTML(),
    sitename = "UnivariateFunctions",
    modules = [UnivariateFunctions],
    pages = Any["Index" => "index.md",
                "API" => "api.md"]
)

deploydocs(
    repo   = "github.com/s-baumann/UnivariateFunctions.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
