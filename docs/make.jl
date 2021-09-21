using Documenter, UnivariateFunctions

makedocs(
    format = Documenter.HTML(),
    sitename = "UnivariateFunctions",
    modules = [UnivariateFunctions],
    pages = Any[
        "Introduction" => "description.md",
        "API" => "api.md"]
)

deploydocs(
    repo   = "github.com/s-baumann/UnivariateFunctions.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
