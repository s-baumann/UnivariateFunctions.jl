using UnivariateFunctions
using Test

# Run tests

println("Test of UnivariateFunctions package.")
@time @test include("basic_tests.jl")
@time @test include("date_tests.jl")
@time @test include("piecewise_tests.jl")
@time @test include("schumaker_test.jl")
@time @test include("interpolation_test.jl")
