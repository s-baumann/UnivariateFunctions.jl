using UnivariateFunctions
using Test

# Run tests

println("Test of UnivariateFunctions package.")
include("1_basic_tests.jl")
include("2_date_tests.jl")
include("3_piecewise_tests.jl")
include("4_schumaker_test.jl")
include("5_interpolation_test.jl")
include("6_test_regressions.jl")
