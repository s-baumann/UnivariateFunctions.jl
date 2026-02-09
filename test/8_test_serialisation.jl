using Test

@testset "Serialisation Tests" begin
    using UnivariateFunctions
    using DataFrames

    x = sort(rand(50)*10)
    function ff(x::Real)
        if x < 5.0
            return 2 * x + 0.5 * x^1.5 + randn(1)[1]*3
        else
            return 1 * x + 0.5 * x^2 + randn(1)[1]*5 + 5
        end
    end
    y = -1 .* ff.(x) .- 10

    uf = Undefined_Function()
    fit = UnivariateFunctions.monotonic_regression(x, y; nbins = 11, equally_spaced_bins=true, increasing=false)

    aa = DataFrames.DataFrame(fit)
    bb = UnivariateFunction(aa)
    cc = bb - fit
    vals = cc.(x)
    @test maximum(abs.(vals)) < 1e-10

    uf2 = UnivariateFunction(DataFrames.DataFrame(uf))
    @test isa(uf2, Undefined_Function)

    # Round-trip test for mixed Piecewise_Function (PE + Sum segments)
    pe1 = PE_Function(3.0, 0.0, 0.0, 0)  # constant 3
    sum1 = PE_Function(1.0, 0.0, 0.0, 1) + PE_Function(2.0, 0.0, 0.0, 0)  # x + 2
    pw_mixed = Piecewise_Function([0.0, 5.0], [pe1, sum1])
    df_mixed = DataFrames.DataFrame(pw_mixed)
    pw_roundtrip = UnivariateFunction(df_mixed)
    test_pts = [0.0, 2.5, 5.0, 7.5, 10.0]
    for p in test_pts
        @test abs(evaluate(pw_mixed, p) - evaluate(pw_roundtrip, p)) < 1e-10
    end

end


