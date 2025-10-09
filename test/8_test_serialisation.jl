using Test

@testset "Serialisation Tests" begin

    using LinearAlgebra, StatsBase
    using NonNegLeastSquares  # for nonneg_lsq
    using Plots
    using UnivariateFunctions
    using DataFrames, UUIDs

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
    
end


