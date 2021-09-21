# UnivariateFunctions.jl

| Build | Coverage |
|-------|----------|



| Build | Coverage | Documentation |
|-------|----------|---------------|
| [![Build status](https://github.com/s-baumann/UnivariateFunctions.jl/workflows/CI/badge.svg)](https://github.com/s-baumann/UnivariateFunctions.jl/actions) | [![codecov](https://codecov.io/gh/s-baumann/UnivariateFunctions.jl/branch/master/graph/badge.svg?token=uO1mDGPfML)](https://codecov.io/gh/s-baumann/UnivariateFunctions.jl) | [![docs-latest-img](https://img.shields.io/badge/docs-latest-blue.svg)](https://s-baumann.github.io/UnivariateFunctions.jl/dev/index.html) |


This implements single algebra and evaluation on simple univariate functions.
There are a few ways in which it can be used.
* UnivariateFunctions can be used in the creation of splines. This has the added
    advantage that a spline implemented as a UnivariateFunction can be manipulated
    easily. It can be differentiated and then added to another function, etc.
* Any continuous interpolation scheme can be done with the added benefit that derivatives/integrals/products etc can be found analytically.
* Basic approximation schemes like OLS regression and chebyshev polynomials can be done with the added benefit that derivatives/integrals/products etc can be found analytically.

This is faster and simpler than the closely related package [MultivariateFunctions.jl](https://github.com/s-baumann/MultivariateFunctions.jl). The cost is that it is restricted to only one dimension.
