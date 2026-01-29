# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

UnivariateFunctions.jl is a Julia package for single-variable function algebra and evaluation. It enables creating, manipulating, and analyzing univariate functions analytically through differentiation, integration, and algebraic operations.

## Common Commands

```bash
# Run all tests
julia --project -e 'using Pkg; Pkg.test()'

# Run a specific test file
julia --project -e 'include("test/1_basic_tests.jl")'

# Development REPL (load package in dev mode)
julia --project
# Then: using UnivariateFunctions

# Build documentation locally
julia --project=docs docs/make.jl
```

## Architecture

### Type Hierarchy

The package uses a 4-level type hierarchy for `UnivariateFunction` (abstract base type), with strict seniority ordering:

1. **`Undefined_Function`** - Behaves like `missing`; any operation with it returns undefined
2. **`PE_Function`** - Basic power-exponential form: `a*exp(b*(x-base))*(x-base)^d`
3. **`Sum_Of_Functions`** - Vector of PE_Functions (enables polynomial approximation)
4. **`Piecewise_Function`** - Different univariate function for each domain segment

**Critical constraint**: Types are strictly nested. A Sum_Of_Functions never contains another Sum_Of_Functions, and a Piecewise_Function never contains another Piecewise_Function.

### Source File Organization

| File | Purpose |
|------|---------|
| `0_structs_and_generic_reversals.jl` | Core type definitions and operator overloads |
| `1_undefined_function.jl` | Undefined type implementation |
| `2_pe_functions.jl` | evaluate(), derivative(), indefinite_integral() for PE_Function |
| `3_sum_of_functions.jl` | Sum_Of_Functions implementation |
| `4_piecewise_functions.jl` | Piecewise_Function implementation |
| `5_calculus.jl` | Definite integration functions |
| `6_splines_and_interpolation.jl` | Spline and interpolation methods |
| `7_regressions_and_approximation.jl` | OLS and Chebyshev approximation |
| `8_isotonic_regressions.jl` | Isotonic/monotonic regression (Pool Adjacent Violators) |
| `9_supersmoother.jl` | Friedman's SuperSmoother algorithm |
| `10_unimodals.jl` | Unimodal regression and cross-validation |
| `date_conversions.jl` | Date/DateTime to float conversion (base: 2000-01-01) |
| `z_plotting.jl` | VegaLite visualization |
| `z_serialisation.jl` | DataFrame conversion |

### Key Design Constraints

1. **No division or negative powers** - Banned to maintain analytical differentiability/integrability
2. **Base date overflow** - When multiplying PE_Functions with very different bases (e.g., 2010 vs 50), numerical overflow can occur. Workaround: change bases on paper before implementation
3. **Dates as year fractions** - All dates convert to year fractions from global base 2000-01-01; no date info stored in PE_Function
4. **No optimization** - No built-in root/optima finding functions

### Main Operations

All calculus operations return `UnivariateFunction`:
- `evaluate(f, x)` - Evaluate at point
- `derivative(f)` - Returns UnivariateFunction
- `indefinite_integral(f)` - Returns UnivariateFunction
- `evaluate_integral(f, a, b)` - Definite integral

Algebraic operations: `+`, `-`, `*` (scalar and function), `^` (non-negative integer powers)

### Regression Functions

- `supersmoother(x, y)` - Adaptive local regression
- `monotonic_regression(x, y; increasing=true)` - Piecewise linear monotonic
- `isotonic_regression(x, y)` - Step function via Pool Adjacent Violators
- `unimodal_regression(x, y)` - Peak/trough constrained
- `cv_shape_regression(x, y)` - Cross-validated shape selection

All accept both arrays and DataFrames with column symbols.

## Test Structure

Tests are in `test/` with files numbered 1-10 matching feature areas. Run individual tests by including the specific file.
