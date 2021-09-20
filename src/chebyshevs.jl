const number_of_chebyshevs_to_compile_into_binaries = 20

function next_chebyshev(previous::Sum_Of_Functions, two_previous::Sum_Of_Functions)
    return PE_Function(2.0,0.0,0.0,1) * previous - two_previous
end

first_kind_chebyshevs = Array{Sum_Of_Functions}(undef, number_of_chebyshevs_to_compile_into_binaries)
first_kind_chebyshevs[1] = Sum_Of_Functions([PE_Function(1.0,0.0,0.0,0)])
first_kind_chebyshevs[2] = Sum_Of_Functions([PE_Function(1.0,0.0,0.0,1)])
for i in 3:number_of_chebyshevs_to_compile_into_binaries
    first_kind_chebyshevs[i] = next_chebyshev(first_kind_chebyshevs[i-1], first_kind_chebyshevs[i-2])
end

second_kind_chebyshevs = Array{Sum_Of_Functions}(undef, number_of_chebyshevs_to_compile_into_binaries)
second_kind_chebyshevs[1] = Sum_Of_Functions([PE_Function(1.0,0.0,0.0,0)])
second_kind_chebyshevs[2] = Sum_Of_Functions([PE_Function(2.0,0.0,0.0,1)])
for i in 3:number_of_chebyshevs_to_compile_into_binaries
    second_kind_chebyshevs[i] = next_chebyshev(second_kind_chebyshevs[i-1], second_kind_chebyshevs[i-2])
end


"""
    get_chevyshevs_up_to(N::Integer, first_kind::Bool = true)

Get the first N chebyshev polynomials returned as a vector of UnivariateFunctions.
The first 20 polynomials of each are precompiled into the binaries for speed. If
you need more than that they will be calculated at runtime.

These can be from either the first kind or second kind polynomial sequence.
### Inputs
* `N` - How many chebyshev polynomials do you want.
* `first_kind` - A Bool. If true you get first kind polynomials. If false you get second kind.
### Returns
* A `Vector` of `UnivariateFunction`s for each polynomial.
"""
function get_chevyshevs_up_to(N::Integer, first_kind::Bool = true)
    chebyshevs = Array{Sum_Of_Functions}(undef, N)
    if N >= number_of_chebyshevs_to_compile_into_binaries
        if first_kind
            chebyshevs[1:number_of_chebyshevs_to_compile_into_binaries] = first_kind_chebyshevs
        else
            chebyshevs[1:number_of_chebyshevs_to_compile_into_binaries] = second_kind_chebyshevs
        end
        for i in (number_of_chebyshevs_to_compile_into_binaries+1):N
            chebyshevs[i] = next_chebyshev(chebyshevs[i-1], chebyshevs[i-2])
        end
    else
        if first_kind
            chebyshevs[1:N] = first_kind_chebyshevs[1:N]
        else
            chebyshevs[1:N] = second_kind_chebyshevs[1:N]
        end
    end
    return chebyshevs
end
