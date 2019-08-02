export
  local_lind,
  evolve_generator

"""
    local_lind(A[; epsilon])

Split the elements of the matrix `A` into a collection of sparse matrices with
exactly one non-zero element. Martices are created if the absolute value of the
nonzero element is there are not smaller than `epsilon`, where `epsilon` should
be nonnegative. The `epsilon` defaults to `eps()` if not specified.

# Examples
```jldoctest; setup = :(using QSWalk)
julia> A = [1. 2.; 3. 4.]
2×2 Array{Float64,2}:
 1.0  2.0
 3.0  4.0

julia> local_lind(A)
4-element Array{SparseArrays.SparseMatrixCSC{Float64,Ti} where Ti<:Integer,1}:

  [1, 1]  =  1.0

  [1, 2]  =  2.0

  [2, 1]  =  3.0

  [2, 2]  =  4.0

julia> local_lind(A, epsilon = 1.5)
3-element Array{SparseArrays.SparseMatrixCSC{Float64,Ti} where Ti<:Integer,1}:

  [1, 2]  =  2.0

  [2, 1]  =  3.0

  [2, 2]  =  4.0

```
"""
function local_lind(A::Matrix{T}; epsilon::Real = eps()) where T<:Number
  @argumentcheck epsilon >= 0 "Epsilon should be nonegative"

  L = SparseMatrixCSC{T}[]
  for i = 1:size(A, 1), j = 1:size(A, 2)
    if abs(A[i, j]) >=  epsilon
        push!(L, A[i, j]*ketbra(i, j, size(A, 1)))
    end
  end
  L
end

function local_lind(A::SparseMatrixCSC{T}; epsilon::Real = eps()) where T<:Number
  @argumentcheck epsilon>= 0 "Epsilon should be nonegative"

  L = SparseMatrixCSC{T}[]
  for i = 1:size(A, 1), j = A[i, :].nzind
    if abs(A[i, j]) >=  epsilon
        push!(L, A[i, j]*ketbra(i, j, size(A, 1)))
    end
  end
  L
end

"""

    evolve_generator_create(H, L, localH, α, β)

Internal function for creating the generator for the evolution superoperator.
Given Hamiltonian `H`, collection of Lindblad operator `L`, local Hamiltonian
`localH` and scaling parameters `α` and `β` the function computes

``-i α (H ⊗ 1 - 1 ⊗ H) + β (-i(localH ⊗ 1-1 ⊗ localH)+∑ (L ⊗ L̄ - 1/2(L^†L ⊗ 1 + 1 ⊗ L^T L̄ )))``

where α and β should sum to one, or both be equal to one. In the later case
there are ignored.

*Note:* The function does not check whether `H` or `localH` are hermitian.

# Arguments
- `H`: Hamiltonian, must be hermitian,
- `L`: collection of Lindblad operators, each must be of the same size as `H`,
- `locH`: local Hamiltonian, suggested for nonmoralized QS walk, must be
hermitian and of the size of `H`,
- `α`: scaling parameter corresponding to Hamiltonian, should be in [0, 1],
- `β`: scaling parameter corresponding to Lindbladian part, should be in [0, 1].

# Return
The function return the generator, which can be used in `evolve` function.
"""
function evolve_generator_create(H::AbstractMatrix{<:Number},
                                 L::AbstractVector{<:AbstractMatrix{<:Number}},
                                 localH::AbstractMatrix{<:Number},
                                 α::Real,
                                 β::Real)
  @argumentcheck size(H) !=  (0, 0) "Matrix H must not be sizeless"
  @argumentcheck size(H, 1) == size(H, 2) "Matrix H must be square"
  @assert all([size(lindbladian) == size(H) for lindbladian in L]) "Lindblad operators must be of the same size as Hamiltonian"
  @assert size(H) == size(localH) "Matrix localH must be of the same size as H"
  @argumentcheck 0 <=  α <=  1 && 0 <=  β <=  1 "Value of ω must be nonngeative and smaller than one"

  F = spzeros(ComplexF64, (size(H).^2)...)
  id = Matrix{ComplexF64}(I, size(H)...)
  for i = 1:length(L)
      F +=  kron(L[i], conj(L[i]))-0.5*kron(L[i]'*L[i], id)-0.5*kron(id, transpose(L[i])*conj(L[i]))
  end
  F +=  im*(kron(id, conj(localH))-kron(localH, id))
  F = α*im*(kron(id, conj(H))-kron(H, id)) + β*F
  F
end

"""
    evolve_generator(H, L[, localH][, ω])

Create the generator for the evolution superoperator. Given Hamiltonian `H`,
collection of Lindblad operator `L`, local Hamiltonian `localH` and scaling
parameter `ω`, the generator is obtained as a sum

``-i(1-ω) (H ⊗ 1 - 1 ⊗ H) + ω (-i(localH ⊗ 1 - 1 ⊗ localH) + ∑(L ⊗ L̄ - 1/2(L^†L ⊗ 1 + 1 ⊗ L^T L̄ )))``

The last two arguments are optional.

If `localH` is not given, it defaults to sparse zero matrix of the size of `H`.

If `ω` is not given, both parts are taken with the same intensity and the global
operator takes the form

``-i(H ⊗ 1 - 1 ⊗ H) + (-i(localH ⊗ 1 - 1 ⊗ localH) + ∑(L ⊗ L̄ - 1/2(L^†L ⊗ 1 + 1 ⊗ L^T L̄ )))``

# Arguments
- `H`: Hamiltonian, must be hermitian,
- `L`: collection of Lindblad operators, each must be of the same size as `H`,
- `localH`: local Hamiltonian, suggested for nonmoralized QS walk, must be hermitian and of the size of `H`,
- `ω`: scaling parameter, should be in [0, 1].

# Return
The generator matrix, which can be used in `evolve` function.


# Examples

```jldoctest; setup = :(using QSWalk)
julia> H, L, localH = [0 1+im; 1-im 0], [0. 1; 0 0], [1. 0.; 0. 1.]
(Complex{Int64}[0+0im 1+1im; 1-1im 0+0im], [0.0 1.0; 0.0 0.0], [1.0 0.0; 0.0 1.0])

julia> evolve_generator(H, [L], localH, 1/2)
4×4 Array{Complex{Float64},2}:
  0.0+0.0im    0.5+0.5im    0.5-0.5im   0.5+0.0im
 -0.5+0.5im  -0.25+0.0im    0.0+0.0im   0.5-0.5im
 -0.5-0.5im    0.0+0.0im  -0.25+0.0im   0.5+0.5im
  0.0+0.0im   -0.5-0.5im   -0.5+0.5im  -0.5+0.0im

```
"""
function evolve_generator(H::AbstractMatrix{<:Number},
                          L::AbstractVector{<:AbstractMatrix{<:Number}},
                          localH::AbstractMatrix{<:Number},
                          ω::Real)
  evolve_generator_create(H, L, localH, 1-ω, ω)
end,

function evolve_generator(H::AbstractMatrix{<:Number},
                          L::AbstractVector{<:AbstractMatrix{<:Number}},
                          localH::AbstractMatrix{<:Number} = spzeros(eltype(H), size(H)...))
  evolve_generator_create(H, L, localH, 1., 1.)
end,

function evolve_generator(H::AbstractMatrix{T},
                          L::AbstractVector{<:AbstractMatrix{<:Number}},
                          ω::Real) where T<:Number
  evolve_generator_create(H, L, spzeros(T, size(H)...), 1-ω , ω)
end
