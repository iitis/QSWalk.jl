export
  classical_lindblad_operators,
  global_operator

"""
    classical_lindblad_operators(A[; epsilon])

The function splits the elements of the matrix `A` into collection of
single-nonzero element sparse matrices. Martices are added if the absolute value
of the nonzero element is not smaller than `epsilon`, hence `epsilon` should be
nonnegative. The `epsilon` defaults to `eps()` if not specified.

# Examples

```jldoctest
julia> A = [1. 2.; 3. 4.]
2×2 Array{Float64,2}:
 1.0  2.0
 3.0  4.0

julia> classical_lindblad_operators(A)
4-element Array{SparseMatrixCSC{Float64,Ti<:Integer},1}:

	[1, 1]  =  1.0

	[1, 2]  =  2.0

	[2, 1]  =  3.0

	[2, 2]  =  4.0

julia> classical_lindblad_operators(A, epsilon=1.5)
3-element Array{SparseMatrixCSC{Float64,Ti<:Integer},1}:

	[1, 2]  =  2.0

	[2, 1]  =  3.0

	[2, 2]  =  4.0
```
"""
function classical_lindblad_operators(A::Matrix{T} where T<:Number;
                                    epsilon::Real=eps())
  @argument epsilon>=0 "epsilon should be nonegative"
  L = SparseMatrixCSC{eltype(A)}[]
  for i=1:size(A,1), j=1:size(A,2)
    if abs(A[i,j]) >= epsilon
        push!(L, A[i,j]*ketbra(eltype(A),i,j,size(A,1)))
    end
  end
  L
end

function classical_lindblad_operators(A::SparseMatrixCSC{T} where T<:Number;
                                    epsilon::Real=eps())
  @argument epsilon>=0 "epsilon should be nonegative"
  L = SparseMatrixCSC{eltype(A)}[]
  for i=1:size(A,1), j=A[i,:].nzind
    if abs(A[i,j]) >= epsilon
        push!(L, A[i,j]*ketbra(eltype(A),i,j,size(A,1)))
    end
  end
  L
end

"""

    global_operator_util(H, L, localH, α, β)

Util function creating global operator for evolution. Given Hamiltonian `H`,
collection of Lindblad operator `L`, local Hamiltonian `localH` and scaling
parameters `α` and `β` the function computes

``-i α (H ⊗ 1 - 1 ⊗ H) + β (-i(localH ⊗ 1-1 ⊗ localH)+∑ (L ⊗ L̄ - 1/2(L\^†L ⊗ 1 + 1 ⊗ L\^T L̄ )))``

α and β should at the same time equal to one or sum to one. In the first case
their values are simply ignored, in the second they correspond to evolution used
in [1].

The function does not check whether `H` or `localH` are Hermitian operators.

[1] Domino, K., Glos, A., & Ostaszewski, M. (2017). Spontaneous moralization
problem in quantum stochastic walk. arXiv preprint arXiv:1701.04624.

# Arguments
- `H`: Hamiltonian, must be hermitian,
- `L`: collection of Lindblad operators, each must be of the same size as `H`,
- `locH`: local Hamiltonian, suggested for nonmoralized QS walk, must be hermitian and of the size of `H`,
- `α`: scaling parameter corresponding to Hamiltonian, should be in [0,1],
- `β`: scaling parameter corresponding to Lindbladian part, should be in [0,1].

# Return
The function return global operator used for evolution.

# Examples
```jldoctest
julia>
```
"""
function global_operator_util(H::SparseDenseMatrix,
                            L::Vector{T} where T<:AbstractArray,
                            localH::SparseDenseMatrix,
                            α::Real,
                            β::Real)
  @argument size(H) != (0,0) "H must not be sizeless"
  @assert all([size(lindbladian) == size(H) for lindbladian in L]) "Lindblad operators must be of the same size as Hamiltonian"
  @argument all([eltype(el)<:Number for el in L]) "Lindblad operators elements must be numbers"
  @argument all([typeof(el)<:SparseDenseMatrix for el in L]) "Lindblad operators must be SparseMatrixCSC or Matrix"
  @assert size(H) == size(localH) "localH must be of the same size as H"
  @argument 0 <= α <= 1 && 0 <= β <= 1 "ω must be nonngeative and smaller than one"

  F = spzeros(Complex128,(size(H).^2)...)
  id = eye(H)
  for i = 1:length(L)
      F += kron(L[i],conj(L[i]))-0.5*kron(L[i]'*L[i],id)-0.5*kron(id,transpose(L[i])*conj(L[i]))
  end
  F += im*(kron(id,conj(localH))-kron(localH,id))
  F = α*im*(kron(id,conj(H))-kron(H,id)) + β*F
  F
end

"""
    global_operator(H, L[, localH][, ω])

The function creates global operator for evolution. Given Hamiltonian `H`,
collection of Lindblad operator `L`, local Hamiltonian `localH` and scaling
parameter `ω`

``-i (1-ω) (H ⊗ 1 - 1 ⊗ H) + ω (-i(localH ⊗ 1-1 ⊗ localH)+∑ (L ⊗ L̄ - 1/2(L\^†L ⊗ 1 + 1 ⊗ L\^T L̄ )))``

The `localH` defaults to sparse zero matrix of the size `H` if not specified. If
`ω` is not given, the global operator takes the form

``-i (H ⊗ 1 - 1 ⊗ H) + (-i(localH ⊗ 1-1 ⊗ localH)+∑ (L ⊗ L̄ - 1/2(L\^†L ⊗ 1 + 1 ⊗ L\^T L̄ )))``

The formulas were given in [1].

[1] K. Domino, A. Glos, M. Ostaszewski, Superdiffusive quantum stochastic walk
definable on arbitrary directed graph, Quantum Information & Computation,
Vol.17 No.11&12, pp. 0973-0986, arXiv:1701.04624.

# Arguments
- `H`: Hamiltonian, must be hermitian,
- `L`: collection of Lindblad operators, each must be of the same size as `H`,
- `localH`: local Hamiltonian, suggested for nonmoralized QS walk, must be hermitian
and of the size of `H`,
- `ω`: scaling parameter, should be in [0,1].

# Examples

```jldoctest
julia> H, L, localH = [0 1+im; 1-im 0], [0. 1; 0 0], eye(2)
(
Complex{Int64}[0+0im 1+1im; 1-1im 0+0im],

[0.0 1.0; 0.0 0.0],

[1.0 0.0; 0.0 1.0])

julia> global_operator(H, [L], localH, 1/2)
4×4 Array{Complex{Float64},2}:
  0.0+0.0im    0.5+0.5im    0.5-0.5im   0.5+0.0im
 -0.5+0.5im  -0.25+0.0im    0.0+0.0im   0.5-0.5im
 -0.5-0.5im    0.0+0.0im  -0.25+0.0im   0.5+0.5im
  0.0+0.0im   -0.5-0.5im   -0.5+0.5im  -0.5+0.0im

```
"""
function global_operator(H::SparseDenseMatrix,
                        L::Vector{T} where T<:AbstractArray,
                        localH::SparseDenseMatrix,
                        ω::Real)
 global_operatorutil(H, L, localH, 1-ω, ω)
end

function global_operator(H::SparseDenseMatrix,
                        L::Vector{T} where T<:AbstractArray,
                        localH::SparseDenseMatrix=spzeros(eltype(H),size(H)...))
 global_operator_util(H, L, localH, 1., 1.)
end

function global_operator(H::SparseDenseMatrix,
                        L::Vector{T} where T<:AbstractArray,
                        ω::Real)
 global_operatorutil(H, L, spzeros(eltype(H),size(H)...), 1-ω , ω)
end
