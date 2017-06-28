export
  classicallindbladoperators,
  globaloperator

"""
classicallindbladoperators(A[, ϵ])
metho

# Examples

```jldoctest
julia>
```
"""
function classicallindbladoperators(A::SparseDenseMatrix, ϵ::Float64=eps(1.))
    L = SparseMatrixCSC{Complex128,Int64}[]
    for i=1:size(A,1),j=1:size(A,1)
        if abs(A[i,j]) >= ϵ
            push!(L, A[i,j]*ketbra(i,j,size(A,1)))
        end
    end
    L
end



"""

    globaloperator(H, L[, locH][, w])

# Arguments
- `H::SparseMatrixCSC{Complex128,Int64}`, `H::Matrix{Complex128}`: Hamiltonian operator,
-
-
- `w::Float64=-1.`: scaling parameter, should be in [0,1].
# Return

# Examples

```jldoctest
julia>
```
"""
function globaloperator{T<:SparseDenseMatrix}(H::SparseDenseMatrix,  L::Vector{T},
  locH::SparseDenseMatrix=spzeros(size(H,1),size(H,1)); w::Real=-1.)
    if 0 <= w <= 1
      α = w
      β = 1.-w
    else
      α = β = 1.
    end
    F = spzeros(size(H,1)^2,size(H,1)^2)
    id = speye(size(H,1),size(H,1))
    for i = 1:length(L)
        F += kron(L[i],conj(L[i]))-0.5*kron(L[i]'*L[i],id)-0.5*kron(id,transpose(L[i])*conj(L[i]))
    end
    F += im*(kron(id,conj(locH))-kron(locH,id))
    F = α*F + β*im*(kron(id,conj(H))-kron(H,id))
    F
end
