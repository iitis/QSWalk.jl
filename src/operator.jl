export
  classicallindbladoperators,
  globaloperator


"""
classicallindbladoperators(A[, ϵ=1e-7])


# Examples

```jldoctest
julia>
```
"""
function classicallindbladoperators(A::SparseDenseMatrix; ϵ::Float64=1e-7)
    L = SparseMatrixCSC{Complex128,Int64}[]
    for i=1:size(A,1),j=1:size(A,1)
        if abs(A[i,j]) >= ϵ
            push!(L, A[i,j]*ketbra(i,j,size(A,1)))
        end
    end
    L
end



"""

    globaloperator(H, L[, locH=spzeros(Complex128,size(H,1),size(H,1))][, w=-1])

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
function globaloperator{T::SparseDenseMatrix}(::Type{T},
  H::SparseDenseMatrix{U},  L::Vector{SparseDenseMatrix},
  locH::SparseDenseMatrix=spzeros(typeof(H[1]),size(H,1),size(H,1)); w::Real=-1.)
    #=if 0 <= w <= 1
      α = w
      β = 1.-w
    else
      α = β = 1.
    end
    F = convert(T,spzeros(size(H,1)^2,size(H,1)^2))
    id = speye(H)
    for i = 1:size(L,1)
        F += kron(L[i],conj(L[i]))-0.5*kron(L[i]'*L[i],id)-0.5*kron(id,transpose(L[i])*conj(L[i]))
    end
    F += im*(kron(id,conj(locH))-kron(locH,id))
    F = α*F + β*im*(kron(id,conj(H))-kron(H,id))
    convert(T,F)=#
end
