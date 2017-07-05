export
  simpleevolve





"""
    simpleevolve(globaloperator, initialstate, timepoint)
    simpleevolve(globaloperator, initialstate, timepoints)

Simulates the GKSL master equation accordin to the equation

``|result⟩⟩ = exp(timepoint*globaloperator)|initialstate⟩⟩``

where ``|⋅⟩⟩`` denotes vectorization. The function return unvectorized `result`.
List of point of time (`timepoints`) can be given. Points of time needs to be
nonnegative (you cannot go back in time). If `globaloperator` is of type `Matrix`,
the exponentation is done by `expm` function. If `globaloperator` is of type
`SparseMatrixCSC`, `expmv` from `Expokit.jl` is used.

# Examples

```jldoctest
julia> H, L = [0 1; 1 0], [[0 1; 0 0],[0 0; 1 0]]
(
[0 1; 1 0],

Array{Int64,2}[
[0 1; 0 0],

[0 0; 1 0]])

julia> simpleevolve(globaloperator(H, L), proj(1,2), 4.)
2×2 Array{Complex{Float64},2}:
 0.499815-0.0im              0.0-0.00127256im
      0.0+0.00127256im  0.500185-0.0im

julia> simpleevolve(globaloperator(H, L), proj(1,2), [1.,2.,3.,4.])
4-element Array{Array{Complex{Float64},2},1}:
 Complex{Float64}[0.433203-0.0im 0.0-0.107605im; 0.0+0.107605im 0.566797-0.0im]
 Complex{Float64}[0.485766-0.0im 0.0+0.0171718im; 0.0-0.0171718im 0.514234-0.0im]
 Complex{Float64}[0.505597-0.0im 0.0+0.00261701im; 0.0-0.00261701im 0.494403-0.0im]
 Complex{Float64}[0.499815-0.0im 0.0-0.00127256im; 0.0+0.00127256im 0.500185-0.0im]

```
"""
function simpleevolve{T<:Number,S<:Number}(globaloperator::Matrix{T},
                                           initialstate::SparseDenseMatrix{S},
                                           timepoint::Real)
  @argument timepoint>=0 "Time needs to be nonnegative"
  unreshuffle(expm(timepoint*globaloperator)*reshuffle(initialstate))
end

function simpleevolve{T<:Number,S<:Number}(globaloperator::SparseMatrixCSC{T},
                                           initialstate::Matrix{S},
                                           timepoint::Real)
  @argument timepoint>=0 "Time needs to be nonnegative"
  unreshuffle(expmv(timepoint, globaloperator, reshuffle(initialstate)))
end

function simpleevolve{T<:Number,S<:Number,R<:Real}(globaloperator::SparseDenseMatrix{T},
                               initialstate::SparseDenseMatrix{S},
                               timepoints::Vector{R})
  @argument all(timepoints.>=0) "All timepoints needs to be nonnegative"

  [simpleevolve(globaloperator,initialstate,t) for t=timepoints]
end
