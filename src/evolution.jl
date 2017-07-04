export
  distributionsummation,
  simpleevolve


"""

    distributionsummation(probability, vertexset)
    distributionsummation(state, vertexset)

Returns joint probabilty of `probability`, which is real-valued probability vector
according to `vertexset`.

Returns joint probabilty of cannonical measurement of density matrix `state`,
according to `vertexset`.

# Examples



```jldoctest
julia>
```
"""
function distributionsummation{T<:Real}(probability::Vector{T},
                                        vertexset::VertexSet)
    [sum(probability[vertex]) for vertex=vertexset()]
end

function distributionsummation{T<:Number}(state::SparseDenseMatrix{T},
                                          vertexset::VertexSet)
    distributionsummation(Vector{Real}(diag(state)),vertexset)
end


"""
    simpleevolve(globaloperator, initialstate, timepoint)
    simpleevolve(globaloperator, initialstate, timepoints)

Simulates the GKSL master equation accordin to the equation

``|result⟩⟩ = exp(timepoint*globaloperator)|initialstate⟩⟩``

where ``|⋅⟩⟩`` denotes vectorization. The function return unvectorize `result`.
List of point of time (`timepoints`) can be given. Points of time needs to be
nonnegative (you cannot go back in time). If `globaloperator` is of type `Matrix`,
the exponentation is done by `expm` function. If `globaloperator` is of type
`SparseMatrixCSC`, `expmv` from `Expokit.jl` is used.

# Examples

```jldoctest
julia>
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

function simpleevolve{T<:Number,S<:Number}(globaloperator::SparseDenseMatrix{T},
                               initialstate::SparseDenseMatrix{S},
                               timepoints::Vector{S})
  @argument all(timepoints.>=0) "All timepoints needs to be nonnegative"

  [simpleevolve(globaloperator,initialstate,t) for t=timepoints]
end
