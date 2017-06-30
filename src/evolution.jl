export
  distributionsummation,
  simpleevolve


"""

    distributionsummation(probability, partition)

# Arguments
-
-

# Return

# Examples

```jldoctest
julia>
```
"""
function distributionsummation{T<:Real}(probability::Vector{T},partition::Vector{Vector{Int}})
    [sum(probability[block]) for block=partition]
end

function distributionsummation{T<:Number}(state::SparseDenseMatrix{T},partition::Vector{Vector{Int}})
    distributionsummation(Vector{Real}(diag(state)),partition)
end


"""
    simpleevolve()

# Examples

```jldoctest
julia>
```
"""
function simpleevolve{T<:Number}(globaloperator::Matrix{T},
  initialstate::SparseDenseMatrix,  timepoint::Real)
  unres(expm(timepoint*globaloperator)*res(initialstate))
end

function simpleevolve{T<:Number,S<:Number}(globaloperator::SparseMatrixCSC{T},
  initialstate::Matrix{S}, timepoint::Real)
  unres(expmv(timepoint, globaloperator, res(initialstate)))
end

function simpleevolve{S<:Real}(globaloperator::SparseDenseMatrix,
  initialstate::SparseDenseMatrix,  timepoints::Vector{S})
  [simpleevolve(globaloperator,initialstate,t) for t=timepoints]
end
