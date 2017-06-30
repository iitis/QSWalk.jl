export
  evolve,
  initialstate,
  partitioning



function evolve{T<:Number,S<:Number}(H::Matrix{T}, L::Vector{Matrix{S}}, timepoint::Real;
  w=1.::Real)

end


"""

    initialstate()

# Examples

```jldoctest
julia>
```
"""
function initialstate(initialvertices::Vector{Int}, partition::Vector{Vector{Int}})
  L = spzeros(Complex128, partitionsize(partition),partitionsize(partition))
  for v=initialvertices
    normalization = length(partition[v])*length(initialvertices)
    L[partition[v],partition[v]] = speye(Complex128,length(partition[v]))/normalization
  end
  L
end

function initialstate{T<:SparseDenseMatrix}(initialstates::Vector{T})
  sizeL = sum([size(A,1) for A in initialstates])
  L = spzeros(Complex128, sizeL, sizeL)
  startindex = 1
  endindex = 0
  for state=initialstates
    endindex += size(state,1)
    L[startindex:endindex,startindex:endindex] = sparse(state)
    startindex += size(state,1)
  end
  L
end
