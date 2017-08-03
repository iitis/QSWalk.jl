export
  SparseDenseMatrix,
  SparseDenseVector,
  Vertex,
  VertexSet

import Base: ==, hash, getindex, length

SparseDenseMatrix = Union{SparseMatrixCSC{T} where T<:Number, Matrix{S} where S<:Number}
SparseDenseVector = Union{SparseVector{T} where T<:Number, Vector{S} where S<:Number}

"""
    type Vertex

The type consists of list of `Int` which describes the labels of vectors from
cannonical basis which correspond to the `Vertex`, see [1].

[1] Domino, K., Glos, A., & Ostaszewski, M. (2017). Spontaneous moralization
problem in quantum stochastic walk. arXiv preprint arXiv:1701.04624.
"""
immutable Vertex
  linspace::Vector{Int}

  Vertex(v::Vector{Int}) = all(v.>0) ? new(v) : throw(ArgumentError("vector should consist of positive elments"))
end

(v::Vertex)() = v.linspace
==(v::Vertex, w::Vertex) = v() == w()
hash(v::Vertex) = hash(v())

==(v::Tuple{Vertex,Vertex}, w::Tuple{Vertex,Vertex}) = [v[1](),v[2]()] == [w[1](),w[2]()]
hash(v::Tuple{Vertex,Vertex}) = hash([v[1](),v[2]()])
getindex(v::Vertex, i::Int) = v.linspace[i]
length(v::Vertex) = length(v())


function checkVertexSet(partition::Vector{Vector{Int}})
  joined = Int[]
  for lin=partition
    append!(joined,lin)
  end
  length(joined) == length(Set(joined))
end

checkVertexSet(vset::Vector{Vertex}) = checkVertexSet([vertex.linspace for vertex=vset])


"""
    type VertexSet

The type consists of list of `Vertex` which describes the partition of
linear subspace, see [1]

[1] Domino, K., Glos, A., & Ostaszewski, M. (2017). Spontaneous moralization
problem in quantum stochastic walk. arXiv preprint arXiv:1701.04624.
"""
type VertexSet
  vertices::Vector{Vertex}

  VertexSet{T<:Vertex}(vset::Vector{T}) = checkVertexSet(vset) ? new(vset) : throw(ArgumentError("Vertices should correspond to orthogonal linear spaces"))
  VertexSet(vset::Vector{Vector{Int}}) = checkVertexSet(vset) ? new([Vertex(v) for v=vset]) : throw(ArgumentError("Vertices should correspond to orthogonal linear spaces"))
end


(v::VertexSet)() = v.vertices

==(v::VertexSet,w::VertexSet) = v.vertices == w.vertices

getindex(vset::VertexSet, i::Int) = vset()[i]
getindex(vset::VertexSet, veci::Vector{Int}) = vset()[veci]
length(vset::VertexSet) = length(vset())
"""

    partitionsize(vertexset)

Return the dimenions of linearspace corresponding to given `vertexset'.

# Examples

```jldoctest
julia> QSWalk.partitionsize(VertexSet([[1,2,3],[4,5]]))
5
```
"""
function vertexsetsize(vertexset::VertexSet)
  sum([length(vertex()) for vertex=vertexset()])
end


macro argument(ex, msgs...)
    msg = isempty(msgs) ? ex : msgs[1]
    if isa(msg, AbstractString)
        msg = msg # pass-through
    elseif !isempty(msgs) && (isa(msg, Expr) || isa(msg, Symbol))
        # message is an expression needing evaluating
        msg = :(Main.Base.string($(esc(msg))))
    elseif isdefined(Main, :Base) && isdefined(Main.Base, :string) && applicable(Main.Base.string, msg)
        msg = Main.Base.string(msg)
    else
        # string() might not be defined during bootstrap
        msg = :(Main.Base.string($(Expr(:quote,msg))))
    end
    return :($(esc(ex)) ? $(nothing) : throw(Main.Base.ArgumentError($msg)))
end
