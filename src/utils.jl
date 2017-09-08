export
  SparseDenseMatrix,
  SparseDenseVector,
  Vertex,
  VertexSet,
  vertexsetsize

import Base: ==, hash, getindex, length

SparseDenseMatrix = Union{SparseMatrixCSC{T} where T<:Number, Matrix{S} where S<:Number}
SparseDenseVector = Union{SparseVector{T} where T<:Number, Vector{S} where S<:Number}

"""
    type Vertex

Type consisting of list of `Int`, describing the labels of vectors from the
cannonical basis which correspond to the `Vertex`. See [1] for the more
information and usage exmaples.

To get the vector label one can use `Vertex()` function, or `Vertex[i]` for
unique label.

[1] K. Domino, A. Glos, M. Ostaszewski, Superdiffusive quantum stochastic walk
definable on arbitrary directed graph, Quantum Information & Computation,
Vol.17 No.11&12, pp. 0973-0986, arXiv:1701.04624.
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

Type consisting of a list of `Vertex` objects. It describes the partition of the
linear subspace. Should be constructed by `make_vertex_set` or by `demoralized_lindbladian`
functions. In order to get list of the vertices of the object `vertexset`, one
should use the function `vertexset()`, of for concrete `Vertex` an getindex function
`vertexset[i]`.

[1] K. Domino, A. Glos, M. Ostaszewski, Superdiffusive quantum stochastic walk
definable on arbitrary directed graph, Quantum Information & Computation,
Vol.17 No.11&12, pp. 0973-0986, arXiv:1701.04624.
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
    vertexsetsize(vertexset)

Return the dimenion of the linearspace corresponding to given `vertexset'.

# Examples

```jldoctest
julia> vertexsetsize(VertexSet([[1,2,3],[4,5]]))
5
```
"""
function vertexsetsize(vertexset::VertexSet)
  sum([length(vertex()) for vertex=vertexset()])
end


macro argument(ex, msgs...)
    msg = isempty(msgs) ? ex : msgs[1]
    return :($(esc(ex)) ? $(nothing) : throw(Main.Base.ArgumentError($msg)))
end
