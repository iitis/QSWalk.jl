export
  SparseDenseMatrix,
  SparseDenseVector,
  Vertex,
  VertexSet,
  vertexsetsize,
  subspace,
  vertices

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
struct Vertex
  subspace::Vector{Int}

  Vertex(v::Vector{Int}) = all(v.>0) ? new(v) : throw(ArgumentError("vector should consist of positive elments"))
end

subspace(v::Vertex) = v.subspace

==(v::Vertex, w::Vertex) = subspace(v) ==  subspace(w)
hash(v::Vertex) = hash(subspace(v))

==(v::Tuple{Vertex, Vertex}, w::Tuple{Vertex, Vertex}) = [subspace(v[1]), subspace(v[2])] ==  [subspace(w[1]), subspace(w[2])]
hash(v::Tuple{Vertex, Vertex}) = hash([subspace(v[1]), subspace(v[2])])
getindex(v::Vertex, i::Int) = v.subspace[i]
length(v::Vertex) = length(subspace(v))


function checkvertexset(partition::Vector{Vector{Int}})
  joined = Int[]
  for lin = partition
    append!(joined, lin)
  end
  length(joined) ==  length(Set(joined))
end

checkvertexset(vset::Vector{Vertex}) = checkvertexset([subspace(v) for v = vset])


"""
    type VertexSet

Type consisting of a list of `Vertex` objects. It describes the partition of the
linear subspace. Should be constructed by `make_vertex_set` or by `nonmoralizing_lindbladian`
functions. In order to get list of the vertices of the object `vertexset`, one
should use the function `vertexset()`, of for concrete `Vertex` an getindex function
`vertexset[i]`.

[1] K. Domino, A. Glos, M. Ostaszewski, Superdiffusive quantum stochastic walk
definable on arbitrary directed graph, Quantum Information & Computation,
Vol.17 No.11&12, pp. 0973-0986, arXiv:1701.04624.
"""
struct VertexSet
  vertices::Vector{Vertex}

  VertexSet(vset::Vector{Vertex}) = checkvertexset(vset) ? new(vset) : throw(ArgumentError("Vertices should correspond to orthogonal linear spaces"))
end

VertexSet(vset::Vector{Vector{Int}}) = VertexSet([Vertex(v) for v = vset])

vlist(vset::VertexSet) = vset.vertices

==(v::VertexSet, w::VertexSet) = v.vertices ==  w.vertices

getindex(vset::VertexSet, i::Int) = vlist(vset)[i]
getindex(vset::VertexSet, veci::Vector{Int}) = vlist(vset)[veci]
length(vset::VertexSet) = length(vlist(vset))

"""
    vertexsetsize(vertexset)

Return the dimenion of the linearspace corresponding to given `vertexset'.

# Examples

```jldoctest
julia> vertexsetsize(VertexSet([[1, 2, 3], [4, 5]]))
5
```
"""
function vertexsetsize(vset::VertexSet)
  sum(length.(vlist(vset)))
end

macro argumentcheck(ex, msgs...)
    msg = isempty(msgs) ? ex : msgs[1]
    if isa(msg, AbstractString)
        msg = msg # pass-through
    elseif !isempty(msgs) && (isa(msg, Expr) || isa(msg, Symbol))
        #message is an expression needing evaluating
        msg = :(Main.Base.string($(esc(msg))))
    #elseif isdefined(Main, :Base) && isdefined(Main.Base, :string) && applicable(Main.Base.string, msg)
        #msg = Main.Base.string(msg)
    #else
        # string() might not be defined during bootstrap
        #msg = :(Main.Base.string($(Expr(:quote, msg))))
    end
    return :($(esc(ex)) ? $(nothing) : throw(Main.Base.ArgumentError($msg)))
end
