export
  initialstate,
  distributionsummation

"""

    distributionsummation(probability, vertexset)
    distributionsummation(state, vertexset)

Returns joint probabilty of `probability`, which is real-valued probability vector
according to `vertexset`.

Returns joint probabilty of cannonical measurement of density matrix `state`,
according to `vertexset`.

# Examples

```jldoctest
julia> probability = [0.05,0.1,0.25,0.3,0.01,0.20,0.04,0.05]
8-element Array{Float64,1}:
^[[A^[[A 0.05
 0.1
 0.25
 0.3
 0.01
 0.2
 0.04
 0.05

julia> distributionsummation(probability, VertexSet([[1,4],[2,3,5],[6],[7,8]]))
4-element Array{Float64,1}:
 0.35
 0.36
 0.2
 0.09

 julia> state = [1/6 0 1/6; 0 2/3 0; 1/6 0 1/6]
 3×3 Array{Float64,2}:
  0.166667  0.0       0.166667
  0.0       0.666667  0.0
  0.166667  0.0       0.166667

 julia> distributionsummation(state, VertexSet([[1,3],[2]]))
 2-element Array{Float64,1}:
  0.333333
  0.666667
```
"""
function distributionsummation{T<:Real}(probability::Vector{T},
                                        vertexset::VertexSet)
    [sum(probability[vertex()]) for vertex=vertexset()]
end

function distributionsummation{T<:Number}(state::SparseDenseMatrix{T},
                                          vertexset::VertexSet)
    distributionsummation(real.(diag(state)),vertexset)
end

"""

    initialstate(initialvertices, vertexset)
    initialstate(initialstates, vertexset)

Function creating initial state in the case of demoralized evolution. It returns
a block diagonal matrix, where each block correspond to vertex from `vertexset`.
If first argument is of type `Vector{Vertex}`, then default block matrix `eye` (see
example for more details).

If first argument is of type `Dict{Vertex,SparseDenseMatrix}`,
then for each given vertex a block from dictionary is used, otherwise zero matrix
is chosen. Each matrix from dictionary should be nonnegative and sum of all traces
should equal one. The keys of `initialstates` shuold be a subset of `vertexset()`.
Note that matrix from `initialstates` corresponding to vertex `v` should be of
size `length(v)`×`length(v)`.

The function returns sparse matrix with `Complex128` field type.

# Examples

```jldoctest
julia> vset = VertexSet([[1],[2,3,4],[5,6,7],[8,9]])
QSWalk.VertexSet(QSWalk.Vertex[QSWalk.Vertex([1]),QSWalk.Vertex([2,3,4]),QSWalk.Vertex([5,6,7]),QSWalk.Vertex([8,9])])

julia> initialstate(vset[[1,3,4]],vset)
9×9 sparse matrix with 6 Complex{Float64} nonzero entries:
	[1, 1]  =  0.333333+0.0im
	[5, 5]  =  0.111111+0.0im
	[6, 6]  =  0.111111+0.0im
	[7, 7]  =  0.111111+0.0im
	[8, 8]  =  0.166667+0.0im
	[9, 9]  =  0.166667+0.0im

julia> A1, A2, A3 = ones(Complex128,1,1)/4, [ 1/5+0im 0 1/5; 0 1/10 0 ; 1/5 0 1/5 ], [0.125 -0.125+0im; -0.125 0.125]
(
Complex{Float64}[0.25+0.0im],

Complex{Float64}[0.2+0.0im 0.0+0.0im 0.2+0.0im; 0.0+0.0im 0.1+0.0im 0.0+0.0im; 0.2+0.0im 0.0+0.0im 0.2+0.0im],

Complex{Float64}[0.125+0.0im -0.125+0.0im; -0.125+0.0im 0.125+0.0im])

julia> initialstate(Dict(vset[1]=>A1, vset[3]=>A2, vset[4]=>A3), vset)
9×9 sparse matrix with 10 Complex{Float64} nonzero entries:
	[1, 1]  =  0.25+0.0im
	[5, 5]  =  0.2+0.0im
	[7, 5]  =  0.2+0.0im
	[6, 6]  =  0.1+0.0im
	[5, 7]  =  0.2+0.0im
	[7, 7]  =  0.2+0.0im
	[8, 8]  =  0.125+0.0im
	[9, 8]  =  -0.125+0.0im
	[8, 9]  =  -0.125+0.0im
	[9, 9]  =  0.125+0.0im
```
"""
function initialstate(initialvertices::Vector{Vertex}, vertexset::VertexSet)
  L = spzeros(Complex128, vertexsetsize(vertexset),vertexsetsize(vertexset))
  for vertex=initialvertices
    normalization = length(vertex)*length(initialvertices)
    L[vertex(),vertex()] = speye(Complex128,length(vertex))/normalization
  end
  L
end

function initialstate{T<:SparseDenseMatrix}(initialstates::Dict{Vertex,T},
                                            vertexset::VertexSet)
  L = spzeros(Complex128, vertexsetsize(vertexset), vertexsetsize(vertexset))
  for vertex=keys(initialstates)
    L[vertex(),vertex()] = initialstates[vertex]
  end
  L
end
