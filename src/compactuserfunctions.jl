export
  init_nonmoralized,
  measurement_nonmoralized

"""

    measurement_nonmoralized(probability, vertexset)
    measurement_nonmoralized(state, vertexset)

Returns joint probability of `probability`, which is real-valued probability vector
according to `vertexset`.

Returns joint probability of cannonical measurement of density matrix `state`,
according to `vertexset`.

*note* It is up to user to provide proper probability vector/density state.

# Examples

```jldoctest
julia> probability = [0.05,0.1,0.25,0.3,0.01,0.20,0.04,0.05]
8-element Array{Float64,1}:
 0.05
 0.1
 0.25
 0.3
 0.01
 0.2
 0.04
 0.05

julia> measurement_nonmoralized(probability, VertexSet([[1,4],[2,3,5],[6],[7,8]]))
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

 julia> measurement_nonmoralized(state, VertexSet([[1,3],[2]]))
 2-element Array{Float64,1}:
  0.333333
  0.666667
```
"""
function measurement_nonmoralized(probability::Vector{T} where T<:Number,
                                  vertexset::VertexSet)
  @assert vertexsetsize(vertexset) == length(probability) "vertexset size and probability vector length do not match"

  [sum(probability[vertex()]) for vertex=vertexset()]
end

function measurement_nonmoralized(state::SparseDenseMatrix,
                                  vertexset::VertexSet)
  @argument size(state, 1) == size(state, 2) "state should be square matrix"
  @assert vertexsetsize(vertexset) == size(state, 1) "vertexset size and state size do not match"

  measurement_nonmoralized(real.(diag(state)),vertexset)
end

"""

    init_nonmoralized(initialvertices, vertexset)
    init_nonmoralized(initial_states, vertexset)

Function creating initial state in the case of nonmoralizing evolution. It returns
a block diagonal matrix, where each block correspond to vertex from `vertexset`.
If first argument is of type `Vector{Vertex}`, then default block matrix `eye` (see
example for more details).

If first argument is of type `Dict{Vertex,SparseDenseMatrix}`,
then for each given vertex a block from dictionary is used, otherwise zero matrix
is chosen. Each matrix from dictionary should be nonnegative and sum of all traces
should equal one. The keys of `initialvertices` should be a subset of `vertexset()`.
Note that matrix from `initial_states` corresponding to vertex `v` should be of
size `length(v)`×`length(v)`.

The function returns sparse matrix with `Complex128` field type.

# Examples

```jldoctest
julia> vset = VertexSet([[1],[2,3,4],[5,6,7],[8,9]])
QSWalk.VertexSet(QSWalk.Vertex[QSWalk.Vertex([1]),QSWalk.Vertex([2,3,4]),QSWalk.Vertex([5,6,7]),QSWalk.Vertex([8,9])])

julia> init_nonmoralized(vset[[1,3,4]],vset)
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

julia> init_nonmoralized(Dict(vset[1]=>A1, vset[3]=>A2, vset[4]=>A3), vset)
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
function init_nonmoralized(initialvertices::Vector{Vertex},
                           vertexset::VertexSet)
  @assert all([v in vertexset() for v=initialvertices]) "initialvertices is not a subset of vertexset"

  L = spzeros(Complex128, vertexsetsize(vertexset),vertexsetsize(vertexset))
  for vertex=initialvertices
    normalization = length(vertex)*length(initialvertices)
    L[vertex(),vertex()] = speye(Complex128,length(vertex))/normalization
  end
  L
end

function init_nonmoralized(initial_states::Dict{Vertex,T} where T,
                           vertexset::VertexSet)
  @argument all([typeof(state)<:SparseDenseMatrix for state=values(initial_states)]) "All elements in `hamiltonians` must be SparseMatrixCSC or Matrix"
  @argument all([eltype(state)<:Number for state=values(initial_states)]) "All elements of elements in `hamiltonians` must be Number"
  @assert all([size(state, 1) == length(k) for (k,state)=initial_states]) "The size of initial state and the vertex do not match"
  @assert all([k in vertexset() for k=keys(initial_states)]) "keys of initial_states is not a subset of vertexset"

  L = spzeros(Complex128, vertexsetsize(vertexset), vertexsetsize(vertexset))
  for vertex=keys(initial_states)
    L[vertex(),vertex()] = initial_states[vertex]
  end
  L
end
