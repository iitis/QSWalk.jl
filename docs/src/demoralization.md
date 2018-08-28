```@meta
DocTestSetup = quote
   using QSWalk
end
```

## Nonmoralizing model

Global interaction quantum stochastic walk suffers for creating additional connections. For
removing such effect, nonmoralizing quantum stochastic walk was introduced, see [here](http://www.rintonpress.com/journals/qiconline.html#v17n1112).
Such model is constructed in several steps. First, the dimensionality is increased, hence to each
vertex multidimensional subspaces is attached. Then, Hamiltonian and Lindblad
operator is increased, furthermore additional Hamiltonian called "local Hamiltonian".
is introduced.

Below we present additional functionalities typical for nonmoralizing quantum
stochastic walk. By default the the operator is generalized as in the original [paper](http://www.rintonpress.com/journals/qiconline.html#v17n1112).

```@index
Order = [:type, :function]
Modules = [QSWalk]
Pages   = ["demoralization.md"]
```

## Full docs

```@docs
Vertex
VertexSet
nm_measurement
nm_loc_ham
nm_init
default_nm_loc_ham
make_vertex_set
vlist
subspace
fourier_matrix
vertexsetsize
```
