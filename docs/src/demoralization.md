## Nonmoralizing model

Global interaction quantum stochastic walk suffers from creating additional connections. This renders it unsuitable for constructing fast quantum walks on directed graphs. To counteract this effect, the nonmoralizing quantum stochastic walk was introduced (see [arXiv preprint](https://arxiv.org/abs/1701.04624) and [its published version](http://www.rintonpress.com/journals/qiconline.html#v17n1112)). Such a model is constructed in several steps. First, the dimensionality of the system is increased by attaching a multidimensional subspace to each vertex. Next, the Hamiltonian and the Lindblad operators are modified, and an additional Hamiltonian - so-called *local Hamiltonian* - is introduced.

*Please note that current definition of `nm_glob_ham` differs from the one presented in the paper.*

Below we present additional functionality useful for analyzing nonmoralizing quantum stochastic walk. By default, the operator is generalized as in the original [paper](http://www.rintonpress.com/journals/qiconline.html#v17n1112).



```@index
Order = [:type, :function]
Modules = [QSWalk]
Pages   = ["demoralization.md"]
```

## Full docs

```@docs
Vertex
VertexSet
nm_glob_ham
nm_lind
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
