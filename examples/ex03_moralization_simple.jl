
# coding: utf-8

# # Spontaneous moralization on simple graph

# ## Loading modules

using QSWalk
using LightGraphs 
using LinearAlgebra


# ## Moralizing evolution

# Here we provide the simple example of graph, for which the spontaneous moralization happens. Note there is no path from vertex 1 to 2.

digraph = DiGraph(3)
add_edge!(digraph, 1, 3)
add_edge!(digraph, 2, 3)


# Here we generate some basic operators. Note in the case of directed graphs we need to transpose adjacency matrix, as *QSWalk.jl* multiplies the state on the right side of evolution operator. Note we choose zero matrix as Hamiltonian of the system. 
# 
# As we deal with the graph of very small size, we choose full-matrix evolution algorithm. In order to do such, we need to provide in opmoral at least one full matrix.

adjacency = Matrix(transpose(adjacency_matrix(digraph)))
timepoint = 100.

opmoral = evolve_generator(zero(adjacency), [adjacency])
println(typeof(opmoral))


# As a result of the evolution we get an stationary state with non-zero probability of measuring vertex 2. Note the state is actually a stationary state.

rho = evolve(opmoral, proj(1, 3), timepoint)
println("Cannonical measurement on stationary state: $(real.(diag(rho)))")
println("Norm of opmoral times rho: $(norm(opmoral*res(rho)))")


# ## Non-moralizing evolution

# In this example we present model, which do not possess unintuitive property.

lnonmoral, vset = nm_lind(adjacency)
hlocal = nm_loc_ham(vset)
opnonmoral = evolve_generator(zero(lnonmoral), [lnonmoral], hlocal);


# *vset* and *hlocal* represent parametrization and additional operator for non-moralizing evolution. Note that vset is actually a partition of the new, increased linear space between vertices.

println("vset: $vset")
println("Subspace dimension: $(vertexsetsize(vset))")
println("Size of new Lindblad operator: $(size(lnonmoral))")


# Since the subspace has different dimension, one should consider using more advanced functions for initial states. The function below creates a state localized in subspace corresponding to vertex 1. Note it coincides with the first element of *vset*. 

rhoinit = nm_init(vset[[1]], vset)
println("state: $rhoinit")
println("equivalence class corresponding to vertex 1: $(vset[[1]])")


# This result in evolution, which prohibit passing amplitude to vertex 2. This can be seen in the example below.

rho = evolve(opnonmoral, rhoinit, timepoint)
println("Natural measurement on state: $(nm_measurement(rho, vset))")

