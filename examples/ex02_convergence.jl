#!/usr/bin/env python
# coding: utf-8

# # Convergence on local interaction model
# 

# ## Loading modules

using QSWalk
using LightGraphs # for graph functions 
using GraphPlot # for ploting graphs 
using LinearAlgebra # for linear algebra utilities


# ## Numerical proof of unique stationary state

# Basic parameters. We use Erdős–Rényi model to generate a directed graph. Strongly connected graphs have unique stationary state. Note the Hamiltonian is chosen to be the adjacency matrix of the underlying graph.

# number of nodes
dim = 10
# smaller vale of 'prob' can be used to geneate graphs which are not strongly connected
prob = 0.5
digraph = erdos_renyi(dim, prob, is_directed=true)
graph = Graph(digraph)

adj_digraph = Matrix(adjacency_matrix(digraph, dir=:in))
adj_graph = Matrix(adjacency_matrix(graph))
time = 100.

println("The graph is strongly connected: $(is_strongly_connected(digraph))")
gplot(digraph)


# Now we can calcuate the lindbladian and the subgroup generator.

lind = local_lind(adj_digraph)
evo_gen = evolve_generator(adj_graph, lind);


# The sufficient and necessary condition for the convergence of quantum stochastic evolution is that the null-space is  one-dimensional. Note for large matrices *eigs* may be a better option.

null_dim = count(x->abs(x)<1e-5, eigvals(evo_gen))
println("Dimensionality of null-space of the evolution operator: $null_dim")


# This allows efficient stationary state generation. Note that the trace may differ from one, as the eigenstate is normalized according to different norm.

eigendecomposition = eigen(evo_gen)
zeroindex = findfirst(x -> abs(x)<=1.e-5, eigendecomposition.values)
stationary_state = unres(vec(eigendecomposition.vectors[:, zeroindex]))

println("Trace of stationary state: $(sum(diag(stationary_state)))")
stationary_state /= sum(diag(stationary_state))
println("Trace of stationary state after the normalization: $(sum(diag(stationary_state)))")


# ## Convergence

# Since the stationary state is unique, all of states converge to it. We can show check this by taking three different states. Note, that for larger density states larger times of evolution might be required to achieve the convergence.

rhoinit1 = proj(1, dim)
rhoinit2 = proj(3, dim)
rhoinit3 = Diagonal(I,dim)/dim


# Since we apply the same evolution for all of the initial states, it is more efficient to calulate exponent once.

U = evolve_operator(evo_gen, time)
rho1 = evolve(U, rhoinit1)
rho2 = evolve(U, rhoinit2)
rho3 = evolve(U, rhoinit3);


# In order to show those states are essentialy the same we can calulate the norm of the difference.

println(norm(rho1-rho2))
println(norm(rho2-rho3))
println(norm(stationary_state-rho3))

