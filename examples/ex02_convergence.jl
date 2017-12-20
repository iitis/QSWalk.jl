# ------------------------------------------------------------------------------
# Example 2: convergence for the local interaction regime
# ------------------------------------------------------------------------------

using QSWalk
using LightGraphs # for graph functions

## operator preparations
# number of nodes
dim = 10

# smaller vale of prob can be used to geneate graphs whcih are not strongly connected
prob = 0.5
digraph = erdos_renyi(dim, prob, is_directed=true)
graph = Graph(digraph)
adj_digraph = full(adjacency_matrix(digraph, dir=:in))
adj_graph = full(adjacency_matrix(graph))
time = 100.

lind = classical_lindbladian(adj_digraph)
evo_gen = evolve_generator(adj_graph, lind)

## dimensionality of the null space
println(count(x->abs(x)<1e-5, eigvals(evo_gen)))

## different initial states
rhoinit1 = proj(1, dim)
rhoinit2 = proj(3, dim)
rhoinit3 = eye(dim)/dim

U = evolve_operator(evo_gen, time)
rho1 = evolve(U, rhoinit1)
rho2 = evolve(U, rhoinit2)
rho3 = evolve(U, rhoinit3)

println(norm(rho1-rho2))
println(norm(rho2-rho3))
