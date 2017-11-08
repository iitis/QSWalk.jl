# ------------------------------------------------------------------------------
# Case 4: convergence for local interaction
# ------------------------------------------------------------------------------

using QSWalk
using LightGraphs # for erdos_renyi

## operator preparations
dim = 10
digraph = erdos_renyi(dim, 0.5, is_directed=true)
graph = Graph(digraph)
adjacencydigraph = full(adjacency_matrix(digraph))
adjacencygraph = full(adjacency_matrix(graph))
time = 100.

lind = classical_lindblad_operators(adjacencydigraph)
globaloperator = evolve_generator(adjacencygraph, lind)

## dimensionality of null-space
println(count(x->abs(x)<1e-5, eigvals(globaloperator)))

## different initial states examples
rhoinit1 = proj(1, dim)
rhoinit2 = proj(3, dim)
rhoinit3 = eye(dim)/dim

U = evolve_operator(globaloperator, time)
rho1 = evolve(U, rhoinit1)
rho2 = evolve(U, rhoinit2)
rho3 = evolve(U, rhoinit3)

println(norm(rho1-rho2))
println(norm(rho2-rho3))
