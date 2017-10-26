using QSWalk
using LightGraphs #for grid graph
using Gadfly #for matrix plot

dim = 50
t = dim/4.
w = 0.5

g = Grid([dim, dim])
adj = adjacency_matrix(g)

initial_state = ketbra(1275,1275,dim^2)
F = global_operator(adj, [adj], w)

ρ = evolve(F, initial_state, t)
probability = real(diag(ρ))
lattice = unres(probability)

spy(lattice)
