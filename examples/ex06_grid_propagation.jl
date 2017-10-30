# ------------------------------------------------------------------------------
# Case 6: propagation on grid
# ------------------------------------------------------------------------------

using QSWalk
using LightGraphs # for grid graph
using Gadfly # for matrix plot

## operators
dim = 31 # needs to be odd for unique middle-point
t = dim/2.
w = 0.5
middlepoint = div(dim^2+1, 2)

g = Grid([dim, dim])
adj = adjacency_matrix(g)

initial_state = ketbra(middlepoint, middlepoint, dim^2)
F = global_operator(adj, [adj], w)

## evolution
ρ = evolve(F, initial_state, t)
probability = real.(diag(ρ))

## plot
spy(unres(sum.(probability)))
