# ------------------------------------------------------------------------------
# Case 3: spontaneous moralization on path graph
# ------------------------------------------------------------------------------

using QSWalk
using LightGraphs # for PathGraph

## nonsymmetric case

dim = 100 #odd for unique middle point
w = 0.5
time = 40.
adjacency = adjacency_matrix(PathGraph(dim))

midpoint = ceil(Int, dim/2)
lind, vset = nonmoralizing_lindbladian(adjacency)
hglobal = global_hamiltonian(adjacency)
hlocal = local_hamiltonian(vset)
opnonsymmetric = global_operator(hglobal, [lind], hlocal, w)
rhoinit = init_nonmoralized(vset[[midpoint]], vset)

rho = evolve(opnonsymmetric, rhoinit, time)
positions = (collect(1:dim)-midpoint)
println(sum(positions .* measurement_nonmoralized(rho, vset)))

## symmetric case

dim = 100 #odd for unique middle point
w = 0.5
time = 40.
adjacency = adjacency_matrix(PathGraph(dim))
midpoint = ceil(Int, dim/2)

linddescription1 = Dict(1 => ones(1,1), 2 => [1 1; 1 -1])
linddescription2 = Dict(1 => ones(1,1), 2 => [1 -1; 1 1])
lind1, vset = nonmoralizing_lindbladian(adjacency, linddescription1)
lind2, vset = nonmoralizing_lindbladian(adjacency, linddescription2)
hglobal = global_hamiltonian(adjacency)
hlocal = local_hamiltonian(vset)
opsymmetric = global_operator(hglobal, [lind1, lind2], hlocal, w)
rhoinit = init_nonmoralized(vset[[midpoint]], vset)

rho = evolve(opsymmetric, rhoinit, time)
positions = (collect(1:dim)-midpoint)
println(sum(positions .* measurement_nonmoralized(rho, vset)))
