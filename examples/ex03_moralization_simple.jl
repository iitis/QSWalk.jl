# ------------------------------------------------------------------------------
# Case 2: spontaneous moralization on simple graph
# ------------------------------------------------------------------------------

using QSWalk

## operators
adjacency = [0 0 0;
             0 0 0;
             1 1 0]
opmoral = evolve_generator(zero(adjacency), [adjacency])
time = 100.

## evolution
rho = evolve(opmoral, proj(1,3), time)
println(diag(rho))

## nonmoralizing evolution case
lnonmoral, vset = nonmoralizing_lindbladian(adjacency)
hlocal = local_hamiltonian(vset)
opnonmoral = evolve_generator(zero(lnonmoral), [lnonmoral], hlocal)

println(vertexsetsize(vset))
println(vset)

rhoinit = init_nonmoralized(vset[[1]], vset)
rho = evolve(opnonmoral, rhoinit, time)
println(measurement_nonmoralized(rho, vset))
