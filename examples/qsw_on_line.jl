# ------------------------------------------------------------------------------
# Initialisation
# ------------------------------------------------------------------------------
using QSWalk

# Dimension of the system
dim = 51

# Initial position on the line
s0 = Int((dim+1)/2)

# Adjency matrix for the line
adjmtx = spdiagm((ones(dim-1),ones(dim-1)),(-1,1))
##

# ------------------------------------------------------------------------------
# Case 1: classical Lindbladian
# ------------------------------------------------------------------------------
ham = adjmtx
lin = classical_lindblad_operators(adjmtx)
evo = global_operator(ham, lin)
init = proj(s0, dim)
time_step = 1.0
time_points = collect(0:10)*time_step
##

# Versions using sparse evolution
@time evolve(evo, full(init), time_step)
@time evolve(evo, init, time_step)

# Versions using dense evolution
#@time evolve(full(evo), full(init), time_step)
#@time evolve(full(evo), init, time_step)

#evolve_operator(full(evo), time_step)
##

# ------------------------------------------------------------------------------
# Case 2: stochastic case with moralization
# ------------------------------------------------------------------------------

ham = adjmtx
lin = [adjmtx]
omg = 0.5

evo = global_operator(ham, lin, omg)

init = proj(s0, dim)
time_point = 1.
##

@time evolve(evo, init, time_point)

@time evolve(full(evo), init, time_point)

##

# ------------------------------------------------------------------------------
# Case 2: stochastic case with demoralization procedure
# ------------------------------------------------------------------------------

lin, vset = demoralized_lindbladian(adjmtx)
ham = (x-> if x!=0 1 else 0 end).(lin) # makes 0-1 matrix
ham_local = local_hamiltonian(vset)

omg = 0.5
evo = global_operator(ham, [lin], ham_local, omg)
init = initialstate([vset[s0]], vset)

timepoint = 1.
