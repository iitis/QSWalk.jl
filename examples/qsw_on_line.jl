using QSWalk

# Dimension of the system
dim = 51

# Initial position on the line
s0 = Int((dim+1)/2)

# Adjency matrix for the line
adjmtx = spdiagm((ones(dim-1),ones(dim-1)),(-1,1))

# Case 1: classical Lindbladian
ham = adjmtx
lin = classicallindbladoperators(adjmtx)
evo = globaloperator(ham, lin)
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

evolve_operator(full(evo), time_step)
