using QSWalk
using PyPlot

n = 51 #needs to be odd
medium = Int((n+1)/2)
lineadjacencymatrix = spdiagm((ones(n-1),ones(n-1)),(-1,1))

##
#classical lindlbadian case
H = lineadjacencymatrix
Lclassical = QSWalk.classicallindbladoperators(lineadjacencymatrix)
F = globaloperator(H, Lclassical)
init = proj(medium, n)

timepoint = 1.
timepoints = collect(0:10)*1.

@time simpleevolve(F, init, timepoint)
Ffull = full(F)
@time simpleevolve(Ffull, init, timepoint)

simpleevolve(F, init, timepoints)

##
#stochastic moralized case
H = lineadjacencymatrix
L = [lineadjacencymatrix]
ω = 0.5
F = globaloperator(H, L, ω)
init = proj(medium, n)

timepoint = 1.

init = proj(medium, n)

@time simpleevolve(F, init, timepoint)
Ffull = full(F)
@time simpleevolve(Ffull, init, timepoint)

##
#stochastic demoralized case, see
#[1] Domino, K., Glos, A., & Ostaszewski, M. (2017). Spontaneous moralization
#problem in quantum stochastic walk. arXiv preprint arXiv:1701.04624.
L, vset = demoralizedlindbladian(lineadjacencymatrix)
Hlocal = localhamiltonian(vset)
ω = 0.5
H = (x-> if x!=0 1 else 0 end).(L) #makes 0-1 matrix

F = globaloperator(H, [L], Hlocal, ω)
init = initialstate([vset[medium]], vset)

timepoint = 1.

#if F is a sparse matrix, init needs to be full matrix
result = simpleevolve(F, full(init), timepoint)
distribution = distributionsummation(result, vset)
plot(distribution)
#the distribution is not symmetric due to theoretical reason, see
#[1] Domino, K., Glos, A., & Ostaszewski, M. (2017). Spontaneous moralization
#problem in quantum stochastic walk. arXiv preprint arXiv:1701.04624.
#below we present symmetrized version

##
#symmetrized stochastic demoralized case, see
#[1] Domino, K., Glos, A., & Ostaszewski, M. (2017). Spontaneous moralization
#problem in quantum stochastic walk. arXiv preprint arXiv:1701.04624.
L, vset = demoralizedlindbladian(lineadjacencymatrix,
                                 Dict(1=>ones(Complex128,1,1), 2=>[Complex128(1) 1; 1 -1]))
Lsym, _ = demoralizedlindbladian(lineadjacencymatrix,
                                 Dict(1=>ones(Complex128,1,1), 2=>[Complex128(1) 1; -1 1]))
Hlocal = localhamiltonian(vset)
ω = 0.5
H = (x-> if x!=0 1 else 0 end).(L) #makes 0-1 matrix

F = globaloperator(H, [L, Lsym], Hlocal, ω)
init = initialstate([vset[medium]], vset)

timepoint = 1.

#if F is a sparse matrix, init needs to be full matrix
result = simpleevolve(F, full(init), timepoint)
distribution = distributionsummation(result, vset)
plot(distribution)
