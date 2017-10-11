using QSWalk
using PyPlot
using LightGraphs

###############################################################################
#scaling on line
##
#evolution calcularion

n = 251 #odd for unique middle point
ω = 0.5
timepoints = collect(0:1:100)
A = adjacency_matrix(PathGraph(n))

Llocal = classical_lindblad_operators(A)
midpoint = ceil(Int, n/2)
Fglobal = global_operator(A, [A], ω)
Flocal = global_operator(A, Llocal, ω)

ρglobal = evolve(Fglobal, proj(midpoint, n), timepoints)
ρlocal = evolve(Flocal, proj(midpoint, n), timepoints)

##
#second moment calculation
secmomentglobal = Float64[]
secmomentlocal = Float64[]
positions = (collect(1:n)-midpoint)
for i=1:length(timepoints)
  push!(secmomentglobal, sum(positions.^2 .* diag(ρglobal[i])))
  push!(secmomentlocal, sum(positions.^2 .* diag(ρlocal[i])))
end

##
#plots

figure(figsize=[2.5, 1.5])
plot(timepoints, secmomentlocal, "k")
xlabel("t")
ylabel(L"\mu_2")
tick_params(labelsize=9)
axis([0, timepoints[end], 0, maximum(secmomentlocal)])
savefig("pics/secondmomentlocal.pdf", bbox_inches="tight")

figure(figsize=[2.5,1.5])
plot(timepoints, secmomentglobal, "k")
xlabel("t")
ylabel(L"\mu_2")
tick_params(labelsize=9)
axis([0, timepoints[end], 0, maximum(secmomentglobal)])
savefig("pics/secondmomentglobal.pdf", bbox_inches="tight")
###############################################################################
#spontaneous moralization simple example
##
A = [0 0 0;
     0 0 0;
     1 1 0]
F = global_operator(zero(A), [A])
time = 100.

ρ = evolve(F, proj(1,3), time)
println(diag(ρ))

##
Ldemoralized, vset = nonmoralizing_lindbladian(L)
Hlocal = local_hamiltonian(vset)
Fdemoralized = global_operator(zero(Ldemoralized), [Ldemoralized], Hlocal)
ρinit = init_nonmoralized(vset[[1]], vset)
ρ = evolve(Fdemoralized, ρinit, time)
println(measurement_nonmoralized(ρ, vset))

###############################################################################
#line graph
##
n = 100 #odd for unique middle point
ω = 0.5
time = 40.
A = adjacency_matrix(PathGraph(n))

midpoint = ceil(Int, n/2)
L, vset = nonmoralizing_lindbladian(A)
H = global_hamiltonian(A)
Hlocal = local_hamiltonian(vset)
F = global_operator(H, [L], Hlocal, ω)
ρinit = init_nonmoralized(vset[[midpoint]], vset)

ρ = evolve(F, ρinit, time)
positions = (collect(1:n)-midpoint)
println(sum(positions .* measurement_nonmoralized(ρ, vset)))
##
n = 100 #odd for unique middle point
ω = 0.5
time = 40.
A = adjacency_matrix(PathGraph(n))

midpoint = ceil(Int, n/2)

lindbladians1 = Dict(1 => ones(1,1), 2 => [1 1; 1 -1])
lindbladians2 = Dict(1 => ones(1,1), 2 => [1 -1; 1 1])
L1, vset = nonmoralizing_lindbladian(A, lindbladians1)
L2, vset = nonmoralizing_lindbladian(A, lindbladians2)
H = global_hamiltonian(A)
Hlocal = local_hamiltonian(vset)
F = global_operator(H, [L1, L2], Hlocal, ω)
ρinit = init_nonmoralized(vset[[midpoint]], vset)

ρ = evolve(F, ρinit, time)
positions = (collect(1:n)-midpoint)
println(sum(positions .* measurement_nonmoralized(ρ, vset)))

###############################################################################
#convergence
n = 5
A = full(adjacency_matrix(PathGraph(n)))
t = 100.

L = classical_lindblad_operators(A)
F = global_operator(A, L)
println(filter(x->abs(x)<1e-5, eigvals(F)))

ρ1 = proj(1, n)
ρ2 = proj(3, n)
ρ3 = eye(n)/5

U = evolve_operator(F, time)
println(round.(evolve(U, ρ1), 2))
println(round.(evolve(U, ρ2), 2))
println(round.(evolve(U, ρ3), 2))
###############################################################################
