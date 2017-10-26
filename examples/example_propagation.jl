# ------------------------------------------------------------------------------
# Case 1: propagation on path graph
# ------------------------------------------------------------------------------

using QSWalk
using PyPlot # for plot
using LightGraphs # for PathGraph

## operators
dim = 251 #odd for unique middle point
w = 0.5
timepoints = collect(0:2:100)
adjacency = adjacency_matrix(PathGraph(dim))
# alternatively adjacency = spdiagm((ones(dim-1),ones(dim-1)),(-1,1))

lindlocal = classical_lindblad_operators(adjacency)
midpoint = ceil(Int, dim/2)
opglobal = global_operator(adjacency, [adjacency], w)
oplocal = global_operator(adjacency, lindlocal, w)

## evolution
rhoglobal = evolve(opglobal, proj(midpoint, dim), timepoints)
rholocal = evolve(oplocal, proj(midpoint, dim), timepoints)

## second moment calculation
secmomentglobal = Float64[]
secmomentlocal = Float64[]
positions = (collect(1:dim)-midpoint)
for i=1:length(timepoints)
  push!(secmomentglobal, sum(positions.^2 .* diag(rhoglobal[i])))
  push!(secmomentlocal, sum(positions.^2 .* diag(rholocal[i])))
end

## plotting results

figure(figsize=[2.5, 1.5])
plot(timepoints, secmomentlocal, "k")
xlabel("t")
ylabel(L"\mu_2")
tick_params(labelsize=9)
axis([0, timepoints[end], 0, maximum(secmomentlocal)])
savefig("secondmomentlocal.pdf", bbox_inches="tight")
# show()

figure(figsize=[2.5,1.5])
plot(timepoints, secmomentglobal, "k")
xlabel("t")
ylabel(L"\mu_2")
tick_params(labelsize=9)
axis([0, timepoints[end], 0, maximum(secmomentglobal)])
savefig("secondmomentglobal.pdf", bbox_inches="tight")
# show()
##
