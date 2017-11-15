# ------------------------------------------------------------------------------
# Example 1: propagation on path graph
# ------------------------------------------------------------------------------

using QSWalk
using PyPlot # for plot
using LightGraphs # for PathGraph

## operators
dim = 251 # odd for unique middle point
w = 0.5
timepoints = collect(0:2:100)
adjacency = adjacency_matrix(PathGraph(dim))
# adjacency = spdiagm((ones(dim-1),ones(dim-1)),(-1,1)) # alternative version

lind_local = classical_lindbladian(adjacency)
midpoint = ceil(Int, dim/2)
op_global = evolve_generator(adjacency, [adjacency], w)
op_local = evolve_generator(adjacency, lind_local, w)

## evolution
rho_global = evolve(op_global, proj(midpoint, dim), timepoints)
rho_local = evolve(op_local, proj(midpoint, dim), timepoints)

## calculation of the second moment
secmoment_global = Float64[]
secmoment_local = Float64[]
positions = (collect(1:dim)-midpoint)
for i=1:length(timepoints)
   push!(secmoment_global, sum(positions.^2 .* diag(rho_global[i])))
   push!(secmoment_local, sum(positions.^2 .* diag(rho_local[i])))
end

## plotting the results

figure(figsize=[2.5, 1.5])
plot(timepoints, secmoment_local, "k")
xlabel("t")
ylabel(L"\mu_2")
tick_params(labelsize=9)
axis([0, timepoints[end], 0, maximum(secmoment_local)])
savefig("secondmomentlocal.pdf", bbox_inches="tight")
# show()

figure(figsize=[2.5,1.5])
plot(timepoints, secmoment_global, "k")
xlabel("t")
ylabel(L"\mu_2")
tick_params(labelsize=9)
axis([0, timepoints[end], 0, maximum(secmoment_global)])
savefig("secondmomentglobal.pdf", bbox_inches="tight")
# show()
##
