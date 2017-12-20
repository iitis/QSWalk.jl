
# coding: utf-8

# # Spontaneous moralization on path graph

# ## Loading modules

using QSWalk
using LightGraphs
using PyPlot


# ## Non-symmetric case

# Here we provide more advanced functions corresponding to sponatneous moralization. Below we start with simple parametrizations. Note *dim* should be odd for unique middle-point.

dim = 101
midpoint = ceil(Int, dim/2)
w = 0.5
timepoint = 40.
adjacency = adjacency_matrix(PathGraph(dim));


# We generate all of the operators needed for the evolution, including initial state.

lind, vset = nm_lindbladian(adjacency)
hglobal = nm_glob_ham(adjacency)
hlocal = nm_loc_ham(vset)
opnonsymmetric = evolve_generator(hglobal, [lind], hlocal, w)

rhoinit = nm_init(vset[[midpoint]], vset);


# Finally we make an evolution.

rho = evolve(opnonsymmetric, rhoinit, timepoint);


# Note that first moment of natural measurement deviates from zero.

positions = (collect(1:dim)-midpoint)
measurement_nonsymmetric = nm_measurement(rho, vset)
println("First moment centralized in midpoint: $(sum(positions .* nm_measurement(rho, vset)))")


# This is because of the non-symmetrices *lind* choice (analysis shows, that even removing *hlocal* and *hglobal* operators does not result in symmetric evolution). To confirm this, we plot natural measurement and its reflection around *n=midpoint*.

plot(positions, measurement_nonsymmetric, "k")
plot(positions, reverse(measurement_nonsymmetric), "b")
xlabel("position")
ylabel("probability")
axis([positions[1], positions[end], 0, maximum(measurement_nonsymmetric)])
vlines(0, 0., maximum(measurement_nonsymmetric), linestyles="--")


# ## Symmetric case

# The way to correct this is to choose another, symmetric Lindblad operator. While standard suage of *nm_lindbladian* will always output the same result, we can choose different basic orthogonal matrices to form different operators. In following example. We choose dictionary, which for different vertex degree chooses different matrix.

linddescription1 = Dict(1 => ones(1, 1), 2 => [1 1; 1 -1])
linddescription2 = Dict(1 => ones(1, 1), 2 => [1 1; -1 1])
lind1, vset = nm_lindbladian(adjacency, linddescription1)
lind2, vset = nm_lindbladian(adjacency, linddescription2);


# We can make similar creation for each vertex. For example one can choose. We restrict ourselves to *lind1* and *lind2*, as those guarantees symmetrization.

vset = make_vertex_set(adjacency)
linddescription3 = Dict(v=>rand(length(v), length(v)) for v = vlist(vset))
lind3, _ = nm_lindbladian(adjacency, linddescription3);


# Other functions should be adjusted to use both *lind1* and *lind2*.

hglobal = nm_glob_ham(adjacency)
hlocal = nm_loc_ham(vset)
opsymmetric = evolve_generator(hglobal, [lind1, lind2], hlocal, w)

rhoinit = nm_init(vset[[midpoint]], vset)
rho = evolve(opsymmetric, rhoinit, timepoint);


# Now both first momement and the distribution confirms symmetric evolution.

positions = (collect(1:dim)-midpoint)
measurement_symmetric = nm_measurement(rho, vset)
println("First moment centralized in midpoint: $(sum(positions .* measurement_symmetric))")

plot(positions, measurement_symmetric, "k")
plot(positions, reverse(measurement_symmetric), "b")
xlabel("position")
ylabel("probability")
axis([positions[1], positions[end], 0, maximum(measurement_symmetric)])
vlines(0, 0., maximum(measurement_symmetric), linestyles="--")

