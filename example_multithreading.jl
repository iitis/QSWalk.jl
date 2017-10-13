# ------------------------------------------------------------------------------
# Case 5: multithreading
# ------------------------------------------------------------------------------

using QSWalk
using JLD
using PyPlot
## basic time measuring functions
function line_evolution_local(dim::Int, t::Real=1., ω::Real=0.5)
  adjacency = spdiagm((ones(dim-1),ones(dim-1)),(-1,1))
  lind = classical_lindblad_operators(adjacency)
  globaloperator = global_operator(adjacency, lind, ω)
  ρinit = proj(ceil(Int, dim/2), dim)
  evolve(globaloperator, ρinit, t)
end

function line_evolution_global(dim::Int, t::Real=1., ω::Real=0.5)
  adjacency = spdiagm((ones(dim-1),ones(dim-1)),(-1,1))
  globaloperator = global_operator(adjacency, [adjacency], ω)
  ρinit = proj(ceil(Int, dim/2), dim)
  evolve(globaloperator, ρinit, t)
end

function comptime(linesizes::Vector{Int},
                                 threadnumbers::Vector{Int},
                                 repeat::Int,
                                 evolvingfunction::Function)
  data = Dict()

  # first run
  line_evolution_global(10)

  # the calculation
  for p=threadnumbers
    data["$p"] = Float64[]
    BLAS.set_num_threads(p)
    addprocs(p)
    println("number of workers = ", nworkers())
    eval(Expr(:toplevel, :(@everywhere using QSWalk)))

    for n=linesizes, m=1:repeat
      t = time_ns()
      evolvingfunction(n)
      push!(data["$p"], (time_ns()-t)/1.0e9) # nanosecond normalization
    end
    rmprocs()
  end

  data
end

## data generation

linesizes = 1000:100:1300
threadnumbers = [4]
filename = "multithreading_data_global.jld"
repeat = 1
t = 100.
ω = 0.5
funtionmeasured = (n->line_evolution_global(n, t, ω))

data = comptime(collect(linesizes), threadnumbers, repeat, funtionmeasured)


data["linesizes"] = linesizes
data["repeat"] = repeat
data["function"] = "line_evolution_global"
data["threadnumbers"] = threadnumbers
data["function_parameters"] = Dict("time" => t, "ω"=>ω)

save(filename, data)

## plotting

ns = Float64[]
for n=data["linesizes"], m=1:d["repeat"]
  push!(ns, n)
end

for p=data["threadnumbers"]
  plot(ns, data["$p"], label="$p threads")
end

legend(loc="upper left",fancybox="true")
show()
