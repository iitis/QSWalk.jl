# ------------------------------------------------------------------------------
# Case 5a: multithreading
# ------------------------------------------------------------------------------

using QSWalk
using JLD
using PyPlot
## basic time measuring functions

function line_evolution_local(dim::Int, t::Real=1., ω::Real=0.5)
  adjacency = spdiagm((ones(dim-1),ones(dim-1)),(-1,1))
  lind = classical_lindblad_operators(adjacency)
  globaloperator = evolve_generator(adjacency, lind, ω)
  ρinit = proj(ceil(Int, dim/2), dim)
  evolve(globaloperator, ρinit, t)
end

function line_evolution_global(dim::Int, t::Real=1., ω::Real=0.5)
  adjacency = spdiagm((ones(dim-1),ones(dim-1)),(-1,1))
  globaloperator = evolve_generator(adjacency, [adjacency], ω)
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
    #addprocs(p)
    println("number of workers = ", nworkers())
    #eval(Expr(:toplevel, :(@everywhere using QSWalk)))

    for n=linesizes
      println("line size: $n")
      for m=1:repeat
        t = time_ns()
        evolvingfunction(n)
        push!(data["$p"], (time_ns()-t)/1.0e9) # nanosecond normalization
      end
    end

    #rmprocs(workers())
  end

  data
end

## data generation

linesizes = 100:100:1300
threadnumbers = [1,2,4,8] # collect(1:8)
filename = "multithreading_data_global_noaddproc.jld"
repeat = 1
t = 100.
ω = 0.5
funtionmeasured = (n->line_evolution_global(n, t, ω))

data = comptime(collect(linesizes), threadnumbers, repeat, funtionmeasured)
println("Calculations done")

data["linesizes"] = linesizes
data["repeat"] = repeat
data["function"] = "line_evolution_global"
data["threadnumbers"] = threadnumbers
data["function_parameters"] = Dict("time" => t, "ω"=>ω)

save(filename, data)
println("Data saved")
