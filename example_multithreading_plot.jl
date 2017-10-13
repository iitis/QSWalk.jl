# ------------------------------------------------------------------------------
# Case 5b: multithreading plot
# ------------------------------------------------------------------------------

using QSWalk
using JLD
using PyPlot

##

filename = "multithreading_data_global.jld"
data = load(filename)

ns = Float64[]
for n=data["linesizes"], m=1:data["repeat"]
  push!(ns, n)
end

for p=data["threadnumbers"]
  plot(ns, data["$p"], label="$p threads")
end
xlabel("n")
ylabel("t")

legend(loc="upper left",fancybox="true")
show()
savefig("time.pdf")
