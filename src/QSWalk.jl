module QSWalk

using Expokit
using SparseArrays
using LinearAlgebra

include("utils.jl")
include("dirac.jl")
include("operator.jl")
include("demoralization.jl")
include("evolution.jl")

"""
QSWalk provides package for Julia programming language which enables
high-performance analysis of quantum stochastic walks. There are two main
advantages of the presented packages over the existing software. First, it can
be use to describe quantum stochastic walks in the local, as well as global
regime. Second, it enables the user to seamlessly utilize parallel computing
capabilities.

"""
QSWalk

end
