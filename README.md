[![Build Status](https://travis-ci.org/QuantumWalks/QSWalk.jl.svg?branch=master)](https://travis-ci.org/QuantumWalks/QSWalk.jl)
[![Coverage Status](https://coveralls.io/repos/github/QuantumWalks/QSWalk.jl/badge.svg?branch=master)](https://coveralls.io/github/QuantumWalks/QSWalk.jl?branch=master)

# QSWalk

## Package description

QSWalk provides package for [Julia programming language](https://julialang.org/) which enables high-performance analysis of quantum stochastic walks. There are two main advantages of the presented packages over the existing software. First, it can be use to describe quantum stochastic walks in the local, as well as global regime. Second, it enables the user to seamlessly utilize parallel computing capabilities.

## Installation

QSWalk requires [Expokit package for Julia](https://github.com/acroy/Expokit.jl), implementing some routines contained in [EXPOKIT](http://www.maths.uq.edu.au/expokit). This package will be installed automatically with `QSWalk` installation

QSWalk can be installed directly form GitHub repository.

```julia
Pkg.add("QSWalk")
```

Package can be updated using ```Pkg.update``` function.

```julia
Pkg.update("QSWalk")
```
## Examples

Examples can be found in `examples` subdirectory. They require [QSWalk](https://github.com/QuantumWalks/QSWalk.jl), [LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl), [PyPlot](https://github.com/JuliaPy/PyPlot.jl), [TikzGraphs](https://github.com/sisl/TikzGraphs.jl) modules.

## Citing QSWalk

The package is described in

Adam Glos, Jarosław Adam Miszczak, Mateusz Ostaszewski, *QSWalk.jl: Julia package for quantum stochastic walks analysis*, [arXiv:1801.01294](https://arxiv.org/abs/1801.01294) (2018).
