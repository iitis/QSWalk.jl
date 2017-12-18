[![Build Status](https://travis-ci.org/ZKSI/QSWalk.jl.svg?branch=master)](https://travis-ci.org/ZKSI/QSWalk.jl)
[![Coverage Status](https://coveralls.io/repos/github/ZKSI/QSWalk.jl/badge.svg?branch=master)](https://coveralls.io/github/ZKSI/QSWalk.jl?branch=master)

# QSWalk

## Package description

QSWalk provides package for [Julia programming language](https://julialang.org/) which enables high-performance analysis of quantum stochastic walks. There are two main advantages of the presented packages over the existing software. First, it can be use to describe quantum stochastic walks in the local, as well as global regime. Second, it enables the user to seamlessly utilize parallel computing capabilities.

## Installation

QSWalk requires [Expokit package for Julia](https://github.com/acroy/Expokit.jl), implementing some routines contained in [EXPOKIT](http://www.maths.uq.edu.au/expokit). This package will be installed automatically with `QSWalk` installation

QSWalk can be installed directly form GitHub repository.

```julia
Pkg.clone("git://github.com/ZKSI/QSWalk.jl.git")
```

Package can be updated using ```Pkg.update``` function.

```julia
Pkg.update("QSWalk")
```
## Examples

Examples can be found in [QSWalkTutorials](https://github.com/ZKSI/QSWalkTutorials) project. Tutorials require [QSWalk](https://github.com/ZKSI/QSWalk.jl), [LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl), [PyPlot](https://github.com/JuliaPy/PyPlot.jl), [TikzGraphs](https://github.com/sisl/TikzGraphs.jl) modules.
