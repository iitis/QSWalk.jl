# QSWalk

## Package description

QSWalk provides package for [Julia programming language](https://julialang.org/) which enables high-performance analysis of quantum stochastic walks. There are two main advantages of the presented packages over the existing software. First, it can be use to describe quantum stochastic walks in the local, as well as global regime. Second, it enables the user to seamlessly utilize parallel computing capabilities.

## Installation

QSWalk requires [Expokit package for Julia](https://github.com/acroy/Expokit.jl), implementing some routines contained in [EXPOKIT](http://www.maths.uq.edu.au/expokit). This package can be installed as:

```julia
Pkg.clone("git://github.com/acroy/Expokit.jl.git")
```

QSWalk can be installed directly form GitHub repository.

```julia
Pkg.clone("git://github.com/ZKSI/QSWalk.jl.git")
```

Package can be updated using ```Pkg.update``` function.

```julia
Pkg.update("QSWalk")
```
