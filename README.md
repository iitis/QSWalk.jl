![Build Status](https://travis-ci.org/ZKSI/QSWalk.jl.svg?branch=master)](https://travis-ci.org/ZKSI/QSWalk.jl)

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
## Examples

There are some examples placed in the ```examples``` directory. Some of the
depand on external Julia packages (eg. ```PyPlot``` for ploting the results). To
run an example issue

```
julia examples/ex01_qsw_on_line.jl
``` 

Occasionaly Julia will complain about some depreciated functions. To get rid of
this type of warning use

```
julia --depwarn=no examples/ex01_qsw_on_line.jl
```
