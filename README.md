

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3457951.svg)](https://doi.org/10.5281/zenodo.3457951)


[![Build Status](https://travis-ci.org/iitis/QSWalk.jl.svg?branch=master)](https://travis-ci.org/iitis/QSWalk.jl)
[![Coverage Status](https://coveralls.io/repos/github/iitis/QSWalk.jl/badge.svg?branch=master)](https://coveralls.io/github/iitis/QSWalk.jl?branch=master)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://iitis.github.io/QSWalk.jl/latest)
[![DOI](https://zenodo.org/badge/100469695.svg)](https://zenodo.org/badge/latestdoi/100469695)


# QSWalk

## Package description

`QSWalk` provides a package for [Julia programming language](https://julialang.org/) which enables high-performance analysis of quantum stochastic walks. There are two main advantages of the presented packages over the existing software. First, it can be use to describe quantum stochastic walks in the local, as well as global regime. Second, it enables the user to seamlessly utilize parallel computing capabilities.

## Installation

The package can be installed simply with Pkg REPL:
```julia-repl
(v1.0) pkg> add QSWalk
```
## Examples

Examples can be found in `examples` subdirectory. They require [QSWalk](https://github.com/iitis/QSWalk.jl), [LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl), and [PyPlot](https://github.com/JuliaPy/PyPlot.jl).

## Citing QSWalk

The package is described in

Adam Glos, Jarosław Adam Miszczak, Mateusz Ostaszewski, *QSWalk.jl: Julia package for quantum stochastic walks analysis*, Computer Physics Communications (2018), DOI:[10.1016/j.cpc.2018.09.001](https://doi.org/10.1016/j.cpc.2018.09.001), [arXiv:1801.01294](https://arxiv.org/abs/1801.01294).

```tex
@article{GLOS2018,
title = "QSWalk.jl: Julia package for quantum stochastic walks analysis",
journal = "Computer Physics Communications",
year = "2018",
issn = "0010-4655",
doi = "https://doi.org/10.1016/j.cpc.2018.09.001",
author = "Glos, Adam and Miszczak, Jarosław Adam and Ostaszewski, Mateusz",
}
```
