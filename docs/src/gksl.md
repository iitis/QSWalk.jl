```@meta
DocTestSetup = quote
   using QSWalk
end
```

## GKSL master equation

The GKSL master equation is general continuous-time open quantum evolution. The master equation base on to form of operators: the Hamiltonian, which describes the evolution of a closed system, and the Lindbladian operators which describes the evolution of an open system. Basic facts in the context of quantum stochastic walks can found [in this preprint](https://arxiv.org/abs/0905.2942) or [its published version](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.81.022323). For local interaction, where each Lindblad operator is a matrix with a single nonzero element, we provide a `local_lind` function which splits the matrix into mentioned operators.

Below we present documentation for essential functions used for simulating GKSL master equation.

```@index
Order = [:type, :function]
Modules = [QSWalk]
Pages   = ["gksl.md"]
```

## Full docs

```@docs
ket
bra
ketbra
proj
res
unres
evolve_generator(::AbstractMatrix{<:Number}, ::AbstractVector{<:AbstractMatrix{<:Number}}, ::AbstractMatrix{<:Number}, ::Real)
evolve_operator
evolve(::AbstractMatrix{<:Number}, ::AbstractMatrix{<:Number})
local_lind
```
