```@meta
DocTestSetup = quote
   using QSWalk
end
```

## GKSL master equation

GKSL master equation is general continuous-time open quantum evolution. The master equation base on to form of operators: the Hamiltonian, which describes the evolution of closed system, and the Lindbladian operators which describes the evolution of open system. Basic facts in context of quantum stochastic walks can found [here](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.81.022323). For local interaction, where each Lindblad operator is a matrix with single nonzero element, we created a `local_lind` function which splits matrix into mentioned operators.

Below we present a documentation of basic function used for simulating GKSL master equation.

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
