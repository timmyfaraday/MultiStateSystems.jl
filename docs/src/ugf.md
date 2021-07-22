# Universal Generating Function

## Introduction

An universal generating function (abbr: ugf) approach approach enables the
combination of components and sources to evaluate their performance with respect
to a network's users.

An universal generator function $ω_x(z)$ of a component/source $x$ combines the
probabilities `p` and associated values `v` of a specific measure `m` in one
polynomial expression:

```math
    ω_{x}(z) = \sum_{o \in \mathcal{O}} p_o \cdot z^{v_o}
```
where 𝓞 is the reduced state-space which only contains unique values for a
specific measure.

## Constructors

Two constructors are implemented:
```@docs
    UGF(msr::Symbol, val::Vector, prb::Vector)

A UGF constructor for a specific measure `msr` based on a given probability
vector `prb` and value vector `val`.

This function automatically reduces the state-space to where only unique values
and associated probabilites remain.

# Example
```julia-repl
julia> ugfᵍᵉⁿ = UGF(:flow, [0.0u"MW",0.0u"MW",2.0u"MW"], [0.1,0.2,0.7])
julia> isequal(ugfᵍᵉⁿ,[0.0u"MW",2.0u"MW"])
true
```
```

```@docs
MultiStateSystems.UGF(msr::Symbol, std::MultiStateSystems.AbstractSTD)
```
