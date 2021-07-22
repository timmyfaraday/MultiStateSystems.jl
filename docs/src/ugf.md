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
where $𝓞$ is the reduced state-space which only contains unique values for a
specific measure.

```@docs
MultiStateSystems.UGF
```

Alternatively, an ugf may be constructed based on a *solved* state-transition diagram:

```@docs
MultiStateSystems.UGF(msr::Symbol, std::MultiStateSystems.AbstractSTD)
```
