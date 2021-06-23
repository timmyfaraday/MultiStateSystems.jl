# Distributions

## Introduction

The package includes a collection of probabilistic distributions and related
functions. The implementation is heavily influenced by the package
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl), but extends
it to allow for weighted distributions as well as the
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) parameters.

<<<<<<< HEAD
## Distributions

Each distribution has a number of constructors, both based on its full name and
abbreviation, e.g., `Weibull()` is equivalent to `𝑾()`. Furthermore,
additional constructors are included for limited input, where all remaining
parameters are set to their default, e.g., `Weibull(10.0u"hr")` is equivalent
to `Weibull(10.0u"hr",1.0,1.0)`.

All input parameters are either of the type `Number` or its subtype `Real`:
```
Number
 |- Unitful.Quantity, e.g., 10.0u"hr"
 |- Real
     |- AbstractFloat, e.g., 10.0
     |- Integer, e.g., 10
```

All distributions may be scaled using a weight parameter `ω`, where
`0.0 < ω ≤ 1.0`.

### Dirac Distribution
=======
## Dirac Distribution
>>>>>>> parent of 078bbdd (Full documentation for distributions.jl)

```@docs
Dirac(o,ω)
```

## Uniform Distribution

```@docs
Uniform(a,b,ω)
```

## Exponential Distribution

```@docs
Exponential(θ,ω)
```

## Weibull Distribution

```@docs
Weibull(θ,α,ω)
```

<<<<<<< HEAD
### Raised Cosine Distribution

```@docs
Cosine
```
=======
## LogNormal Distribution

```@docs
LogNormal(μ,σ,ω)
```
>>>>>>> parent of 078bbdd (Full documentation for distributions.jl)
