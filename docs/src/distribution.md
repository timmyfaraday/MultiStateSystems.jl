# Distributions

## Introduction

The package includes a collection of probabilistic distributions and related
functions. The implementation is heavily influenced by the package
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl), but extends
it to allow for weighted distributions as well as the
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) parameters.

## Dirac Distribution

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

## LogNormal Distribution

```@docs
LogNormal(μ,σ,ω)
```
