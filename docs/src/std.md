# State-Transition Diagram

## Quick Links

```@index
Pages = ["std.md"]
```

## Introduction

In its essence, a state-transition diagram (abbr: std) is a collections of
states and transitions mapped onto a
[directed graph](https://en.wikipedia.org/wiki/Directed_graph) through metadata
on its vertices and edges, respectively.

Any state-transition diagram has four attributes:
1. a directed graph `std.graph::SimpleDiGraph`
2. state-transition diagram properties `std.props::PropDict`
3. state properties `std.sprops::Dict{Int,PropDict}`
4. transition properties `std.tprops::Dict{Edge,PropDict}`

## Constructors

Two constructors are implemented:
```@docs
MultiStateSystems.STD()
```

```@docs
MultiStateSystems.STD(;prob::Array, kwargs...)
```

## Info

A number of properties are defined for state-transition diagrams, states and
transitions, respectively captured by (a) `STDInfo`, (b) `StateInfo` and (c)
`TransInfo`.

* `markovian` [(a),(c)]:
    the [Markov property](https://en.wikipedia.org/wiki/Markov_property) refers
    to the memoryless property of a stochastic process. It entails that the
    conditional probability distribution of future states only depends on the
    present state, including calendar time, however, they are independent of the
    present state's sojourn time (see: strong Markov property).

* `renewal` [(a), (b), (c)]:
    the renewal property entails that the stochastic process probabilistically start
    over at each arrival epoch. Consequently, a non-renewal transition entails that
    its to-state is not necessarily entered with a zero sojourn time.

* `time_homogeneous` [(a), (c)]:
    the [time-homogeneous property](https://stats.oecd.org/glossary/detail.asp?ID=3674)
    entails that the transition probability between two given states at any two
    times depends only on the difference between those times

* `trapping` [(b)]:
    the trapping property entails that a state is only partially/never exited upon
    entering.

* `solved` [(a)]:
    the solved property describes whether the state probabilities of a
    state-transition diagram have been determined.

## State

Any state is represented by a `PropDict` mapped onto the vertices of the
directed graph `std.graph`. The collection of all states is stored in
`std.sprops` which is a dictionary indexed using the vertex ids [Int].

Any property may be added the `PropDict` representing a state, however certain
properties are reserved for specific functionality of the tool, each linked to
a specific key [Symbol].

* `:info` is reserved for the StateInfo
* `:init` is reserved for the initial state probability ``p(0)`` [-]
* `:prob` is reserved for the state probability ``p(t)`` [-]
* `:φinit` is reserved for the initial state sojourn time ``φ(0)`` [hr]
* `:flow` is reserved for the state flow measure [m³/hr]
* `:power` is reserved for the state power measure [MW]

States may either be added to state-transition diagram individually or grouped using, respectively:

```@docs
MultiStateSystems.add_state!(std::MultiStateSystems.AbstractSTD; kwargs...)
```

```@docs
MultiStateSystems.add_states!(std::MultiStateSystems.AbstractSTD; kwargs...)
```

## Transition

Any transition is represented by a `PropDict` mapped onto the edges of the
directed graph `std.graph`. The collection of all transitions is stored in
`std.tprops` which is a dictionary indexed using the edge ids [Edge].

Any property may be added to the `PropDict` representing a transition, however
certain properties are reserved for specific functionality of the tool, each
linked to a specific key [Symbol].

* `:states` is reserved for the tuple (fr,to) of the from- and to-state
* `:rate` is reserved for the transition rate ``\rho(t,φ)`` [1/hr]
* `:distr` is reserved for the transition distribution
* `:type` is reserved for the transition type
    - `:f`   - failure
    - `:p`   - preventive maintenance decision
    - `:r`   - recoverability action
    - `:mcm` - minimal corrective maintenance
    - `:pcm` - perfect corrective maintenance
    - `:ppm` - perfect preventive maintenance

Transitions may either be added to a state-transition diagram individually or
grouped using, respectively:

```@docs
MultiStateSystems.add_transition!(std::MultiStateSystems.AbstractSTD; kwargs...)
```

```@docs
MultiStateSystems.add_transitions!(std::MultiStateSystems.AbstractSTD; kwargs...)
```
