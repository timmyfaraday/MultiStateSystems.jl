# Network

## Quick Links

```@index
Pages = ["network.md"]
```

## Introduction

In its essence, a network (abbr: ntw) is a collection of sources, components and
users mapped onto a directed [multigraph](https://en.wikipedia.org/wiki/Multigraph).

Any network has eight attributes:
1. a multigraph `ntw.graph::DiMultigraph`
2. network properties `ntw.props::PropDict`
3. components `ntw.cmp::Vector{PropDict}`
4. sources `ntw.src::Vector{PropDict}`
5. users `ntw.usr::Vector{PropDict}`
6. component library `ntw.clib::LibDict`
7. source library `ntw.slib::LibDict`
8. users library `ntw.ulib::LibDict`

## Constructors

One constructor is implemented:
```@docs
MultiStateSystems.Network()
```

## Info

A number of properties are defined for networks captured by `NetworkInfo`.

* `solved`:
    the solved property describes whether the ugf of a network's users have been
    determined.

* `dependent_source`:
    the dependent source property flags that the output of all sources in the
    network are dependent on an uniform source, e.g., wind turbines in a single
    windfarm.

## Component

Any component (abbr: cmp) is represented by a `PropDict` mapped onto either the
vertices or edges of the multigraph `ntw.graph`. The collection of all
components is stored in `ntw.cmp`, which is a vector.

The link between the component id, i.e., its index in `ntw.cmp`, and the
vertex/edge of `ntw.graph` to which it is mapped is stored in `ntw.clib`
(key: vertex/edge, value: [cmp...]).

Any property may be added to the `Propdict` representing the component, however
certain properties are reserved for specific functionality of the tool, each
linked to a specific key [Symbol].

* `:node` is reserved for the component's node in the multigraph `ntw.graph`
* `:edge` is reserved for the component's edge in the multigraph `ntw.graph`
* `:std` is reserved for the state-transition diagram of the component

!!! note
    In order to link the specific component to either a vertex or an edge of the
    network, it is obligatory to provide them through the named argument `:node`
    or `:edge`, respectively.

Components may either be added to a network individually or grouped using,
respectively:

```@docs
add_component!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)
```

```@docs
add_components!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)
```

## Source

Any source (abbr: src) is represented by a `PropDict` mapped onto either the
vertices of the multigraph `ntw.graph`. The collection of all sources is stored
in `ntw.src`, which is a vector.

The link between the source id, i.e., its index in `ntw.src`, and the vertex of
`ntw.graph` to which it is mapped is stored in `ntw.slib`
(key: vertex, value: [cmp...]).

Any property may be added to the `Propdict` representing the source, however
certain properties are reserved for specific functionality of the tool, each
linked to a specific key [Symbol].

* `:node` is reserved for the source's node in the multigraph `ntw.graph`
* `:std` is reserved for the state-transition diagram of the source
* `:ntw` is reserved for the tuple (ntw,usr) representing the source, where usr is the user-id [Int] of the network ntw.
* `:dep_source` is reserved for the `dependent_source` property.
* `:dep_eval` is reserved for the `dependent_evaluation` property.

!!! note
    In order to link the specific source to a vertex of the network, it is
    obligatory to provide them through the named argument `:node`.

Sources may either be added to a network individually or grouped using,
respectively:

```@docs
add_source!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)
```

```@docs
add_sources!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)
```

## Users

Any user (abbr: usr) is represented by a `PropDict` mapped onto either the
vertices of the multigraph `ntw.graph`. The collection of all users is stored
in `ntw.usr`, which is a vector.

The link between the user id, i.e., its index in `ntw.usr`, and the vertex of
`ntw.graph` to which it is mapped is stored in `ntw.ulib`
(key: vertex, value: [cmp...]).

Any property may be added to the `Propdict` representing the user, however
certain properties are reserved for specific functionality of the tool, each
linked to a specific key [Symbol].

* `:node` is reserved for the user's vertex in the multigraph `ntw.graph`
* `:std` is reserved for the state-transition diagram of the user
* `:ind` is reserved for the user's indices

!!! note
    In order to link the specific user to a vertex of the network, it is
    obligatory to provide them through the named argument `:node`.

Users may either be added to a network individually or grouped using,
respectively:

```@docs
add_user!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)
```

```@docs
add_users!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)
```
