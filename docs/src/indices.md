# Indices

## Introduction

Indices attempt to quantitatively assess the reliability of a system using a single numerical value. 

## Expected Energy Not Served (EENS) [MWh]

This is the amount of electricity demand of a user, measured in MWh, that is
expected not to be met in a given year. 

```@docs
MultiStateSystems.EENS(usr::MultiStateSystems.PropDict)
```

## Generation Ratio Availability (GRA) [-]

This is the probability of at least transferring a specific percentage of the generation to a user through the network.

```@docs
MultiStateSystems.GRA(usr::MultiStateSystems.PropDict)
```