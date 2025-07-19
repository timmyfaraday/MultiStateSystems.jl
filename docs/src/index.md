# MultiStateSystems.jl

## Overview

MultiStateSystems.jl is a Julia package (v1.8+) to solve multi-state systems,
often found in reliability engineering.

## Installation

The latest stable release of MultiStateSystems can be installed using the Julia
package manager:

```
] add MultiStateSystems
```

In order to test whether the package works, run:

```
] test MultiStateSystems
```

## Documentation Structure

This documentation is organized as follows:
- **Getting Started**: Quick introduction and basic examples
- **DSL Manual**: Domain-specific language components (STDs, distributions, networks)
- **Models Manual**: Advanced modeling features (dependence, stochastic processes)
- **Integrations**: External tool integrations (DCIDE integration)
- **Output**: Reliability indices and result interpretation
- **API Reference**: Complete function reference with examples

## Acknowledgements

The primary developer is Tom Van Acker ([@timmyfaraday](https://github.com/timmyfaraday)), Elia, 
with support from the following contributors:
- Glenn Emmers ([@glenn-sergej](https://github.com/Glenn-sergej)), KU Leuven, semi-Markov implementation, DCIDE implementation, LVDC examples
- Gayan Abeynayake ([@gayan86](https://github.com/gayan86)), Ørsted, Anholt wind farm example

## License

This code is provided under a BSD license.