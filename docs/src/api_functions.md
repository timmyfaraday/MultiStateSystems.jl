# MultiStateSystems.jl API Documentation Summary

This file provides an overview of the documented functions in the MultiStateSystems.jl package.

## Exported Functions and Types

### Constants
- `BASE_DIR` - Base directory path of the package

### Distribution Types and Constructors
- `Exponential()` / `𝑬()` - Exponential distribution constructors
- `LogNormal()` / `𝑳()` - Log-normal distribution constructors  
- `Normal()` / `𝑵()` - Normal (Gaussian) distribution constructors
- `Weibull()` / `𝑾()` - Weibull distribution constructors

### Universal Generating Functions
- `UGF()` - Universal generating function constructor

### Network Types and Functions
- `Network()` - Network constructor
- `add_component!()` - Add single component to network
- `add_components!()` - Add multiple components to network
- `add_bidirectional_component!()` - Add single bidirectional component
- `add_bidirectional_components!()` - Add multiple bidirectional components
- `add_source!()` - Add single source to network
- `add_sources!()` - Add multiple sources to network
- `add_user!()` - Add single user to network
- `add_users!()` - Add multiple users to network

### State Transition Diagrams
- `STD()` - State transition diagram constructor
- `solvedSTD()` - Solved STD constructor with probabilities
- `add_state!()` - Add single state to STD
- `add_states!()` - Add multiple states to STD
- `add_transition!()` - Add single transition to STD
- `add_transitions!()` - Add multiple transitions to STD

### Stochastic Process Types
- `SteadyStateProcess` - Steady-state Markov process
- `MarkovProcess` - Time-dependent Markov process
- `SemiMarkovProcess` - Semi-Markov process with general distributions

### Solver Functions
- `solve!()` - Main solver function (works with STD and Network types)
- `state_conv()` - Convolution function for semi-Markov processes

### Reliability Indices (for Users)
- `EENS()` - Expected Energy Not Served
- `GRA()` - Generation Ratio Availability