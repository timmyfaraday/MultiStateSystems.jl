# Quick Start Guide

Welcome to **MultiStateSystems.jl**! This guide will get you up and running with reliability analysis of multi-state systems in just a few minutes.

## Installation

Install MultiStateSystems.jl using the Julia package manager:

```julia
using Pkg
Pkg.add("MultiStateSystems")
```

Or in the Julia REPL package mode:
```julia
] add MultiStateSystems
```

Then load the package:
```julia
using MultiStateSystems
```

## Core Concepts

MultiStateSystems.jl provides tools for modeling and analyzing the reliability of complex systems with multiple operational states. The main building blocks are:

- **State Transition Diagrams (STDs)**: Model component reliability behavior
- **Universal Generating Functions (UGFs)**: Represent system performance distributions  
- **Networks**: Connect components to model system topology
- **Distributions**: Define failure and repair time distributions
- **Stochastic Processes**: Solve reliability models over time

## Basic Workflow

### 1. Define Component Reliability with STDs

Create a simple two-state component (Available/Failed):

```julia
# Create a component with exponential failure and repair
component = STD(2)  # 2 states: Available, Failed

# Add states with power capacity and initial probabilities
add_states!(component, 
    name = ["Available", "Failed"],
    power = [100.0, 0.0],  # MW capacity in each state
    init = [1.0, 0.0]      # Start in Available state
)

# Add transitions with failure and repair rates
add_transitions!(component,
    states = [(1, 2), (2, 1)],  # Available→Failed, Failed→Available
    distr = [Exponential(8760.0), Exponential(24.0)]  # MTTF=1yr, MTTR=1day
)

# Solve the reliability model
solve!(component, SteadyStateProcess())
```

### 2. Work with Distributions

MultiStateSystems.jl supports various probability distributions:

```julia
# Exponential distribution (constant failure rate)
λ = Exponential(1000.0)  # Mean time = 1000 hours

# Weibull distribution (varying failure rate)  
weibull = Weibull(1000.0, 2.0)  # Scale=1000, Shape=2 (wear-out)

# Normal distribution
normal = Normal(100.0, 15.0)  # Mean=100, StdDev=15

# Log-normal distribution
lognormal = LogNormal(5.0, 0.5)  # Log-mean=5, Log-stddev=0.5

# Short notation available
dst = 𝑬(500.0)  # Same as Exponential(500.0)
```

### 3. Create Universal Generating Functions

UGFs represent performance probability distributions:

```julia
# Create UGF for a generator
generator = UGF(:power, 
    [0.0, 50.0, 100.0],     # Power levels (MW)
    [0.05, 0.15, 0.80]      # Probabilities
)

# Create UGF for a load
load = UGF(:demand,
    [30.0, 50.0, 70.0],     # Demand levels (MW)  
    [0.3, 0.5, 0.2]         # Probabilities
)
```

### 4. Build System Networks

Connect components to model system topology:

```julia
# Create a network
network = Network()

# Add sources (generators)
add_sources!(network,
    name = ["Gen1", "Gen2"],
    std = [gen1_std, gen2_std],
    node = [1, 2]
)

# Add users (loads)  
add_users!(network,
    name = ["Load1", "Load2"],
    type = ["critical", "non-critical"],
    node = [3, 4]
)

# Add components (transmission lines, transformers)
add_components!(network,
    name = ["Line1", "Transformer1"],
    std = [line_std, transformer_std],
    edge = [(1, 3), (2, 4)]
)

# Solve network reliability
solve!(network)
```

### 5. Calculate Reliability Indices

Extract reliability metrics from solved systems:

```julia
# For individual components
availability = component.props[:prob][1]  # Steady-state availability

# For network users (loads)
eens = EENS(network.usr[1])  # Expected Energy Not Served
gra = GRA(network.usr[1])    # Generation Ratio Availability

println("System availability: $(availability)")
println("EENS: $(eens) MWh/year")
println("GRA: $(gra)")
```

## Common Use Cases

### Power System Reliability

```julia
# Model a simple power system
network = Network()

# Add generation sources
add_sources!(network,
    name = ["Nuclear", "Gas", "Wind"],
    std = [nuclear_std, gas_std, wind_std],
    node = [1, 1, 2]  # Connected to buses 1 and 2
)

# Add transmission components
add_components!(network,
    name = ["Transmission_Line"],
    std = [line_std],
    edge = [(1, 3)]  # Connect bus 1 to load bus 3
)

# Add loads
add_users!(network,
    name = ["City_Load"],
    type = ["residential"],
    node = [3]
)

solve!(network)
```

### Manufacturing System

```julia
# Model manufacturing line reliability
production_line = STD(3)  # Normal, Degraded, Failed

add_states!(production_line,
    name = ["Normal", "Degraded", "Failed"],
    capacity = [1000.0, 500.0, 0.0],  # Units per hour
    init = [1.0, 0.0, 0.0]
)

# Add state transitions
add_transitions!(production_line,
    states = [(1, 2), (2, 3), (2, 1), (3, 1)],
    distr = [𝑬(720.0), 𝑬(168.0), 𝑬(8.0), 𝑬(4.0)]  # Hours
)

solve!(production_line, SteadyStateProcess())
```

### Series-Parallel System

```julia
# Model redundant system with 2-out-of-3 configuration
network = Network()

# Three parallel components
for i in 1:3
    component_std = STD(2)
    add_states!(component_std, 
        name = ["Up", "Down"],
        capacity = [1.0, 0.0],
        init = [1.0, 0.0]
    )
    add_transitions!(component_std,
        states = [(1, 2), (2, 1)],
        distr = [𝑬(8760.0), 𝑬(24.0)]
    )
    solve!(component_std, SteadyStateProcess())
    
    add_components!(network,
        name = ["Component_$i"],
        std = [component_std],
        edge = [(1, 2)]  # All in parallel
    )
end

solve!(network)
```

## Advanced Features

### Time-Dependent Analysis

```julia
# Use Markov processes for time-dependent analysis
solve!(component, MarkovProcess(t_end = 8760.0))  # Analyze over 1 year

# Semi-Markov processes for general distributions
solve!(component, SemiMarkovProcess(t_end = 8760.0))
```

### DCIDE Integration

```julia
# Load DCIDE JSON file for automatic model creation
using JSON
dcide_data = JSON.parsefile("system.json")
solve!(dcide_data)  # Automatically builds and solves the model
```

### Custom Performance Measures

```julia
# Define custom measures beyond power/capacity
availability_ugf = UGF(:availability, [0.0, 1.0], [0.05, 0.95])
cost_ugf = UGF(:cost, [0.0, 100.0, 200.0], [0.1, 0.6, 0.3])
```

## Tips and Best practices

1. **Start Simple**: Begin with basic two-state components before adding complexity
2. **Check Units**: Use consistent units throughout your model (e.g., all times in hours)
3. **Validate Models**: Compare results with analytical solutions for simple cases
4. **Documentation**: Use descriptive names for states and components
5. **Steady-State First**: Use `SteadyStateProcess()` for initial analysis before time-dependent methods

## Next Steps

- Explore the [API Reference](api_functions.md) for complete function documentation
- See [Network Manual](network.md) for advanced network modeling
- Check [Distributions Manual](distribution.md) for probability distribution details
- Review [Stochastic Processes](processes.md) for time-dependent analysis

## Getting Help

- Check the documentation for detailed explanations
- Look at example systems in the `examples/` directory
- Open issues on GitHub for bugs or feature requests
- Join discussions in the Julia community forums