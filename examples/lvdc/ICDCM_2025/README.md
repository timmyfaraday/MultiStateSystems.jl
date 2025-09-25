# LVDC Short-Circuit Analysis Examples

This directory contains comprehensive examples demonstrating the capabilities of the **MultiStateSystems.jl** package for analyzing short-circuit behavior and reliability in Low Voltage DC (LVDC) distribution systems.

## Overview

The examples showcase:
- **Short-circuit fault current calculations** using RLC circuit analysis
- **Protection device modeling** (Solid State Circuit Breakers, MCCBs, Fuses)
- **State transition diagram generation** for system components
- **Network reliability analysis** for single and double bus configurations
- **Combined failure scenario analysis** for system resilience assessment

## Project Structure

```
ICDCM_2025/
├── main.jl                     # Main entry point
├── config/
│   └── system_parameters.jl    # System configuration and parameters
├── src/
│   ├── fault_analysis.jl       # DC fault current calculations
│   ├── std_generation.jl       # State transition diagram creation
│   ├── network_builder.jl      # Network topology construction
│   └── reliability_analysis.jl # Reliability metrics and analysis
├── examples/
│   ├── single_bus_example.jl   # Single bus network example
│   └── double_bus_example.jl   # Double bus network example
├── utils/
│   └── utilities.jl            # Helper functions and utilities
├── results/                    # Output directory for analysis results
└── README.md                   # This file
```

## Quick Start

### Running the Examples

```julia
# Load the package and run examples
include("main.jl")

# Run individual examples
single_bus_results = run_single_bus_example()
double_bus_results = run_double_bus_example()

# Run complete analysis
complete_results = run_complete_analysis()

# Demonstrate package capabilities
demonstrate_package_capabilities()
```

### Basic Usage

```julia
using MultiStateSystems
using Unitful

# Get default system configuration
sys_params = get_default_system_parameters()
int_params = get_default_integration_parameters()

# Run single bus analysis
results = single_bus_analysis(sys_params, int_params)

# Print summary
print_summary(results)
```

## System Configuration

The system is configured through structured parameters:

### Electrical Parameters
- **DC Voltage**: 750V nominal
- **Maximum Current**: 200A
- **Minimum Voltage**: 637.5V (85% of nominal)
- **Bus Capacitance**: 25mF (double bus configuration)
- **Cable Parameters**: 2.94µH/m inductance, 3.08mΩ/m resistance

### Protection Devices
- **SSCB**: Solid State Circuit Breaker (11µs mean clearing time)
- **MCCB**: Molded Case Circuit Breaker (5.1ms mean clearing time)  
- **Fuse**: Current-limiting fuse (I²t = 5000 A²s)

### Network Topologies
- **Single Bus**: Radial configuration with common DC bus
- **Double Bus**: Split configuration with optional bridge connection

## Analysis Features

### 1. Short-Circuit Analysis
- RLC circuit modeling for fault current calculation
- Underdamped and overdamped response analysis
- Critical time calculations (voltage drop, overcurrent, I²t energy)

### 2. Protection Device Modeling
- Circuit breaker clearing probability based on current magnitude
- Fuse clearing probability based on I²t energy integral
- Combined device coordination analysis

### 3. Reliability Analysis
- State transition diagrams with exponential, Weibull, and log-normal distributions
- Semi-Markov and Markov process modeling
- System availability and failure rate calculations

### 4. Network Analysis
- Single and double bus network topologies
- Cross-zone failure propagation analysis
- Bridge connection impact assessment

## Example Results

### Single Bus Network
Typical availability results for different protection devices:
- **SSCB**: >99.99% availability
- **MCCB**: >99.95% availability  
- **Fuse**: >99.90% availability

### Double Bus Network
- **Zone separation** reduces cross-zone failure propagation
- **Bridge connections** improve overall system resilience
- **Combined failure scenarios** demonstrate redundancy benefits

## Advanced Usage

### Custom Configuration

```julia
# Create custom system parameters
custom_params = SystemParameters(
    V_DC = 1000,           # 1kV system
    I_max = 500,           # 500A maximum current
    V_min = 850,           # 85% voltage limit
    I²t = 10000,           # Higher fuse rating
    L_p = 0.0,             # No parasitic inductance
    C_b = 50e-3,           # 50mF bus capacitance
    λᶜ = 0.0001u"1/yr/m",  # Higher cable failure rate
    P_zone = [0.99, 1.0]   # Zone protection factors
)

# Run analysis with custom parameters
results = single_bus_analysis(custom_params, int_params)
```

### Detailed Analysis

```julia
# Access detailed results
fault_probs = results["fault_probabilities"]
networks = results["networks"]
reliability = results["reliability_metrics"]

# Extract specific metrics
user_availability = reliability["user_availability"]
failure_rates = reliability["failure_rates_1.0"]

# Calculate performance metrics
performance = calculate_performance_metrics(user_availability)
```

### Exporting Results

```julia
# Save results to file
save_results(results, "results/my_analysis.jld2")

# Load previous results
previous_results = load_results("results/previous_analysis.jld2")

# Export to CSV (basic data only)
export_to_csv(performance[1.0]["SSCB"], "results/sscb_performance.csv")
```

## Mathematical Background

### Short-Circuit Current Calculation

The fault current is calculated using RLC circuit analysis:

**Underdamped case** (R² < 4L/C):
```
i(t) = (V_DC)/(ωL) * exp(-R/(2L) * t) * sin(ωt)
where ω = √(1/(LC) - (R/(2L))²)
```

**Overdamped case** (R² ≥ 4L/C):
```
i(t) = (V_DC)/(2ωL) * exp(-R/(2L) * t) * sinh(ωt)
where ω = √((R/(2L))² - 1/(LC))
```

### Protection Device Modeling

**Circuit Breaker Clearing Probability**:
- Modeled using log-normal distribution
- Parameters based on current magnitude and device characteristics
- Clearing time depends on fault severity

**Fuse Clearing Probability**:
- Based on I²t energy integral
- Time to clear calculated from energy accumulation
- Deterministic clearing once energy threshold reached

### State Transition Diagrams

Components modeled with multiple states:
- **A**: Available (normal operation)
- **V1**: Voltage fault (recoverable)
- **U2**: Short unavailable (corrective maintenance)
- **U3**: Long unavailable (replacement required)

Transitions use various distributions:
- **Exponential**: Short-circuit faults
- **Weibull**: Aging-related failures
- **Log-normal**: Repair and maintenance times

## Validation and Testing

The examples have been validated against:
- Industry standard LVDC protection coordination practices
- Published reliability data for power electronic components
- Fault current calculations verified with circuit simulation

## Contributing

To extend the examples:

1. **Add new protection devices**: Extend `get_protection_devices()` in `config/system_parameters.jl`
2. **Modify network topologies**: Update functions in `src/network_builder.jl`
3. **Add analysis metrics**: Extend `src/reliability_analysis.jl`
4. **Create new examples**: Add files to `examples/` directory

## References

1. "Reliability Assessment of Power Electronic Systems" - IEEE Standards
2. "LVDC Distribution Systems for DC Loads" - IEC Technical Specifications
3. "Multi-State System Reliability Analysis" - Academic Literature
4. "Protection Coordination in DC Distribution" - Industry Guidelines

## Dependencies

- **MultiStateSystems.jl**: Core framework for reliability analysis
- **Unitful.jl**: Physical units support
- **JLD2.jl**: Data serialization and storage
- **Statistics.jl**: Statistical functions
- **QuadGK.jl**: Numerical integration
- **Roots.jl**: Root finding algorithms

## License

Copyright 2025, Tom Van Acker, Glenn Emmers

This example is part of the MultiStateSystems.jl package and follows the same licensing terms.

---

For questions or support, please refer to the main MultiStateSystems.jl documentation or contact the package maintainers.
