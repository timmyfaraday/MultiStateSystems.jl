################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
# LVDC Short-Circuit Analysis Examples

This example demonstrates the capabilities of the MultiStateSystems.jl package
for analyzing short-circuit behavior and reliability in Low Voltage DC (LVDC) systems.

The example includes:
- Single bus network configuration
- Double bus network configuration  
- Various protection device types (SSCB, MCCB, Fuse)
- Short-circuit fault analysis
- System reliability calculations

## Usage

### Basic Usage (solves STDs from scratch and saves them):
```julia
# Run single bus example
run_single_bus_example()

# Run double bus example  
run_double_bus_example()

# Run complete analysis
run_complete_analysis()
```

### Advanced Usage (load pre-solved STDs):
```julia
# Generate STD files first
generate_example_stds()

# Then load from STDs (much faster)
run_single_bus_example(load_stds_from="MultiStateSystems.jl/examples/lvdc/ICDCM_2025 copy/results/single_bus_stds.jld2")
run_double_bus_example(load_stds_from="MultiStateSystems.jl/examples/lvdc/ICDCM_2025 copy/results/double_bus_stds.jld2")

# Demonstrate the loading approach
demo_load_from_stds()
```

### STD Management:
```julia
# Only save STDs (default behavior)
run_single_bus_example(save_stds=true)

# Don't save STDs
run_single_bus_example(save_stds=false)

# Load specific STD file
solved_stds = load_solved_stds("results/single_bus_stds.jld2")
```
"""

using MultiStateSystems
using Unitful
using JLD2

const _MSS = MultiStateSystems

# Load configuration and core functions
include("ICDCM_2025.jl")
include("config/system_parameters.jl")
include("src/fault_analysis.jl")
include("src/std_generation.jl") 
include("src/network_builder.jl")
include("src/reliability_analysis.jl")
include("utils/utilities.jl")
include("examples/single_bus_example.jl")
include("examples/double_bus_example.jl")

results = run_complete_analysis()

