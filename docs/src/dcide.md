# DCIDE Integration

## Overview

The DCIDE integration module provides functionality to import, process, and solve electrical system models from DCIDE (DC Integrated Design Environment) JSON format. This integration enables seamless workflow between DCIDE electrical system design and MultiStateSystems.jl reliability analysis.

## JSON Format Support

The module supports DCIDE JSON files containing:
- **Component instances**: Electrical components with specifications and availability data
- **Connections**: Cable connections between components  
- **System topology**: Network structure definition

## Supported Component Types

### Sources (`dcide_src_types`)
- `"utility"`: Grid connection/utility supply
- `"battery"`: Battery energy storage systems
- `"PV"`: Photovoltaic solar panels

### Users/Loads (`dcide_usr_types`) 
- `"motor"`: Motor loads
- `"panel"`: Electrical panels/switchboards
- `"bus"`: Bus bars/distribution buses

### Components (`dcide_cmp_types`)
- `"panel"`: Electrical panels/switchboards
- `"switch"`: Switching devices
- `"protection"`: Protection devices (fuses, breakers)
- `"converter"`: Power electronic converters
- `"transformer"`: Transformers
- `"bus"`: Bus bars/distribution buses

### Special Component Types
- `"cable"`: Generated automatically for connections

## Data Structure

### Component Data Keys
Each component contains the following information:

```julia
dcide_cmp_keys = [:id, :name, :std, :type, :node, :edge, :bidirectional]
dcide_src_keys = [:id, :name, :std, :type, :node]  
dcide_usr_keys = [:id, :name, :type, :node]
```

Where:
- `:id`: Unique identifier in the system
- `:name`: Component name from DCIDE
- `:std`: State-transition diagram for reliability modeling
- `:type`: Component type classification
- `:node`: Node number for nodal components
- `:edge`: Edge tuple `(from, to)` for edge components  
- `:bidirectional`: Boolean for bidirectional components

## Availability Specifications

Components can specify availability data in two formats:

1. Direct Availability (`A`)
2. Failure/Repair Distributions
Example:
```json
"failure": {
    "distr": "Exponential",
    "λ": 0.001,
    "k": null },
"repair": {
    "distr": "LogNormal", 
    "μ": 24,
    "σ": 0.5 }

```

## Power Ratings

Power ratings are extracted from electrical specifications of the components in the DCIDE JSON when present.

## Main Functions

### `solve!(json::Dict{String,Any})`

Main function to solve a complete DCIDE system model.

**Process:**
1. Parse JSON and extract components (`init_elements`)
2. Connect components based on topology (`connect_elements!`)
3. Solve individual component reliability models
4. Build system network model (`build_network`)
5. Solve system-level reliability (`solve!(ntw)`)
6. Extend JSON with results (`extend_json!`)

**Usage:**
```julia
using JSON
using MultiStateSystems

# Load DCIDE JSON file
json_data = JSON.parsefile("system_model.json")

# Solve system reliability
solve!(json_data)

# Results are now embedded in json_data
```

## Examples

### Basic DCIDE System
```julia
using MultiStateSystems
using JSON

# Load and solve DCIDE model
dcide_model = JSON.parsefile("power_system.json")
solve!(dcide_model)

# Extract system availability
# Results are embedded back into the JSON structure
```

### Custom Component Reliability
```julia
# Components without availability specs get default perfect reliability
# Components with specs get appropriate STD models automatically
```