################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models, often found in           #
# reliability engineering.                                                     #
# See https://github.com/timmyfaraday/MultiStateSystems.jl                     #
################################################################################
# Authors: Tom Van Acker, Glenn Emmers                                         #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

# constants ####################################################################
const dcide_bus_types = ["bus", "panel"]
const dcide_cmp_types = ["panel", "switch", "protection", "converter", "transformer", "bus"]
const dcide_src_types = ["utility", "battery", "PV"]
const dcide_usr_types = ["motor"]

const dcide_cmp_keys = [:id, :name, :std, :type, :node, :edge, :bidirectional]
const dcide_src_keys = [:id, :name, :std, :type, :node]
const dcide_usr_keys = [:id, :name, :type, :node]

# functions ####################################################################
## filter
"""
    filter_nodal_cmp(cmp)

Filter components to return only nodal components (components connected to network nodes).
Excludes edge information and returns components that have non-nothing node assignments.

# Arguments
- `cmp`: Component data structure with node and edge information

# Returns
- Named tuple containing filtered component data excluding edge information
"""
filter_nodal_cmp(cmp) =
    (;zip(  [key for key in keys(cmp) if key ≠ :edge],
            [cmp[key][.!isnothing.(cmp.node)] for key in keys(cmp) if key ≠ :edge])...)

"""
    filter_bidirectional_edge_cmp(cmp)

Filter components to return only bidirectional edge components (e.g., cables, switches).
Returns components that have edge connections and are marked as bidirectional.

# Arguments
- `cmp`: Component data structure with edge and bidirectional information

# Returns
- Named tuple containing filtered bidirectional edge components
"""
filter_bidirectional_edge_cmp(cmp) =
    (;zip(  [key for key in keys(cmp) if key ≠ :node],
            [cmp[key][.!isnothing.(cmp.edge) .&& cmp.bidirectional] for key in keys(cmp) if key ≠ :node])...)

"""
    filter_unidirectional_edge_cmp(cmp)

Filter components to return only unidirectional edge components (e.g., diodes, transformers).
Returns components that have edge connections but are not bidirectional.

# Arguments
- `cmp`: Component data structure with edge and bidirectional information

# Returns
- Named tuple containing filtered unidirectional edge components
"""
filter_unidirectional_edge_cmp(cmp) =
    (;zip(  [key for key in keys(cmp) if key ≠ :node],
            [cmp[key][.!isnothing.(cmp.edge) .&& .!cmp.bidirectional] for key in keys(cmp) if key ≠ :node])...)

"""
    filter_panel_cmp(cmp)

Filter components to return only panel and bus type components.
These components act as connection points in the electrical network.

# Arguments
- `cmp`: Component data structure with type and node information

# Returns
- Named tuple containing filtered panel/bus components
"""
filter_panel_cmp(cmp) =
    (;zip(  [key for key in keys(cmp) if key ≠ :edge],
            [cmp[key][.!isnothing.(cmp.node) .&& in.(cmp.type, Ref(dcide_bus_types))] for key in keys(cmp) if key ≠ :edge])...)

## std
"""
    get_std(component)

Extract or create a State Transition Diagram (STD) from a DCIDE component specification.

This function parses the availability specifications from a DCIDE component and creates
an appropriate STD model. If no availability data is provided, returns a perfect 
reliability model (always available with infinite capacity).

# Arguments
- `component`: DCIDE component dictionary containing specifications

# Returns
- `STD` or `solvedSTD`: State transition diagram representing the component's reliability behavior

# Examples
```julia
# Component with availability specification
component = Dict("specifications" => Dict("availability" => Dict(:A => 0.99)))
std = get_std(component)

# Component without availability data (perfect reliability)
basic_component = Dict()
std = get_std(basic_component)  # Returns perfect reliability STD
```

# Notes
- If component has direct availability (:A), creates a two-state solved STD
- If component has failure/repair distributions, creates an unsolved STD with transitions
- Power ratings are extracted from electrical port specifications when available
- Defaults to infinite power capacity if no electrical specifications found
"""
function get_std(component)
    specs = get(component, "specifications", nothing)
    if specs === nothing || !haskey(specs, "availability")
        return solvedSTD(prob = [1], power = [(Inf)u"MW"])
    end

    ports = get(get(specs, "electrical", Dict()), "ports", [])
    power = begin
        if !isempty(ports) && !isempty(ports[1])
            ac = get(get(get(ports[1], "AC", Dict()), "power", Dict()), "nom", nothing)
            dc = get(get(get(ports[1], "DC", Dict()), "power", Dict()), "nom", nothing)
            ac !== nothing ? [ac * u"W", 0.0u"W"] :
            dc !== nothing ? [dc * u"W", 0.0u"W"] :
            [Inf * u"W", 0.0u"W"]
        else
            [Inf * u"W", 0.0u"W"]
        end
    end

    av = specs["availability"]
    if get(av, :A, nothing) !== nothing
        return solvedSTD(prob = [av[:A], 1 - av[:A]], power = power)
    end

    fail = av[:failure]
    rep = av[:repair]
    distr = [
        fail[:k] === nothing ? fail[:distr](fail[:λ]) : fail[:distr](fail[:λ], fail[:k]),
        rep[:σ] === nothing ? rep[:distr](rep[:μ]) : rep[:distr](rep[:μ], rep[:σ])
    ]

    std = STD()
    add_states!(std, name=["Available", "Unavailable"], power=power, init=[1.0, 0.0])
    add_transitions!(std, states=[(1, 2), (2, 1)], distr=distr)
    return std
end
# TODO - add λ, μ, nom. power, etc.

## init 
"""
    init_elements(json)

Initialize component, source, and user elements from DCIDE JSON specification.

Parses the DCIDE JSON structure and categorizes all elements into three main types:
- Components (cmp): Electrical components like switches, converters, buses
- Sources (src): Power generation sources like utilities, batteries, PV systems  
- Users (usr): Electrical loads like motors

# Arguments
- `json`: DCIDE JSON dictionary containing componentInstances and connections

# Returns
- `Tuple{NamedTuple, NamedTuple, NamedTuple}`: (cmp, src, usr) containing categorized elements

# Examples
```julia
# Load DCIDE JSON file
json_data = JSON.parsefile("dcide_system.json")
cmp, src, usr = init_elements(json_data)

# Access component names
println(cmp.name)  # Vector of component names
println(src.type)  # Vector of source types
```

# Notes
- Automatically determines component topology (nodal vs edge) based on port count
- Creates STDs for each component using get_std()
- Assigns initial node/edge values (0 for unconnected, nothing for N/A)
- Processes both componentInstances and connections from JSON
"""
function init_elements(json)
    # Initialize the named tuples to store the components
    cmp = (; zip(dcide_cmp_keys, [[] for key in dcide_cmp_keys])...)
    src = (; zip(dcide_src_keys, [[] for key in dcide_src_keys])...)
    usr = (; zip(dcide_usr_keys, [[] for key in dcide_usr_keys])...)

    for (nc,component) in enumerate(json["componentInstances"])
        if component["componentType"] in dcide_src_types
            for key in dcide_src_keys
                if      key == :id              push!(src[key], nc)
                elseif  key == :name            push!(src[key], component["id"])
                elseif  key == :std             push!(src[key], get_std(component))
                elseif  key == :type            push!(src[key], component["componentType"])
                elseif  key == :node            push!(src[key], 0)
                else    println("ERROR: the key $key is not supported for src yet")
            end end
        elseif component["componentType"] in dcide_usr_types
            for key in dcide_usr_keys
                if      key == :id              push!(usr[key], nc)
                elseif  key == :name            push!(usr[key], component["id"])
                elseif  key == :type            push!(usr[key], component["componentType"])
                elseif  key == :node            push!(usr[key], 0)
                else    println("ERROR: the key $key is not supported for usr yet")
            end end
        elseif component["componentType"] in dcide_cmp_types
            for key in dcide_cmp_keys
                if      key == :id              push!(cmp[key], nc)
                elseif  key == :name            push!(cmp[key], component["id"])
                elseif  key == :std             push!(cmp[key], get_std(component))
                elseif  key == :type            push!(cmp[key], component["componentType"])
                elseif  key == :node
                    if length(component["configuration"]["ports"]) == 2
                                                push!(cmp[key], nothing)
                    else
                                                push!(cmp[key], 0)
                    end
                elseif  key == :edge
                    if length(component["configuration"]["ports"]) == 2
                                                push!(cmp[key], (0,0))
                    else
                                                push!(cmp[key], nothing)
                    end
                elseif  key == :bidirectional   push!(cmp[key], true)
    end end end end

    for (nc, connection) in enumerate(json["connections"]), key in dcide_cmp_keys
        if      key == :id              push!(cmp[key], nc)
        elseif  key == :name            push!(cmp[key], connection["id"])
        elseif  key == :std             push!(cmp[key], get_std(connection))
        elseif  key == :type            push!(cmp[key], "cable")
        elseif  key == :node            push!(cmp[key], nothing)
        elseif  key == :edge            push!(cmp[key], (0,0))
        elseif  key == :bidirectional   push!(cmp[key], true)
    end end

    return cmp, src, usr
end

## connect
"""
    connect_elements!(cmp, src, usr, json)

Establish electrical connections between components based on DCIDE connection specifications.

This function processes the connections defined in the DCIDE JSON and assigns proper
node numbers to create the network topology. It handles different connection scenarios:
- Direct connections to sources and users (single node assignment)
- Connections through edge components (cable connections)
- Multi-port nodal components (buses, panels)

# Arguments
- `cmp`: Component data structure (modified in-place)
- `src`: Source data structure (modified in-place)  
- `usr`: User data structure (modified in-place)
- `json`: DCIDE JSON dictionary containing connections

# Side Effects
- Modifies the node and edge fields of cmp, src, and usr structures
- Assigns unique node numbers for network topology
- Updates edge tuples for cable and connection components

# Algorithm
1. Iterate through each connection in JSON
2. Assign node numbers to 'from' and 'to' components
3. Handle special cases for multi-port components (buses/panels)
4. Update edge information for cable connections
5. Increment node counter appropriately

# Notes
- Node numbering starts from 1
- Sources and users get single node assignments
- Edge components get (from_node, to_node) tuples
- Nodal components can connect to multiple nodes
"""
function connect_elements!(cmp, src, usr, json)
    node_cntr = 1
    stop_cntr = false
    
    for conn in json["connections"]
        # get the fr- and to- elements of a cable connection
        name = conn["id"]
        fr, to = conn["from"]["componentInstanceId"], conn["to"]["componentInstanceId"]
    
        # any src, usr or one side of an edge-cmp can only be associated with one connection
        src.node[src.name .== fr] .= node_cntr
        usr.node[usr.name .== fr] .= node_cntr
        if !isempty(cmp.edge[.!isnothing.(cmp.edge) .&& cmp.name .== fr])
            idx = findfirst(.!isnothing.(cmp.edge) .&& cmp.name .== fr)
            cmp.edge[idx] = (cmp.edge[idx][1], node_cntr)
        end
    
        # a nodal-cmp can be associated with multiple connections
        if !isempty(cmp.node[.!isnothing.(cmp.node) .&& cmp.name .== fr])
            cmp_node = cmp.node[.!isnothing.(cmp.node) .&& cmp.name .== fr][1]
            if iszero(cmp_node)
                cmp.node[.!isnothing.(cmp.node) .&& cmp.name .== fr] .= node_cntr
            else
                idx = findfirst(cmp.name .== name)
                cmp.edge[idx] = (cmp_node, cmp.edge[idx][2])
                stop_cntr = true
            end
        end

        if !stop_cntr
            idx = findfirst(cmp.name .== name)
            cmp.edge[idx] = (node_cntr, cmp.edge[idx][2])
        end

        node_cntr += stop_cntr ? 0 : 1
        stop_cntr = false

        # any src, usr or one side of an edge-cmp can only be associated with one connection
        src.node[src.name .== to] .= node_cntr
        usr.node[usr.name .== to] .= node_cntr
        if !isempty(cmp.edge[.!isnothing.(cmp.edge) .&& cmp.name .== to])
            idx = findfirst(.!isnothing.(cmp.edge) .&& cmp.name .== to)
            cmp.edge[idx] = (node_cntr, cmp.edge[idx][2])
        end

        # a nodal-cmp can be associated with multiple connections
        if !isempty(cmp.node[.!isnothing.(cmp.node) .&& cmp.name .== to])
            cmp_node = cmp.node[.!isnothing.(cmp.node) .&& cmp.name .== to][1]
            if iszero(cmp_node)
            cmp.node[.!isnothing.(cmp.node) .&& cmp.name .== to] .= node_cntr
            else
            idx = findfirst(cmp.name .== name)
            cmp.edge[idx] = (cmp.edge[idx][1], cmp_node)
            stop_cntr = true
            end
        end

        if !stop_cntr
            idx = findfirst(cmp.name .== name)
            cmp.edge[idx] = (cmp.edge[idx][1], node_cntr)
        end

        node_cntr += stop_cntr ? 0 : 1
        stop_cntr = false
    end
end

"""
    solve!(cmp, src)

Solve all State Transition Diagrams (STDs) for components and sources.

Iterates through all STDs in the component and source structures and solves any
that haven't been solved yet. Uses SteadyStateProcess as the default solver.

# Arguments
- `cmp`: Component data structure containing STDs
- `src`: Source data structure containing STDs

# Side Effects
- Modifies STDs by solving them and adding solution information
- Adds timing and probability information to STD properties

# Notes
- Skips STDs that already have time information (already solved)
- Uses SteadyStateProcess solver for reliability calculations
- Could be extended to use different solvers based on STD characteristics
"""
function solve!(cmp, src)
    # Check if the stds are already solved
    for std in Iterators.flatten((cmp.std, src.std))
        if haskey(std.props, :time)
        # elseif typeof(std.tprops[Graphs.Edge(1, 2)][:distr]) == typeof(Exponential(0.0u"yr")) && typeof(std.tprops[Graphs.Edge(2, 1)][:distr]) == typeof(Exponential(0.0u"yr"))
        #     solve!(std, SteadyStateProcess());
        #     println("Steady state process solved for $std")
        else
            solve!(std, SteadyStateProcess());
            # println("Semi-Markov process solved for $std")
        end
    end
end

## build
"""
    build_network(cmp, src, usr)

Construct a MultiStateSystems Network from categorized DCIDE elements.

Creates a Network object and populates it with sources, users, and components
in the proper categories based on their topology and characteristics.

# Arguments
- `cmp`: Component data structure with network topology information
- `src`: Source data structure 
- `usr`: User data structure

# Returns
- `Network`: Complete network model ready for reliability analysis

# Network Construction Process
1. Add all sources to the network
2. Add users and panel components (treated as special users)
3. Add nodal components (buses, switches at nodes)
4. Add unidirectional edge components (transformers, diodes)
5. Add bidirectional edge components (cables, bidirectional switches)

# Examples
```julia
# After initializing and connecting elements
cmp, src, usr = init_elements(json_data)
connect_elements!(cmp, src, usr, json_data)
solve!(cmp, src)

# Build the network
network = build_network(cmp, src, usr)
solve!(network)  # Solve network reliability
```

# Notes
- Panel components are treated as users in the network model
- Different component types are added using appropriate Network functions
- Network topology is preserved through node and edge assignments
"""
function build_network(cmp, src, usr)
    ntw = Network()
    
    add_sources!(ntw; src...)

    add_users!(ntw; usr...)
    add_users!(ntw; filter_panel_cmp(cmp)...)

    println(filter_bidirectional_edge_cmp(cmp))

    add_components!(ntw; filter_nodal_cmp(cmp)...)
    add_components!(ntw; filter_unidirectional_edge_cmp(cmp)...)
    add_bidirectional_components!(ntw; filter_bidirectional_edge_cmp(cmp)...)

    return ntw
end

## extend
"""
    extend_json!(json::Dict{String,Any}, ntw::AbstractNetwork)

Extend the original DCIDE JSON with computed reliability results.

Takes the solved network and adds the computed STDs and UGFs back to the original
JSON structure, enabling export of reliability results in DCIDE format.

# Arguments
- `json`: Original DCIDE JSON dictionary (modified in-place)
- `ntw::AbstractNetwork`: Solved network containing reliability results

# Side Effects
- Adds "std" fields to component and connection entries in JSON
- Adds "ugf" fields to user entries in JSON
- Preserves original JSON structure while augmenting with results

# Examples
```julia
# After solving network
network = build_network(cmp, src, usr)
solve!(network)

# Extend JSON with results
extend_json!(json_data, network)

# Now json_data contains reliability results
# Can be exported back to DCIDE format
```

# Notes
- Cable components are stored in connections array
- Other components are stored in componentInstances array  
- Users get UGF (Universal Generating Function) results
- Components and sources get STD (State Transition Diagram) results
"""
function extend_json!(json::Dict{String,Any}, ntw::AbstractNetwork)
    for cmp in ntw.cmp
        if cmp[:type] == "cable"
            json["connections"][cmp[:id]]["std"] = cmp[:std]
        else
            json["componentInstances"][cmp[:id]]["std"] = cmp[:std]
    end end
    for src in ntw.src
        json["componentInstances"][src[:id]]["std"] = src[:std]
    end
    for usr in ntw.usr 
        json["componentInstances"][usr[:id]]["ugf"] = usr[:ugf]
end end

## solve
"""
    solve!(json::Dict{String,Any})

Complete DCIDE integration workflow: parse, solve, and extend JSON with reliability results.

This is the main entry point for DCIDE integration. It performs the complete workflow:
1. Parse DCIDE JSON to extract system components  
2. Establish network connectivity
3. Solve component reliability models
4. Build and solve system-level network model
5. Export results back to JSON format

# Arguments
- `json::Dict{String,Any}`: DCIDE JSON specification (modified in-place)

# Side Effects
- Modifies the input JSON by adding reliability analysis results
- Computes steady-state availability and reliability indices
- Adds STD and UGF results to component specifications

# Examples
```julia
using JSON
using MultiStateSystems

# Load DCIDE system specification
json_data = JSON.parsefile("dcide_system.json")

# Perform complete reliability analysis
solve!(json_data)

# Results are now available in json_data
# Export enhanced JSON with results
JSON.print(open("results.json", "w"), json_data, 2)
```

# Workflow Steps
1. **Initialize Elements**: Parse JSON and categorize components
2. **Connect Elements**: Establish network topology 
3. **Solve Components**: Calculate individual component reliabilities
4. **Build Network**: Create system-level network model
5. **Solve Network**: Compute system reliability indices
6. **Extend JSON**: Add results back to original specification

# Notes
- Input JSON is modified in-place with results
- Supports all standard DCIDE component types
- Handles both DC and AC electrical systems
- Results include availability, reliability, and performance indices
"""
function solve!(json::Dict{String,Any})
    # extract cmp, src and usr from json
    cmp, src, usr = init_elements(json);
    connect_elements!(cmp, src, usr, json)
    # TODO - Extract the final value from the solved state transition diagram in case it is solve by a timeseries process.
    solve!(cmp,src)

    # build and solve ntw
    ntw = build_network(cmp, src, usr);
    solve!(ntw)

    # extend the json with availability results
    extend_json!(json, ntw)
end