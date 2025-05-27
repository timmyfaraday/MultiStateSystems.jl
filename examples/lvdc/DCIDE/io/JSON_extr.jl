using JSON
using MultiStateSystems
using Unitful

const _MSS = MultiStateSystems

"""
    extract_component_ids_and_types(json_data::Dict{String, Any})::Dict{String, String}

Extracts component IDs and their corresponding types from a JSON-like dictionary.

# Arguments
- `json_data::Dict{String, Any}`: A dictionary containing JSON data with a key `"componentInstances"`, 
  which is expected to be an array of dictionaries. Each dictionary should have keys `"id"` and `"componentType"`.

# Returns
- `Dict{String, String}`: A dictionary where the keys are component IDs (`"id"`) and the values are their 
  corresponding component types (`"componentType"`).
"""
function extract_component_ids_and_types(json_data::Dict{String, Any})::Dict{String, String}
    return Dict(component["id"] => component["componentType"] for component in json_data["componentInstances"])
end

"""
    extract_component_details(json_data::Dict{String, Any})::Dict{String, Dict{String, Any}}

Extracts details about components from a JSON-like dictionary structure and organizes them into a dictionary.

# Arguments
- `json_data::Dict{String, Any}`: A dictionary containing component data, expected to have a key `"componentInstances"` 
  which is an array of component dictionaries. Each component dictionary may contain an `"id"`, `"componentType"`, 
  and optionally a `"configuration"` with `"ports"`.

# Returns
- `Dict{String, Dict{String, Any}}`: A dictionary where each key is a component ID, and the value is another dictionary 
  containing:
  - `"type"`: The type of the component (from `"componentType"`).
  - `"voltageTypes"`: A formatted string representing the voltage types of the component's ports, based on their 
    `powerFlowDirection` and `voltageType`. If no valid voltage types are found, it defaults to `"Unknown"`.

# Details
- The function iterates over the `"componentInstances"` array in the input dictionary.
- For each component, it extracts the `"id"` and `"componentType"`.
- If the component has a `"configuration"` with `"ports"`, it processes the ports to determine the voltage types 
  based on their `powerFlowDirection` (`"input"`, `"output"`, or `"bidirectional"`) and `voltageType`.
- The voltage types are combined into a formatted string:
  - If both `"input"` and `"output"` voltage types are present, they are concatenated with a `/`.
  - If only one type is present, it is used directly.
  - If no valid types are found, it defaults to `"Unknown"`.
- The resulting details for each component are stored in the output dictionary.
"""
function extract_component_details(json_data::Dict{String, Any})::Dict{String, Dict{String, Any}}
    component_details = Dict{String, Dict{String, Any}}()

    for component in json_data["componentInstances"]
        id = component["id"]
        component_type = component["componentType"]
        voltage_types = Dict("input" => "", "output" => "", "bidirectional" => "")

        if haskey(component, "configuration") && haskey(component["configuration"], "ports")
            for port in component["configuration"]["ports"]
                voltage_type = get(port, "voltageType", "Unknown")
                direction = port["powerFlowDirection"]
                if direction in ["input", "output", "bidirectional"]
                    voltage_types[direction] = isempty(voltage_types[direction]) ? voltage_type : "$(voltage_types[direction])/$voltage_type"
                end
            end
        end

        formatted_voltage_types = if !isempty(voltage_types["input"]) && !isempty(voltage_types["output"])
            "$(voltage_types["input"])/$(voltage_types["output"])"
        elseif !isempty(voltage_types["input"])
            voltage_types["input"]
        elseif !isempty(voltage_types["output"])
            voltage_types["output"]
        elseif !isempty(voltage_types["bidirectional"])
            voltage_types["bidirectional"]
        else
            "Unknown"
        end

        component_details[id] = Dict("type" => component_type, "voltageTypes" => formatted_voltage_types)
    end

    return component_details
end

"""
    extract_connections(json_data::Dict{String, Any}, component_details::Dict{String, Dict{String, Any}}) -> Dict{String, Vector{Tuple{String, String}}}

Extracts and organizes connection data from a JSON-like structure into a dictionary format.

# Arguments
- `json_data::Dict{String, Any}`: A dictionary containing connection data, where each connection specifies a `from` and `to` component instance ID.
- `component_details::Dict{String, Dict{String, Any}}`: A dictionary containing details about each component, including its type.

# Returns
- `Dict{String, Vector{Tuple{String, String}}}`: A dictionary where keys are unique connection identifiers (strings) and values are vectors of tuples representing connections between components. Connections involving panels are grouped under the same key.

# Details
- Connections involving components of type `"panel"` are grouped together under a shared key, which is determined by the panel's unique identifier.
- Non-panel connections are assigned unique keys based on a counter.
- The function ensures that all connections are organized into a structured format for further processing.
"""
function extract_connections(json_data::Dict{String, Any}, component_details::Dict{String, Dict{String, Any}})::Dict{Int64, Vector{Tuple{String, String}}}
    connections = Dict{Int64, Vector{Tuple{String, String}}}()
    panel_connections = Dict()
    connection_counter = 1

    for connection in json_data["connections"]
        from_component = connection["from"]["componentInstanceId"]
        to_component = connection["to"]["componentInstanceId"]

        if component_details[from_component]["type"] == "panel" || component_details[to_component]["type"] == "panel"
            panel_id = component_details[from_component]["type"] == "panel" ? from_component : to_component

            if !haskey(panel_connections, panel_id)
                panel_connections[panel_id] = connection_counter
                connections[panel_connections[panel_id]] = Vector{Tuple{String, String}}()
                connection_counter += 1
            end

            push!(connections[panel_connections[panel_id]], (from_component, to_component))
        else
            connection_key = connection_counter
            connections[connection_key] = [(from_component, to_component)]
            connection_counter += 1
        end
    end

    return connections
end

"""
    component_connections(connections::Dict{Int64, Vector{Tuple{String, String}}})::Dict{String, Vector{String}}

Returns a dictionary mapping each component ID to a vector of component IDs it is directly connected to.

# Arguments
- `connections::Dict{Int64, Vector{Tuple{String, String}}}`: Dictionary of connection keys to vectors of (from, to) tuples.

# Returns
- `Dict{String, Vector{String}}`: Dictionary where each key is a component ID and the value is a vector of connected component IDs.
"""
function component_connections(connections::Dict{Int64, Vector{Tuple{String, String}}})::Dict{String, Vector{String}}
    conn_dict = Dict{String, Set{String}}()
    for connlist in values(connections)
        for (from, to) in connlist
            if !haskey(conn_dict, from)
                conn_dict[from] = Set{String}()
            end
            if !haskey(conn_dict, to)
                conn_dict[to] = Set{String}()
            end
            push!(conn_dict[from], to)
            push!(conn_dict[to], from)
        end
    end
    # Convert sets to vectors for output
    return Dict(k => collect(v) for (k, v) in conn_dict)
end

"""
    count_component_connections(connections::Dict{String, Vector{Tuple{String, String}}})::Dict{String, Int}

Counts the number of connections for each unique component in the connections dictionary.

# Arguments
- `connections::Dict{String, Vector{Tuple{String, String}}}`: A dictionary where keys are connection identifiers and values are vectors of tuples representing connections between components.

# Returns
- `Dict{String, Int}`: A dictionary where keys are component IDs and values are the number of times each component is connected to other components.
"""
function count_component_connections(connections::Dict{Int64, Vector{Tuple{String, String}}})::Dict{String, Int}
    connection_counts = Dict{String, Int}()

    for connection_list in values(connections)
        for (from_component, to_component) in connection_list
            connection_counts[from_component] = get(connection_counts, from_component, 0) + 1
            connection_counts[to_component] = get(connection_counts, to_component, 0) + 1
        end
    end

    return connection_counts
end

## Create networks from the connections and component details

# Find sources and loads in the JSON file
function find_sources_and_loads(connections::Dict{Int64, Vector{Tuple{String, String}}}, component_details::Dict{String, Dict{String, Any}})
    sources = String[]
    loads = String[]
    connection_counts = count_component_connections(connections)

    for (component_id, count) in connection_counts
        if count == 1
            if component_details[component_id]["type"] in ["battery", "utility"]
                push!(sources, component_id)
            else
                push!(loads, component_id)
            end
        end
    end

    return sources, loads
end


"""
    find_component_connections(connections::Dict{Int64, Vector{Tuple{String, String}}}, sources::Vector{String}, loads::Vector{String}, component_details::Dict{String, Dict{String, Any}}) -> Dict{Int64, Vector{Tuple{String, String}}}

Finds all connections that are not connected to any source, load, or panel.

# Arguments
- `connections`: Dictionary of connection keys to vectors of (from, to) tuples.
- `sources`: Vector of source component IDs.
- `loads`: Vector of load component IDs.
- `component_details`: Dictionary of component details keyed by component ID.

# Returns
- `Dict{Int64, Vector{Tuple{String, String}}}`: Dictionary of connection keys and their tuples that are not attached to any source, load, or panel.
"""
function find_component_connections(connections::Dict{Int64, Vector{Tuple{String, String}}}, sources::Vector{String}, loads::Vector{String}, component_details::Dict{String, Dict{String, Any}})
    unattached = Dict{Int64, Vector{Tuple{String, String}}}()
    important_ids = Set(sources) ∪ Set(loads)
    for (cid, connlist) in connections
        keep = true
        for (from, to) in connlist
            if from in important_ids || to in important_ids ||
               component_details[from]["type"] == "panel" || component_details[to]["type"] == "panel"
                keep = false
                break
            end
        end
        if keep
            unattached[cid] = connlist
        end
    end
    return unattached
end


"""
    find_panels(component_details::Dict{String, Dict{String, Any}})::Vector{String}

Returns a vector of component IDs where the component type is "panel".

# Arguments
- `component_details::Dict{String, Dict{String, Any}}`: Dictionary of component details keyed by component ID.

# Returns
- `Vector{String}`: Vector of component IDs that are panels.
"""
function find_panels(component_details::Dict{String, Dict{String, Any}})::Union{Vector{String}, Nothing}
    return [id for (id, details) in component_details if details["type"] == "panel"]
end

"""
    collect_availability_data(id::Union{String, Int}, param::Symbol, must_have::Bool)

Collects the availability data (λ or μ) for a given component or connection.

# Arguments
- `id::Union{String, Int}`: The component or connection identifier.
- `param::Symbol`: Either `:λ` or `:μ` indicating which parameter to collect.
- `must_have::Bool`: If true, a value must be returned (fallback to another source if not found); if false, return `nothing` if not found.

# Returns
- The value of the requested parameter, or `nothing` if not found and `must_have` is false.
"""
function collect_availability_data(id, param::Symbol, must_have::Bool, json_data::Dict{String, Any} = json_data)
    # Try to find data in the JSON file
    # For components
    if isa(id, String)
        for comp in json_data["componentInstances"]
            if comp["id"] == id && haskey(comp, "availability")
                avail = comp["availability"]
                if param == :λ && haskey(avail, "lambda")
                    return avail["lambda"]
                elseif param == :μ && haskey(avail, "mu")
                    return avail["mu"]
                end
            end
        end
    # For connections (by Int key)
    elseif isa(id, Int)
        # If your JSON has connection availability, implement here
        # Example (if available):
        # for conn in json_data["connections"]
        #     if conn["id"] == id && haskey(conn, "availability")
        #         avail = conn["availability"]
        #         if param == :λ && haskey(avail, "lambda")
        #             return avail["lambda"]
        #         elseif param == :μ && haskey(avail, "mu")
        #             return avail["mu"]
        #         end
        #     end
        # end
        # Otherwise, skip for now
    end

    if must_have
        if param == :λ
            return 0.01u"yr"
        elseif param == :μ
            return 0.2u"hr"
        else
            error("Unknown parameter $(param) for availability data fallback.")
        end
    else
        return nothing
    end
end

function node_dictionary(sources, loads, panels, unattached_connections, component_connections, json_data)
    node_dict = Dict{String, Any}()
    node_count = 1
    node_dict["sources"] = Dict{String, Any}()
    for src in sources
        node_dict["sources"][src] = Dict(
            "node" => node_count,
            "type" => "source",
            "λ" => collect_availability_data(src, :λ, true, json_data),
            "μ" => collect_availability_data(src, :μ, true, json_data),
            "connected" => component_connections[src]
        )
        node_count += 1
    end
    node_dict["panels"] = Dict{String, Any}()
    for pnls in panels
        node_dict["panels"][pnls] = Dict(
            "node" => node_count,
            "type" => "panel",
            "λ" => collect_availability_data(pnls, :λ, false, json_data),
            "μ" => collect_availability_data(pnls, :μ, false, json_data),
            "connected" => component_connections[pnls]
        )
        node_count += 1
    end
    node_dict["connections"] = Dict{Any, Any}()
    for (key,conn) in unattached_connections
        node_dict["connections"][conn] = Dict(
            "node" => node_count,
            "type" => "unattached",
            "λ" => collect_availability_data(conn, :λ, false, json_data),
            "μ" => collect_availability_data(conn, :μ, false, json_data),
            "connected" => [c for c in conn[1]]
        )
        node_count += 1
    end
    node_dict["loads"] = Dict{String, Any}()
    for lds in loads
        node_dict["loads"][lds] = Dict(
            "node" => node_count,
            "type" => "load",
            "λ" => collect_availability_data(lds, :λ, false, json_data),
            "μ" => collect_availability_data(lds, :μ, false, json_data),
            "connected" => component_connections[lds]
        )
        node_count += 1
    end
    return node_dict    
end

function edge_dictionary(nodes, component_connection_counts)
    edge_dict = Dict{String, Dict{String, Any}}()
    node_keys = ("sources", "panels", "connections", "loads")

    for (component, connection_amount) in component_connection_counts
        if connection_amount == 2
            node_indices = Int[]
            for key in node_keys
                for node_info in values(nodes[key])
                    if component in node_info["connected"]
                        push!(node_indices, node_info["node"])
                        if length(node_indices) == 2
                            break
                        end
                    end
                end
                if length(node_indices) == 2
                    break
                end
            end
            if length(node_indices) == 2
                edge_dict[component] = Dict(
                    "edge" => (minimum([node_indices[1] node_indices[2]]), maximum([node_indices[1] node_indices[2]])),
                    "λ" => collect_availability_data(component, :λ, true),
                    "μ" => collect_availability_data(component, :μ, true)
                )
            end
        end
    end
    return edge_dict
end

# Solve the different state transition diagrams of the components found in nodes and edges, by collecting the data
"""
    has_lambda_or_mu(dict::Dict{String, Any})::Bool

Checks if the given dictionary contains a key "λ" or "μ" with a non-nothing value.

# Arguments
- `dict::Dict{String, Any}`: The dictionary to check.

# Returns
- `Bool`: `true` if "λ" or "μ" is present and not `nothing`, otherwise `false`.
"""
function has_lambda_or_mu(dict::Dict{String, Any})::Bool
    return (haskey(dict, "λ") && dict["λ"] !== nothing) && (haskey(dict, "μ") && dict["μ"] !== nothing)
end

function solve!(nodes::Dict{String, Any}, edges::Dict{String, Dict{String, Any}}, cls)
    # This function would contain the logic to solve the state transition diagrams
    # and add to the dictionary of the nodes or edges based on the nodes and edges provided.

    for type in keys(nodes)
        for component in keys(nodes[type])
            if has_lambda_or_mu(nodes[type][component])
                std = STD()
                # add the states to the std
                add_states!(std, name  = ["available", "unavailable"],
                                power = [(Inf)u"MW", 0.0u"MW"],
                                init  = [1.0, 0.0])

                # add the transitions to the std
                add_transitions!(std,   states = [(1,2), (2,1)],
                                        distr = [Exponential(nodes[type][component]["λ"]),
                                                 Exponential(nodes[type][component]["μ"])])

                _MSS.solve!(std, cls)
                nodes[type][component]["std"] = std                
            end
        end
    end
    for component in keys(edges)
        if has_lambda_or_mu(edges[component])
            std = STD()
                # add the states to the std
                add_states!(std, name  = ["available", "unavailable"],
                                power = [(Inf)u"MW", 0.0u"MW"],
                                init  = [1.0, 0.0])

                # add the transitions to the std
                add_transitions!(std,   states = [(1,2), (2,1)],
                                        distr = [Exponential(edges[component]["λ"]),
                                                 Exponential(edges[component]["μ"])])

                _MSS.solve!(std, cls)
                edges[component]["std"] = std   
        end
    end
end

# # Order components starting from sources and assign numbers
# function order_components(sources::Vector{String}, connections::Dict{Int64, Vector{Tuple{String, String}}}, component_ids_and_types::Dict{String, String})::Dict{Int, String}
#     ordered, visited, num = Dict{Int, String}(), Set{String}(), 1

#     function dfs(comp::String, stop_at_panel::Bool)
#         if comp in visited
#             return
#         end
#         push!(visited, comp)
#         ordered[num] = comp
#         num += 1
#         if stop_at_panel && component_ids_and_types[comp] == "panel"
#             return
#         end
#         for conn in values(connections)
#             for (from, to) in conn
#                 if comp in (from, to) && !(to == comp ? from : to in visited)
#                     dfs(to == comp ? from : to, stop_at_panel)
#                 end
#             end
#         end
#     end

#     foreach(src -> dfs(src, true), sources)
#     foreach(comp -> dfs(comp[1], false), filter(x -> x[2] == "panel" && !(x[1] in visited), component_ids_and_types))
#     foreach(comp -> dfs(comp[1], false), filter(x -> !(x[1] in visited), component_ids_and_types))
#     return ordered
# end



# # Order components starting from the sources and attribute numbers
# ordered_components_with_numbers = order_components(sources, connections, component_ids_and_types)

# # Print the ordered components with numbers
# println("Ordered Components with Numbers:")
# println(ordered_components_with_numbers)

function find_component(connections::Dict{Int64, Vector{Tuple{String, String}}}, component::String)
    matching_keys = Int64[]

    for (key, connection_list) in connections
        for (from_component, to_component) in connection_list
            if from_component == component || to_component == component
                push!(matching_keys, key)
                break
            end
        end
    end

    if length(matching_keys) == 1
        return matching_keys[1]
    elseif length(matching_keys) == 2
        return (matching_keys[1], matching_keys[2])
    else
        return nothing
    end
end

