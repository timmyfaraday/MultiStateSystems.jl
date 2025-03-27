using JSON
using MultiStateSystems

const _MSS = MultiStateSystems
# Load the JSON data from the file
json_file = joinpath(_MSS.BASE_DIR,"examples/lvdc/DCIDE/system_files/two_feeder.json")
json_data = JSON.parsefile(json_file)

# Extract component IDs and types
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
function extract_connections(json_data::Dict{String, Any}, component_details::Dict{String, Dict{String, Any}})::Dict{String, Vector{Tuple{String, String}}}
    connections = Dict{String, Vector{Tuple{String, String}}}()
    panel_connections = Dict{String, String}()
    connection_counter = 1

    for connection in json_data["connections"]
        from_component = connection["from"]["componentInstanceId"]
        to_component = connection["to"]["componentInstanceId"]

        if component_details[from_component]["type"] == "panel" || component_details[to_component]["type"] == "panel"
            panel_id = component_details[from_component]["type"] == "panel" ? from_component : to_component

            if !haskey(panel_connections, panel_id)
                panel_connections[panel_id] = string(connection_counter)
                connections[panel_connections[panel_id]] = Vector{Tuple{String, String}}()
                connection_counter += 1
            end

            push!(connections[panel_connections[panel_id]], (from_component, to_component))
        else
            connection_key = string(connection_counter)
            connections[connection_key] = [(from_component, to_component)]
            connection_counter += 1
        end
    end

    return connections
end


component_ids_and_types = extract_component_ids_and_types(json_data)

component_details = extract_component_details(json_data)

connections = extract_connections(json_data, component_details)


# Print the extracted information
println("Component IDs and Types:")
println(component_ids_and_types)

println("Component Details (IDs, Types, and Voltage Types):")
println(component_details)

println("Connections:")
println(connections)
