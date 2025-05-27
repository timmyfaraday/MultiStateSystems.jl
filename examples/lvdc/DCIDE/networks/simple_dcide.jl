################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# This Julia script takes a .json file containing a low voltage DC system generated
# by the DCIDE tool and solves it using the MultiStateSystems.jl package.

# First, we load the necessary packages. Then, we read the JSON file containing the system data.
# After that, the necessary elements are loaded, structures are created and the STDs are solved.
# Finally, the network is extracted from the data and solved.

### Package importans 
using Unitful
using MultiStateSystems
using JSON

const _MSS = MultiStateSystems

### Load the functions to read the JSON file and extract the data
include(joinpath(_MSS.BASE_DIR,"examples/lvdc/DCIDE/io/JSON_extr.jl"))
json_file = joinpath(_MSS.BASE_DIR,"examples/lvdc/DCIDE/system_files/parallel_feeder.json")
json_data = JSON.parsefile(json_file)

### Extract the data from the JSON file and save it in the right dictionary
component_details = extract_component_details(json_data)
connections = extract_connections(json_data, component_details)
comp_connections = component_connections(connections)
component_connection_counts = count_component_connections(connections)

### Find the necessary components that will create the nodes of the network
sources, loads = find_sources_and_loads(connections, component_details)
unattached_connections = find_component_connections(connections, sources, loads, component_details)
panels = find_panels(component_details)

### Create the nodes and edges and store them in dictionaries
nodes = node_dictionary(sources, loads, panels, unattached_connections, comp_connections, json_data)
edges = edge_dictionary(nodes, component_connection_counts)

### Solve the STDs in the dictionaries
solve!(nodes, edges, SteadyStateProcess())

### Create the network and solve it

netw = Network()
for src in sources
    add_source!(netw, node = nodes["sources"][src]["node"] , std = nodes["sources"][src]["std"])
end
for pnl in panels 
    add_user!(netw, node = nodes["panels"][pnl]["node"])
end
for lds in loads 
    add_user!(netw, node = nodes["loads"][lds]["node"])
end

add_components!(netw, edge = [comp["edge"] for (key, comp) in edges], 
                      name = [component_details[key]["type"] for (key, comp) in edges], 
                      std = [comp["std"] for (key, comp) in edges])

_MSS.solve!(netw)