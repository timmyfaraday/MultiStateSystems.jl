using MultiStateSystems
using JSON
using Unitful

const _MSS = MultiStateSystems

json_file = joinpath(_MSS.BASE_DIR,"examples/lvdc/DCIDE/system_files/test_2.json")
data = JSON.parsefile(json_file)

const dcide_src_types = ["utility", "battery", "PV"]
const dcide_usr_types = ["motor"]
const dcide_cmp_types = ["panel", "switch", "protection", "converter", "transformer", "bus"]
const dcide_bus_types = ["bus", "panel"]

function init_elements(json)
    # Initialize the named tuples to store the components
    src = (name = [], std = [], node = [], type = [])
    cmp = (name = [], std = [], edge = [], node = [], bidirectional = [], type = [])
    usr = (name = [], std = [], node = [], type = [])

    for component in json["componentInstances"]
        if component["componentType"] in dcide_src_types
            push!(src.name, component["id"])
            push!(src.std, get_std(component))
            push!(src.node, 0)
            push!(src.type, component["componentType"])

        elseif component["componentType"] in dcide_usr_types
            push!(usr.name, component["id"])
            push!(usr.std, get_std(component))
            push!(usr.node, 0)
            push!(usr.type, component["componentType"])

        elseif component["componentType"] in dcide_cmp_types            
            push!(cmp.name, component["id"])
            push!(cmp.std, get_std(component))
            push!(cmp.bidirectional, true) # Assuming all components are bidirectional for now
            push!(cmp.type, component["componentType"])

            if length(component["configuration"]["ports"]) == 2 
                push!(cmp.node, nothing)
                push!(cmp.edge, (0,0))
            else
                push!(cmp.node, 0)
                push!(cmp.edge, nothing)
            end
        end
    end

    for connection in json["connections"]
        push!(cmp.name, connection["id"])
        push!(cmp.std, get_std(connection))
        push!(cmp.bidirectional, true) # Assuming all components are bidirectional for now
        push!(cmp.node, nothing)
        push!(cmp.edge, (0, 0))
        push!(cmp.type, "cable")
    end

    return cmp, src, usr #named tuples
end

# "component-X0vbSybQ"

# CB4X
# xMYz

# function solve!(json) #function to solve the STD
#     cmp, src, usr = init_elements(json)
#     connect_elements!(cmp, src, usr, json)

#     for nx in 
#     solve!(nx, SteadyStateProcess())
#     end
# end

filter_nodal_cmp(cmp) =
    (;zip(  [key for key in keys(cmp) if key ≠ :edge],
            [cmp[key][.!isnothing.(cmp.node)] for key in keys(cmp) if key ≠ :edge])...)
filter_bidirectional_edge_cmp(cmp) =
    (;zip(  [key for key in keys(cmp) if key ≠ :node],
            [cmp[key][.!isnothing.(cmp.edge) .&& cmp.bidirectional] for key in keys(cmp) if key ≠ :node])...)
filter_unidirectional_edge_cmp(cmp) =
    (;zip(  [key for key in keys(cmp) if key ≠ :node],
            [cmp[key][.!isnothing.(cmp.edge) .&& .!cmp.bidirectional] for key in keys(cmp) if key ≠ :node])...)
filter_panel_cmp(cmp) =
    (;zip(  [key for key in keys(cmp) if key ≠ :edge],
            [cmp[key][.!isnothing.(cmp.node) .&& in.(cmp.type, Ref(dcide_bus_types))] for key in keys(cmp) if key ≠ :edge])...)

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

function get_std(component)
    solvedSTD(prob = [1], power = [(Inf)u"MW"])
end

function build_network(cmp, src, usr)
    # Create a network from the components, sources, and users
    netw = _MSS.Network()
    add_sources!(netw, node = src.node, name = src.name, std = src.std)
    add_users!(netw, node = usr.node, name = usr.name, std = usr.std)
    add_users!(netw, node = filter_panel_cmp(cmp).node, 
               name = filter_panel_cmp(cmp).name, 
               std = filter_panel_cmp(cmp).std)
    add_bidirectional_components!(netw, edge = filter_bidirectional_edge_cmp(cmp).edge, 
                                  name = filter_bidirectional_edge_cmp(cmp).name, 
                                  std = filter_bidirectional_edge_cmp(cmp).std)
    add_components!(netw, edge = filter_unidirectional_edge_cmp(cmp).edge,
                   name = filter_unidirectional_edge_cmp(cmp).name, 
                   std = filter_unidirectional_edge_cmp(cmp).std)
    return netw
end

comp, source, user = init_elements(data)
connect_elements!(comp, source, user, data)

nodal_cmp = filter_nodal_cmp(comp)
panels = filter_panel_cmp(comp)

network = build_network(comp, source, user)

solve!(network)