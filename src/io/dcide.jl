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

## std
# get_std(component) = solvedSTD(prob = [1], power = [(Inf)u"MW"]) 
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
""
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
""
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

function solve!(cmp, src)
    # Check if the stds are already solved
    for std in Iterators.flatten((cmp.std, src.std))
        if get_info(std, :solved)
            # If the std is already solved, skip it
            continue       
        else
            solve!(std, SteadyStateProcess());
        end
    end
end

## build
""
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
""
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
""
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