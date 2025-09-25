################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# using pkgs
using JSON
using MultiStateSystems
using Unitful

const _MSS = MultiStateSystems

# path to and read-in json
json_file = joinpath(_MSS.BASE_DIR,"examples/lvdc/dcide/data/test.json")
json_data = JSON.parsefile(json_file)

# const dcide_bus_types = ["bus", "panel"]
# const dcide_cmp_types = ["panel", "switch", "protection", "converter", "transformer", "bus"]
# const dcide_src_types = ["utility", "battery", "PV"]
# const dcide_usr_types = ["motor"]

# dcide_av_keys = [:name, :type, :λ, :μ, :distr, :A]

# av_data = (; zip(dcide_av_keys, [[] for key in dcide_av_keys])...)

# Create component dictionary
Components = Dict(    
    "GRID-1" => Dict(
        :type => "source",
        :failure => Dict(
            :λ => nothing,
            :k => nothing,
            :distr => nothing
        ),
        :repair => Dict(
            :μ => nothing,
            :σ => nothing,
            :distr => nothing
        ),
        :A => 0.9999
    ),
    "U-1" => Dict(
        :type => "converter",
        :failure => Dict(
            :λ => 20.0u"yr",
            :k => nothing,
            :distr => Exponential
        ),
        :repair => Dict(
            :μ => 5.0u"hr",
            :σ => nothing,
            :distr => Exponential
        ),
        :A => nothing
    ),
    "U-2" => Dict(
        :type => "converter",
        :failure => Dict(
            :λ => 21.0u"yr",
            :k => nothing,
            :distr => Exponential
        ),
        :repair => Dict(
            :μ => 5.0u"hr",
            :σ => nothing,
            :distr => Exponential
        ),
        :A => nothing
    ),
    "U-3" => Dict(
        :type => "converter",
        :failure => Dict(
            :λ => 22.0u"yr",
            :k => nothing,
            :distr => Exponential
        ),
        :repair => Dict(
            :μ => 5.0u"hr",
            :σ => nothing,
            :distr => Exponential
        ),
        :A => nothing
    ),
    "U-4" => Dict(
        :type => "converter",
        :failure => Dict(
            :λ => 23.0u"yr",
            :k => nothing,
            :distr => Exponential
        ),
        :repair => Dict(
            :μ => 5.0u"hr",
            :σ => nothing,
            :distr => Exponential
        ),
        :A => nothing
    ),
    "U-5" => Dict(
        :type => "converter",
        :failure => Dict(
            :λ => 19.0u"yr",
            :k => nothing,
            :distr => Exponential
        ),
        :repair => Dict(
            :μ => 5.0u"hr",
            :σ => nothing,
            :distr => Exponential
        ),
        :A => nothing
    ),
    "BUS-1" => Dict(
        :type => "bus",
        :failure => Dict(
            :λ => 50.0u"yr",
            :k => nothing,
            :distr => Exponential
        ),
        :repair => Dict(
            :μ => 10.0u"hr",
            :σ => nothing,
            :distr => Exponential
        ),
        :A => nothing
    ),
)



# # Read out Components and push values into av_data
# for (name, comp) in Components
#     push!(av_data[:name], name)
#     push!(av_data[:type], comp[:type])
#     push!(av_data[:λ], get(comp[:failure], :λ, nothing))
#     push!(av_data[:μ], get(comp[:repair], :μ, nothing))
#     push!(av_data[:distr], get(comp[:failure], :distr, nothing))
#     push!(av_data[:A], get(comp, :A, nothing))
# end


function extend_json_av!(json::Dict{String,Any}, data::Dict)
    for (name, cmp) in data
        if cmp[:type] == "cable"
        else
            for (id, json_cmp) in enumerate(json["componentInstances"])
                if name == json_cmp["designator"]
                    println(id)
                    json["componentInstances"][id]["specifications"]["availability"] = data[name]
                end
            end
        end
    end
end

# Test to check whether the panel and buscomponents are nodal components

# Test to check that all nodal components have a node assigned to them

# Test to check that all edge components have two nodes assigned to them

# Test to check if the outcome of a single line example network is correct by comparing to a hard coded example

# Test to check if the outcome of a 


using Graphs
using GraphPlot

function visualize_components(cmp, src, usr)
    # Collect all nodes from cmp, src, usr
    nodes = Set{Int64}()
    edges = []

    # Add nodes from cmp.node, src.node, usr.node
    for n in cmp.node
        if n !== nothing
            push!(nodes, n)
        end
    end
    for n in src.node
        push!(nodes, n)
    end
    for n in usr.node
        push!(nodes, n)
    end

    # Add edges from cmp.edge (assume each edge is a tuple or pair of node names)
    for e in cmp.edge
        if e !== nothing
            push!(nodes, e[1])
            push!(nodes, e[2])
            push!(edges, (e[1], e[2]))
        end
    end

    # Map node names to indices for Graphs.jl
    node_list = collect(nodes)
    node_idx = Dict(n => i for (i, n) in enumerate(node_list))

    g = SimpleGraph(length(node_list))
    for (n1, n2) in edges
        add_edge!(g, node_idx[n1], node_idx[n2])
    end

    # Assign colors: sources=red, users=blue, others=gray
    node_colors = fill(:gray, length(node_list))
    for (i, n) in enumerate(node_list)
        if n in src.node
            node_colors[i] = :red
        elseif n in usr.node
            node_colors[i] = :blue
        end
    end

    # Visualize the graph
    gplot(g, nodelabel=node_list, nodefillc=node_colors)
end

function get_std(component)
    haskey(component["specifications"], "availability") || return nothing
    ports = get(component["specifications"]["electrical"], "ports", [])
    power = if !isempty(ports[1])
        ac = get(get(get(ports[1], "AC", Dict()), "power", Dict()), "nom", nothing)
        dc = get(get(get(ports[1], "DC", Dict()), "power", Dict()), "nom", nothing)
        ac !== nothing ? [ac * u"W", 0.0u"W"] :
        dc !== nothing ? [dc * u"W", 0.0u"W"] :
        [Inf * u"W", 0.0u"W"]
    else
        [Inf * u"W", 0.0u"W"]
    end
    
    # Determine the distribution vector for transitions
    av = component["specifications"]["availability"]
    if av[:A] !== nothing
        return solvedSTD(prob = [av[:A], 1 - av[:A]], power = power)
    end
    distr = [ av[:failure][:k] == nothing ? 
              av[:failure][:distr](av[:failure][:λ]) :
              av[:failure][:distr](av[:failure][:λ], av[:failure][:k]),
              av[:repair][:σ] == nothing ?
              av[:repair][:distr](av[:repair][:μ]) : av[:repair][:distr](av[:repair][:μ], av[:repair][:σ])]

    std = STD()
    add_states!(std, name=["Available", "Unavailable"], power=power, init=[1.0, 0.0])
    add_transitions!(std, states=[(1, 2), (2, 1)],
                       distr=distr)
    return std
end

function solve!(cmp, src, usr)
    # Check if the stds are already solved
    for std in Iterators.flatten((cmp.std, src.std, usr.std))
        if haskey(std.props, :time)
        elseif typeof(std.tprops[Graphs.Edge(1, 2)][:distr]) == typeof(Exponential(0.0u"yr")) && typeof(std.tprops[Graphs.Edge(2, 1)][:distr]) == typeof(Exponential(0.0u"yr"))
            solve!(std, SteadyStateProcess())
        else
            solve!(std, SemiMarkovProcess())
        end
    end
end
