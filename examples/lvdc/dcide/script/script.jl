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

# Create component dictionary to add reliability data to the json data
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

extend_json_av!(json_data, Components)

# solve the problem
solve!(json_data)

# Extract the availability data of the network
# users: motor, bus and panel
# motor M-2:
json_data["componentInstances"][8]["ugf"].val # Shows the values of the potential power output of the motor
json_data["componentInstances"][8]["ugf"].prb # Shows the probability of the potential power output of the motor

# Extract the availability of each component in the system
# Front-end rectifier U-2:
_MSS.get_sprop(json_data["componentInstances"][3]["std"], :name) # Shows the names of the states the converter can be in.
_MSS.get_sprop(json_data["componentInstances"][3]["std"], :power) # Shows the power output of the converter in each state.