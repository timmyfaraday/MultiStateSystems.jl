################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# load pkgs
using Unitful
using MultiStateSystems

# pkg const
const _MSS = MultiStateSystems

# setting for a specific analysis
cls = SteadyStateProcess()

# include the stochastic models for the elements
stdᵃᵈ = include(joinpath(_MSS.BASE_DIR,"examples/lvdc/elements/acdc.jl"))
stdᵃᶜ = include(joinpath(_MSS.BASE_DIR,"examples/lvdc/elements/ac.jl"))
stdᵈᵈ = include(joinpath(_MSS.BASE_DIR,"examples/lvdc/elements/dcdc.jl"))
stdᵇ  = include(joinpath(_MSS.BASE_DIR,"examples/lvdc/elements/battery.jl"))
stdˢʷ = include(joinpath(_MSS.BASE_DIR,"examples/lvdc/elements/switch.jl"))
stdᶜ = include(joinpath(_MSS.BASE_DIR,"examples/lvdc/elements/switch.jl"))
# solve the state-transitions diagrams
solve!(stdᵃᵈ, cls)
solve!(stdᵃᶜ, cls)
solve!(stdᵈᵈ, cls)
solve!(stdᵇ, cls)
solve!(stdˢʷ, cls)
solve!(stdᶜ,cls)

# initialize the network corresponding to the low voltage dc system 

ntwˡᵛ = Network()


add_users!(ntwˡᵛ, node = [5])
add_sources!(ntwˡᵛ, node = [1,2], 
                    name = ["ac grid","ac grid"], 
                    std  = [stdᵃᶜ, stdᵃᶜ])
add_components!(ntwˡᵛ, edge = [(1,3),(1,4),(2,3),(2,4),(3,5), (4,5)],
                       name = ["cable", "cable", "cable", "cable", "cable", "cable"],
                       std  = [stdᶜ, stdᶜ, stdᶜ, stdᶜ, stdᶜ, stdᶜ])
add_bidirectional_component((ntwˡᵛ, edge = (3,4), std = stdᶜ))


solve!(ntwˡᵛ)

ntwˡᵛ.usr[1][:ugf]
ntwˡᵛ.usr[2][:ugf]