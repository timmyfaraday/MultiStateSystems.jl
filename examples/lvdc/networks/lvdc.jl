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

# solve the state-transitions diagrams
solve!(stdᵃᵈ, cls)
solve!(stdᵃᶜ, cls)
solve!(stdᵈᵈ, cls)
solve!(stdᵇ, cls)

# initialize the network corresponding to the low voltage dc system 
ntwˡᵛ = Network()

# add the users to the ntw
add_users!(ntwˡᵛ, node = [3,4])
add_sources!(ntwˡᵛ, node = [1,3], 
                    name = ["ac grid", "battery"], 
                    std  = [stdᵃᶜ,stdᵇ])
add_components!(ntwˡᵛ, edge = [(1,2),(2,3),(2,4)],
                       name = ["ac/dc", "dc/dc", "dc/ac"],
                       std  = [stdᵃᵈ, stdᵈᵈ, stdᵃᵈ])

# solve the network
solve!(ntwˡᵛ)