################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# Load Pkgs
using XLSX
using Unitful
using MultiStateSystems

# Read in the wind power output data, and determine the clustered outputs and
# corresponding rates.
path = joinpath(BASE_DIR,"examples/data/Anholt.xlsx")
wind_power = XLSX.readdata(path,"Sheet1!D2:D52215")
output, rate = cluster_wind_power(wind_power,number_of_clusters = 7)

# Initialize the state-transition diagrams corresponding to the output (wto) and
# reliability (wtr) of a wind turbine.
stdʷᵗᵒ = STD()
stdʷᵗʳ = STD()

# Add the states to the std's
add_states!(stdʷᵗᵒ, flow = (output)u"MW",
                    init = [0.0,0.0,0.0,0.0,0.0,0.0,1.0])
add_states!(stdʷᵗʳ, flow = [(Inf)u"MW",0.0u"MW",0.0u"MW",0.0u"MW",0.0u"MW",
                            0.0u"MW",0.0u"MW",0.0u"MW",0.0u"MW",0.0u"MW"],
                    init = [1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])

# Add the transitions to the std's
add_transitions!(stdʷᵗᵒ, rate = (rate)u"1/yr")
add_transitions!(stdʷᵗʳ, states = [(1,2),(2,1)],
                         rate = [0.0590u"1/yr",0.0132u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,3),(3,1)],
                         rate = [0.0070u"1/yr",0.1695u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,4),(4,1)],
                         rate = [0.0770u"1/yr",0.0158u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,5),(5,1)],
                         rate = [0.0420u"1/yr",0.0361u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,6),(6,1)],
                         rate = [0.0240u"1/yr",0.3704u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,7),(7,1)],
                         rate = [0.3380u"1/yr",0.0442u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,8),(8,1)],
                         rate = [0.4320u"1/yr",0.0752u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,9),(9,1)],
                         rate = [0.4370u"1/yr",0.0625u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,10),(10,1)],
                         rate = [0.5380u"1/yr",0.0515u"1/hr"])

# Solve the problems
solve!(stdʷᵗᵒ, 1.0u"yr")
solve!(stdʷᵗʳ, 1.0u"yr")

# Plots.plot(ustrip.(stdʷᵗᵒ.props[:time]),[stdʷᵗᵒ.sprops[ns][:prob] for ns in 1:7],
#                        labels = reshape([string(stdʷᵗᵒ.sprops[ns][:flow]) for ns in 1:7],1,7))
