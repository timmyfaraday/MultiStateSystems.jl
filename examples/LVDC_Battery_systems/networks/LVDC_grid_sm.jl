################################################################################
#  Copyright 2022, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
This is an example to demonstrate the effect of wear-out failure of power electronic converters 
on the system wide availability in an LVDC distribution.

Pluto.jl environment is also available on at http://github.com/timmyfaraday/MultiStateSystems.jl
"""

# Load pkgs
using Plots
using Unitful
using MultiStateSystems

# Pkg const
const _MSS = MultiStateSystems

# Settings for a specific analysis
cls = SemiMarkovProcess()               
tsim = 40.0u"yr";                       # Total simulation time
dt = 8.0u"hr";                          # Simulation delta time
t   = zero(dt):dt:tsim .|> u"yr"        # Creating a timeseries of tsim and dt.
T = 2u"d";                              #Battery reserve time

# Function to determine probability of repair before battery reserve time is reached
function battery_system_availability(i, T, ntw_av, μ, σ)
        # Calculate the probability that the battery does not run out of charge before
        # the power sources are repaired at time t. With a lognormal distribution for 
        # the repair time with mean μ and standard deviation σ.
        p_failure = 1 - ntw_av[i]
        p_repair_time_ge_T = _MSS.ccdf(LogNormal(μ, σ),T)
        return 1-p_failure*p_repair_time_ge_T
end

# Include the stochastic models for the wind turbines and cables
stdᵃᶜᵈᶜ         = include(joinpath("C:/Users/gemmers/Documents/GitHub/MultiStateSystems.jl/examples/LVDC_Battery_systems/elements/acdc.jl"));
stdᵈᶜᵈᶜ         = include(joinpath("C:/Users/gemmers/Documents/GitHub/MultiStateSystems.jl/examples/LVDC_Battery_systems/elements/dcdc.jl"));
stdⁱⁿᵛ¹         = include(joinpath("C:/Users/gemmers/Documents/GitHub/MultiStateSystems.jl/examples/LVDC_Battery_systems/elements/inv1.jl"));
stdⁱⁿᵛ²         = include(joinpath("C:/Users/gemmers/Documents/GitHub/MultiStateSystems.jl/examples/LVDC_Battery_systems/elements/inv2.jl"));
stdᶜᵃᵖᵃᶜ        = include(joinpath("C:/Users/gemmers/Documents/GitHub/MultiStateSystems.jl/examples/LVDC_Battery_systems/elements/CAP_ac.jl"));
stdᶜᵃᵖᵈᶜ        = include(joinpath("C:/Users/gemmers/Documents/GitHub/MultiStateSystems.jl/examples/LVDC_Battery_systems/elements/CAP_dc.jl"));
stdᶜᵃᵖᴸ¹        = include(joinpath("C:/Users/gemmers/Documents/GitHub/MultiStateSystems.jl/examples/LVDC_Battery_systems/elements/CAP_L1.jl"));
stdᶜᵃᵖᴸ²        = include(joinpath("C:/Users/gemmers/Documents/GitHub/MultiStateSystems.jl/examples/LVDC_Battery_systems/elements/CAP_L2.jl"));
stdᵇᵃᵗ          = include(joinpath("C:/Users/gemmers/Documents/GitHub/MultiStateSystems.jl/examples/LVDC_Battery_systems/elements/battery.jl"));

# Solve the state transition diagrams
solve!(stdᵃᶜᵈᶜ, cls; tsim, dt)
solve!(stdᵈᶜᵈᶜ, cls; tsim, dt)
solve!(stdⁱⁿᵛ¹, cls; tsim, dt)
solve!(stdⁱⁿᵛ², cls; tsim, dt)
solve!(stdᶜᵃᵖᵃᶜ, cls; tsim, dt)
solve!(stdᶜᵃᵖᵈᶜ, cls; tsim, dt)
solve!(stdᶜᵃᵖᴸ¹, cls; tsim, dt)
solve!(stdᶜᵃᵖᴸ², cls; tsim, dt)
solve!(stdᵇᵃᵗ, cls; tsim, dt)

"""
The DC/DC converter can take a lot of memory to solve.
The solved std is provided as data and can be accessed by uncommenting the following two lines:
"""
# dcdc_data = CSV.read("Z:/MyOwn/Journal/semi-Markov/std_DCDC_40.csv", DataFrame)
# stdᵈᶜᵈᶜ = solvedSTD(prob = [dcdc_data[!, :dcdc_state1], dcdc_data[!, :dcdc_state2], dcdc_data[!, :dcdc_state3]], time = collect(t), power = [10.0u"MW", 0.0u"MW", 0.0u"MW"]);

# Set up the state transition diagram for the AC-grid with a constant availability of four nines.
ac_av = 0.9999*ones(length(t)); 
stdᵃᶜ = solvedSTD(prob = [ac_av, 1 .- ac_av], time = collect(t), power = [10.0u"MW", 0.0u"MW"]);

# Set up a perfectly available state transition diagram as source in some of the networks.
stdᵖ = solvedSTD(prob = [ones(length(t))], time = collect(t), power = [10.0u"MW"]);

# Solve the availability of the AC side of the network, including AC grid and AC/DC converter.
ntwᵃᶜ = Network()
add_user!(ntwᵃᶜ, node = 2)
add_sources!(ntwᵃᶜ, node = 1, std = stdᵃᶜ)
add_components!(ntwᵃᶜ, edge = [(1,2)],
                std = [stdᵃᶜᵈᶜ]);
solve!(ntwᵃᶜ)
ntw_ac = ntwᵃᶜ.usr[1][:ugf].prb[2];

# Determining the probability that the AC side network availability has been restored before
# the battery has drained.
bat_av = Float64[]
append!(bat_av, ones(Int64(round(T/dt))))
append!(bat_av,[battery_system_availability(i, T, ntw_ac, log(4)u"d", 0.3u"d") for i in 1:1:length(t)-Int64(round(T/dt))])
stdᵇᵃᵛ = solvedSTD(prob = [bat_av, 1 .- bat_av], time = collect(t), power = [10.0u"MW", 0.0u"MW"])

#Total system availability, contains the availability at all essential users in the sytem.
ntwᴸⱽᴰᶜ = Network()
add_users!(ntwᴸⱽᴰᶜ, node = [4,5,8,9,10])
add_sources!(ntwᴸⱽᴰᶜ, node = 1, std = stdᵖ)
add_components!(ntwᴸⱽᴰᶜ, edge = [(1,2),(2,4),(1,3),(3,4),(4,5),(5,6),(5,6),(6,7),(7,8),(8,9),(8,10)],
                std = [stdᵃᶜ, stdᵃᶜᵈᶜ, stdᵇᵃᵗ, stdᵈᶜᵈᶜ, stdᵇᵃᵛ, stdᶜᵃᵖᴸ¹, stdᶜᵃᵖᴸ², stdᶜᵃᵖᵃᶜ, stdᶜᵃᵖᵈᶜ, stdⁱⁿᵛ¹, stdⁱⁿᵛ²])

solve!(ntwᴸⱽᴰᶜ)