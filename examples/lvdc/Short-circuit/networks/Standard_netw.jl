################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
Run io.jl first, to load the solved state transition diagrams (STD) from the file std_s_data.dat.
Then run this script to calculate the results.
"""

# Set up the state transition diagram for the AC-grid with a constant availability of four nines.
ac_av = 0.9999*ones(length(time)); 
stdᵃᶜ = solvedSTD(prob = [ac_av, 1 .- ac_av], time = collect(time), power = [50.0u"MW", 0.0u"MW"]);

bat_av  = 0.999999*ones(length(time)); # Battery availability
std_s["Source"]["Bat"] = solvedSTD(prob = [bat_av, 1 .- bat_av], time = collect(time), power = [50.0u"MW", 0.0u"MW"]);

# Set up a perfectly available state transition diagram as source in some of the networks.
stdᵖ = solvedSTD(prob = [ones(length(time))], time = collect(time), power = [50.0u"MW"]);

# Solve the availability of the AC side of the network, including AC grid and AC/DC converter.
ntwᵃᶜ = Network()
add_user!(ntwᵃᶜ, node = 3)
add_sources!(ntwᵃᶜ, node = 1, std = stdᵃᶜ)
add_sources!(ntwᵃᶜ, node = 2, std = stdᵃᶜ)
add_components!(ntwᵃᶜ, edge = [(1,3), (2,3)],
                std = [std_s["Source"]["ACDC"], std_s["Source"]["ACDC"]]);
_MSS.solve!(ntwᵃᶜ)
ntw_ac = ntwᵃᶜ.usr[1][:ugf].prb[2] .+ ntwᵃᶜ.usr[1][:ugf].prb[3];

# Determining the probability that the AC side network availability has been restored before
# the battery has drained.
T = 5u"hr"
bat_av = Float64[]
append!(bat_av, ones(Int64(round(T/dt))))
append!(bat_av,[battery_system_availability(i, T, ntw_ac, log(4)u"d", 0.3u"d") for i in 1:1:length(time)-Int64(round(T/dt))])
stdᵇᵃᵛ = solvedSTD(prob = [bat_av, 1 .- bat_av], time = collect(time), power = [50.0u"MW", 0.0u"MW"])

# Determining the availability of the DC bus as a result of the untimely clearing of the circuit breakers

std_bus = calculate_bus_availability()
std_bus_C1 = calculate_bus_availability("C1")
std_bus_C2 = calculate_bus_availability("C2")
std_bus_C3 = calculate_bus_availability("C3")
std_bus_C4 = calculate_bus_availability("C4")
std_bus_C5 = calculate_bus_availability("C5")
    

#Total system availability
Total_bus_av = Dict()
for (key, L_c) in L_tot
    Total_bus_av[key] = Network()
    add_users!(Total_bus_av[key], node = 7)
    add_sources!(Total_bus_av[key], node = 1, std = stdᵖ)
    add_components!(Total_bus_av[key], edge = [(1,2),(2,5),(1,3),(3,5),(1,4),(4,5),(5,6),(6,7)],
                    std = [stdᵃᶜ, std_s["Source"]["ACDC"],stdᵃᶜ, std_s["Source"]["ACDC"], std_s["Source"]["Bat"], std_s["Source"]["DCDC"], stdᵇᵃᵛ, std_bus[key]])
    _MSS.solve!(Total_bus_av[key])
end

# Availability of the loads connected to the DC bus:
# C1
C1_av = Dict()
for (key, L_c) in L_tot
    C1_av[key] = Network()
    add_users!(C1_av[key], node = 8)
    add_sources!(C1_av[key], node = 1, std = stdᵖ)
    add_components!(C1_av[key], edge = [(1,2),(2,5),(1,3),(3,5),(1,4),(4,5),(5,6),(6,7),(7,8)],
                    std = [stdᵃᶜ, std_s["Source"]["ACDC"],stdᵃᶜ, std_s["Source"]["ACDC"], std_s["Source"]["Bat"], std_s["Source"]["DCDC"], stdᵇᵃᵛ, std_bus_C1[key], std_s[key]["C1"]])
    _MSS.solve!(C1_av[key])
end
C1 = C1_av["HCB"].usr[1][:ugf].prb[3]

# C2
C2_av = Dict() 
for (key, L_c) in L_tot
    C2_av[key] = Network()
    add_users!(C2_av[key], node = 8)
    add_sources!(C2_av[key], node = 1, std = stdᵖ)
    add_components!(C2_av[key], edge = [(1,2),(2,5),(1,3),(3,5),(1,4),(4,5),(5,6),(6,7),(7,8)],
                    std = [stdᵃᶜ, std_s["Source"]["ACDC"],stdᵃᶜ, std_s["Source"]["ACDC"], std_s["Source"]["Bat"], std_s["Source"]["DCDC"], stdᵇᵃᵛ, std_bus_C2[key], std_s[key]["C2"]])
    _MSS.solve!(C2_av[key])
end
C2 = C2_av["HCB"].usr[1][:ugf].prb[3]

# C3
C3_av = Dict()
for (key, L_c) in L_tot
    C3_av[key] = Network()
    add_users!(C3_av[key], node = 8)
    add_sources!(C3_av[key], node = 1, std = stdᵖ)
    add_components!(C3_av[key], edge = [(1,2),(2,5),(1,3),(3,5),(1,4),(4,5),(5,6),(6,7),(7,8)],
                    std = [stdᵃᶜ, std_s["Source"]["ACDC"],stdᵃᶜ, std_s["Source"]["ACDC"], std_s["Source"]["Bat"], std_s["Source"]["DCDC"], stdᵇᵃᵛ, std_bus_C3[key], std_s[key]["C3"]])
    _MSS.solve!(C3_av[key])
end
C3 = C3_av["HCB"].usr[1][:ugf].prb[3]

# C4
C4_av = Dict()
for (key, L_c) in L_tot
    C4_av[key] = Network()
    add_users!(C4_av[key], node = 8)
    add_sources!(C4_av[key], node = 1, std = stdᵖ)
    add_components!(C4_av[key], edge = [(1,2),(2,5),(1,3),(3,5),(1,4),(4,5),(5,6),(6,7),(7,8)],
                    std = [stdᵃᶜ, std_s["Source"]["ACDC"],stdᵃᶜ, std_s["Source"]["ACDC"], std_s["Source"]["Bat"], std_s["Source"]["DCDC"], stdᵇᵃᵛ, std_bus_C4[key], std_s[key]["C4"]])
    _MSS.solve!(C4_av[key])
end
C4 = C4_av["HCB"].usr[1][:ugf].prb[3]

# C5
C5_av = Dict()
for (key, L_c) in L_tot
    C5_av[key] = Network()
    add_users!(C5_av[key], node = 8)
    add_sources!(C5_av[key], node = 1, std = stdᵖ)
    add_components!(C5_av[key], edge = [(1,2),(2,5),(1,3),(3,5),(1,4),(4,5),(5,6),(6,7),(7,8)],
                    std = [stdᵃᶜ, std_s["Source"]["ACDC"],stdᵃᶜ, std_s["Source"]["ACDC"], std_s["Source"]["Bat"], std_s["Source"]["DCDC"], stdᵇᵃᵛ, std_bus_C5[key], std_s[key]["C5"]])
    _MSS.solve!(C5_av[key])
end

C5 = C5_av["HCB"].usr[1][:ugf].prb[3]