################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

using MultiStateSystems
using Unitful
using Plots
using DataFrames

const _MSS = MultiStateSystems

include(joinpath(_MSS.BASE_DIR,"examples/lvdc/Short-circuit/functions/DC_Faults.jl"))

# Grid characteristics:
V_DC    = 750; # V
I_max   = 200; # A
<<<<<<< Updated upstream
V_min   = V_DC*0.5; # V
=======
V_min   = V_DC*0.85; # V
>>>>>>> Stashed changes
I²t     = 5000; #A²t
L_p     = 0.0; # H
C_b     = 5e-2;# F
λᶜ = 0.0000743u"1/yr/m"; # Cable failure rate

# Define the characteristics of the protection devices
λ = Dict("SSCB" => [0.05, "CB"], 
         "HCB"  => [0.05, "CB"],
         "MCCB" => [0.05, "CB"],
         "Fuse" => [0.05, "Fuse"])
         
μ = Dict("SSCB" => 11e-6, 
         "HCB"  => 1.1e-3,
         "MCCB" => 5.1e-3,
         "Fuse" => 0)

# Define the feeder lengths
L_c = Dict("C1" => 1u"m":1u"m":100u"m", 
           "C2" => 1u"m":1u"m":200u"m",
           "C3" => 1u"m":1u"m":150u"m",
           "C4" => 1u"m":1u"m":50u"m",
           "C5" => 1u"m":1u"m":100u"m")

# Define the protection device used for each feeder
L_tot = Dict("SSCB" => L_c,
             "MCCB" => L_c,
             "HCB"  => L_c,
<<<<<<< Updated upstream
             "Fuse" => L_c)
=======
             "Fuse" => L_c,
             "Fuse_MCCB" => L_c)
>>>>>>> Stashed changes

# Write function and add option to combine different types of protection devices.

# Integration parameters
t_max   = 0.1u"s"
n       = 100000

# Calculate the probability of clearing a fault for each fault location and protection device
P = Dict()
<<<<<<< Updated upstream
for (key1, L_c) in L_tot
    P_i = Dict()
    for (key, value) in L_c
        if λ[key1][2] == "CB"
            P_i[key] = fault_clear_prob(value, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[key1], λ[key1][1])
            println("Protection device is of CB type")
        else
            P_i[key] = fault_clear_prob_fuse(value, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[key1][1])
            println("Protection device is a fuse")
        end
    end
    P[key1] = P_i
end

Pc = Dict()
for (key1, L_c) in L_tot
    Pc_i = Dict()
    for (key, value) in L_c
        if λ[key1][2] == "CB"
            Pc_i[key] = fault_clear_prob(maximum(value), L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[key1], λ[key1][1])
            println("Protection device is of CB type")
        else
            Pc_i[key] = fault_clear_prob_fuse(maximum(value), L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[key1][1])
            println("Protection device is a fuse")
        end
    Pc[key1] = Pc_i
    end
=======
for (device_type, feeder_lengths) in L_tot
    probabilities = Dict()

    if occursin("_", device_type)
        device_a, device_b = split(device_type, "_")
        for (feeder, lengths) in feeder_lengths
            prob_a = λ[device_a][2] == "CB" ? fault_clear_prob(lengths, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_a], λ[device_a][1]; return_mean = false) :
                                              fault_clear_prob_fuse(lengths, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_a][1]; return_mean = false)
            prob_b = λ[device_b][2] == "CB" ? fault_clear_prob(lengths, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_b], λ[device_b][1]; return_mean = false) :
                                              fault_clear_prob_fuse(lengths, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_b][1]; return_mean = false)
            probabilities[feeder] = mean(max.(prob_a, prob_b))
        end
    else
        for (feeder, lengths) in feeder_lengths
            probabilities[feeder] = λ[device_type][2] == "CB" ? fault_clear_prob(lengths, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_type], λ[device_type][1]) :
                                                                fault_clear_prob_fuse(lengths, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_type][1])
        end
    end

    P[device_type] = probabilities
end

Pc = Dict()
for (device_type, feeder_lengths) in L_tot
    probabilities = Dict()

    if occursin("_", device_type)
        device_a, device_b = split(device_type, "_")
        for (feeder, lengths) in feeder_lengths
            prob_a = λ[device_a][2] == "CB" ? fault_clear_prob(maximum(lengths), L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_a], λ[device_a][1]; return_mean = false) :
                                              fault_clear_prob_fuse(maximum(lengths), L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_a][1]; return_mean = false)
            prob_b = λ[device_b][2] == "CB" ? fault_clear_prob(maximum(lengths), L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_b], λ[device_b][1]; return_mean = false) :
                                              fault_clear_prob_fuse(maximum(lengths), L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_b][1]; return_mean = false)
            probabilities[feeder] = mean(max.(prob_a, prob_b))
        end
    else
        for (feeder, lengths) in feeder_lengths
            probabilities[feeder] = λ[device_type][2] == "CB" ? fault_clear_prob(maximum(lengths), L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_type], λ[device_type][1]) :
                                                                fault_clear_prob_fuse(maximum(lengths), L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_type][1])
        end
    end

    Pc[device_type] = probabilities
>>>>>>> Stashed changes
end

std = Dict()
h = Dict()
tsim = 25.0u"yr";  #25
dt = 0.5u"d";
time = 0.0u"yr":dt:tsim .|>u"yr"
cls = SemiMarkovProcess()

function solve_CB(P, Pc, λᶜ,n)    
    stdᶜᵇ = STD()
    # add the states to the std
    add_states!(stdᶜᵇ, name  = ["A/A", "A/U_1", "A/U_2","A/U_1", "A/U_2","A/U_1", "A/U_2"],
        power = [(Inf)u"MW", 10.0u"MW", 5.0u"MW", 10.0u"MW", 5.0u"MW", 10.0u"MW", 5.0u"MW"],
        init  = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    add_transitions!(stdᶜᵇ, states = [(1,2),(1,3),(2,1),(3,1),(1,4),(1,5),(4,1),(5,1),(1,6),(1,7),(6,1),(7,1)],
        distr = [   Exponential(1/(λᶜ), P/n),
                    Exponential(1/(λᶜ), (1-P)/n),
                    LogNormal(log(14.0)u"d", 0.1u"d"),
                    LogNormal(log(14.0)u"d", 0.1u"d"),
                    Weibull(21.0u"yr", 3.51, P/n), # Cable
                    Weibull(21.0u"yr", 3.51, (1-P)/n), # Cable
                    LogNormal(log(14.0)u"d", 0.1u"d"),
                    LogNormal(log(14.0)u"d", 0.1u"d"),
                    Weibull(15.0u"yr", 3.81, Pc/n), # Converter
                    Weibull(15.0u"yr", 3.81, (1-Pc)/n), # Converter
                    LogNormal(log(10.0)u"d", 0.2u"d"),
                    LogNormal(log(10.0)u"d", 0.2u"d")])
    return stdᶜᵇ
end

#Setting up the state transition diagrams 
for (key, L_c) in L_tot
    std_i = Dict()
    for (key_a, value_a) in L_c
        std_i[key_a] = solve_CB(P[key][key_a], Pc[key][key_a], λᶜ*maximum(value_a), 3)
    end
    std[key] = std_i
end

# Solving the state transition diagrams and extracting the transition rates.
for (key, std_i) in std
    h_i = Dict()
    @time for (cb, std_sol) in std_i
        H = _MSS.solve!(std_sol, cls; tsim, dt)
        println("solved for ", cb)
        Nt = length(_MSS.get_prop(std_sol, :time))
        Ns = _MSS.ns(std_sol)
        h2 = Dict{Int, Vector}()
        for st in _MSS.states(std_sol)
            h2[st] = [H[Ns * (x-1) + st] for x in 1:Nt]
        end
        h_i[cb] = h2
    end
    h[key] = h_i
end

# Calculate the probability of being in state 3, which is derived from the transition probability density.
# This step is taken because the recovery action is significantly faster than other maintenance actions.
prb_state3 = Dict()
for (key, value) in h
    prb_state3_i = Dict()
    for(cb,h_col) in value
        prb_cb = Dict()
        for i in [3, 5, 7]
            @time prb_cb[i] = state_conv(LogNormal(log(2.0)u"hr", 0.25u"hr"), h_col[i], time, 10000)
        end
        prb_state3_i[cb] = prb_cb
    end
    prb_state3[key] = prb_state3_i
end

# Set up the new state transition diagram, also including state 3.
std_s = Dict()
for (key, value) in std
    std_s_i = Dict()
    for (cb, std_sol) in value
        std_s_i[cb] = solvedSTD(prob = [_MSS.get_sprop(std[key][cb], :prob)[1], 
                                    vcat(0,_MSS.get_sprop(std[key][cb], :prob)[2][2:end] .+ _MSS.get_sprop(std[key][cb], :prob)[4][2:end] .+ _MSS.get_sprop(std[key][cb], :prob)[6][2:end]),
                                    vcat(0,prb_state3[key][cb][3][2:end] .+ prb_state3[key][cb][5][2:end] .+ prb_state3[key][cb][7][2:end]),
                                    vcat(0,_MSS.get_sprop(std[key][cb], :prob)[3][2:end] .+ _MSS.get_sprop(std[key][cb], :prob)[5][2:end] .+ _MSS.get_sprop(std[key][cb], :prob)[7][2:end] .- prb_state3[key][cb][3][2:end] .- prb_state3[key][cb][5][2:end] .- prb_state3[key][cb][7][2:end])],
                            time = collect(time),
                            power = [(Inf)u"MW", 10.0u"MW", 0.0u"MW", 10.0u"MW"])
    end
    std_s[key] = std_s_i
end

