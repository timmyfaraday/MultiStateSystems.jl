################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
System configuration parameters for LVDC short-circuit analysis examples.

This module defines all the system parameters, cable configurations, and protection
device specifications used in the MultiStateSystems.jl LVDC examples.
"""

using Unitful

"""
    SystemParameters

Structure containing all system-level electrical parameters.

# Fields
- `V_DC::Int`: DC voltage level (V)
- `I_max::Int`: Maximum current (A) 
- `V_min::Float64`: Minimum voltage threshold (V)
- `I²t::Int`: Energy let-through rating (A²s)
- `L_p::Float64`: System parasitic inductance (H)
- `C_b::Float64`: DC bus capacitance (F)
- `λᶜ`: Cable failure rate (1/yr/m)
- `P_zone::Vector{Float64}`: Zone protection probabilities
"""
struct SystemParameters
    V_DC::Int
    I_max::Int  
    V_min::Float64
    I²t::Int
    L_p::Float64
    C_b::Float64
    λᶜ
    P_zone::Vector{Float64}
end

"""
    IntegrationParameters

Structure containing numerical integration parameters.

# Fields
- `t_max`: Maximum simulation time
- `n::Int`: Number of integration points
- `tsim`: Total simulation time for STD analysis
- `dt`: Time step for STD analysis
"""
struct IntegrationParameters
    t_max
    n::Int
    tsim
    dt
end

"""
    ProtectionDevice

Structure representing a protection device with its characteristics.

# Fields
- `λ::Float64`: Failure rate parameter
- `type::String`: Device type ("CB" for circuit breaker, "Fuse" for fuse)
- `μ::Float64`: Mean time parameter
"""
struct ProtectionDevice
    λ::Float64
    type::String
    μ::Float64
end

"""
    CableConfiguration

Structure defining cable length configurations for different feeders.

# Fields
- `load_feeders::Dict`: Cable lengths for load feeders
- `source_feeders::Dict`: Cable lengths for source feeders  
- `bridge_feeders::Dict`: Cable lengths for bridge connections
"""
struct CableConfiguration
    load_feeders::Dict
    source_feeders::Dict
    bridge_feeders::Dict
end

"""
    get_default_system_parameters() -> SystemParameters

Returns the default system parameters for the LVDC analysis.
"""
function get_default_system_parameters()
    return SystemParameters(
        750,                    # V_DC
        200,                    # I_max
        750 * 0.85,            # V_min
        5000,                   # I²t
        0.0,                    # L_p
        2.5e-2,                # C_b (double bus capacitance)
        0.0000743u"1/yr/m",    # λᶜ
        [0.999, 1.0]           # P_zone
    )
end

"""
    get_default_integration_parameters() -> IntegrationParameters

Returns the default integration parameters for the analysis.
"""
function get_default_integration_parameters()
    return IntegrationParameters(
        0.1u"s",               # t_max
        100000,                # n
        25.0u"yr",             # tsim
        0.5u"d"                # dt
    )
end

"""
    get_protection_devices() -> Dict{String, ProtectionDevice}

Returns a dictionary of protection devices with their characteristics.
"""
function get_protection_devices()
    return Dict(
        "SSCB" => ProtectionDevice(0.05, "CB", 11e-6),
        "MCCB" => ProtectionDevice(0.05, "CB", 5.1e-3),
        "Fuse" => ProtectionDevice(0.05, "Fuse", 0.0)
    )
end

"""
    get_cable_configuration() -> CableConfiguration

Returns the default cable configuration with length ranges for different feeders.
"""
function get_cable_configuration()
    load_feeders = Dict(
        "C1" => 1u"m":1u"m":100u"m",
        "C2" => 1u"m":1u"m":200u"m", 
        "C3" => 1u"m":1u"m":150u"m",
        "C4" => 1u"m":1u"m":50u"m",
        "C5" => 1u"m":1u"m":100u"m"
    )
    
    source_feeders = Dict(
        "S1" => 1u"m":1u"m":20u"m",
        "S2" => 1u"m":1u"m":20u"m"
    )
    
    bridge_feeders = Dict(
        "C" => 1u"m":1u"m":200u"m"
    )
    
    return CableConfiguration(load_feeders, source_feeders, bridge_feeders)
end

"""
    get_device_cable_mapping() -> Dict{String, Dict}

Returns the mapping between protection devices and cable configurations.
"""
function get_device_cable_mapping()
    cable_config = get_cable_configuration()
    
    return Dict(
        "load" => Dict(
            "SSCB" => cable_config.load_feeders,
            "MCCB" => cable_config.load_feeders,
            "Fuse" => cable_config.load_feeders,
            "Fuse_MCCB" => cable_config.load_feeders
        ),
        "source" => Dict(
            "SSCB" => cable_config.source_feeders
        ),
        "bridge" => Dict(
            "SSCB" => cable_config.bridge_feeders
        )
    )
end

"""
    get_network_zones() -> Dict{String, Vector{String}}

Returns the definition of network zones for split-bus configurations.
"""
function get_network_zones()
    return Dict(
        "Zone_1" => ["C1", "C2"],
        "Zone_2" => ["C3", "C4", "C5"]
    )
end

"""
    get_ac_grid_availability() -> Float64

Returns the AC grid availability parameter.
"""
function get_ac_grid_availability()
    return 0.99999
end
