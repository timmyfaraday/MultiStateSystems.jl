################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
    Configuration

Configuration module for protection device reliability analysis.
Contains all simulation parameters, network configuration, and data paths.
"""

using Unitful

"""
    SimulationConfig

Structure containing simulation parameters.
"""
Base.@kwdef struct SimulationConfig
    # Time parameters
    total_time::typeof(1.0u"yr") = 25.0u"yr"
    time_step::typeof(1.0u"d") = 0.5u"d"
    
    # Battery parameters
    battery_backup_time::typeof(1.0u"hr") = 10.0u"hr"
    battery_availability::Float64 = 0.999999
    
    # AC grid parameters
    ac_grid_availability::Float64 = 0.9999
    
    # Repair time parameters (lognormal distribution)
    repair_mean::typeof(1.0u"d") = log(10)u"d"
    repair_std::typeof(1.0u"d") = 0.3u"d"
    
    # Power ratings
    source_power::typeof(1.0u"MW") = 50.0u"MW"
    infinite_power::typeof(1.0u"MW") = Inf*u"MW"
    zero_power::typeof(1.0u"MW") = 0.0u"MW"
end

"""
    NetworkConfig

Structure containing network topology and feeder configuration.
"""
Base.@kwdef struct NetworkConfig
    # Feeder specifications
    feeders::Dict{String, UnitRange} = Dict(
        "C1" => 1u"m":1u"m":100u"m",
        "C2" => 1u"m":1u"m":200u"m", 
        "C3" => 1u"m":1u"m":150u"m",
        "C4" => 1u"m":1u"m":50u"m",
        "C5" => 1u"m":1u"m":100u"m"
    )
    
    # Protection device types
    protection_devices::Vector{String} = ["SSCB", "MCCB", "HCB", "Fuse", "Fuse_MCCB"]
    
    # Network topology
    ac_sources::Int = 2
    dc_bus_nodes::Int = 3
end

"""
    DataPaths

Structure containing all file paths for data input/output.
"""
Base.@kwdef struct DataPaths
    base_dir::String = dirname(dirname(@__FILE__))
    data_dir::String = joinpath(base_dir, "data")
    results_dir::String = joinpath(base_dir, "results")
    
    # Input files
    std_data_file::String = joinpath(data_dir, "std_s_data_data.dat")
    source_data_file::String = joinpath(data_dir, "source_data.dat")
    
    # Output files
    results_file::String = joinpath(results_dir, "analysis_results.dat")
    
    function DataPaths(args...)
        paths = new(args...)
        # Ensure directories exist
        for dir in [paths.data_dir, paths.results_dir]
            if !isdir(dir)
                mkpath(dir)
            end
        end
        return paths
    end
end

"""
    AnalysisConfig

Main configuration structure combining all configuration types.
"""
struct AnalysisConfig
    simulation::SimulationConfig
    network::NetworkConfig
    paths::DataPaths
    
    function AnalysisConfig(;
        simulation_config::SimulationConfig = SimulationConfig(),
        network_config::NetworkConfig = NetworkConfig(),
        data_paths::DataPaths = DataPaths()
    )
        new(simulation_config, network_config, data_paths)
    end
end

"""
    get_default_config() -> AnalysisConfig

Get default configuration for the protection device analysis.

# Returns
- `AnalysisConfig`: Default configuration object
"""
function get_default_config()
    return AnalysisConfig()
end

"""
    validate_config(config::AnalysisConfig) -> Bool

Validate configuration parameters.

# Arguments
- `config::AnalysisConfig`: Configuration to validate

# Returns
- `Bool`: True if valid, false otherwise
"""
function validate_config(config::AnalysisConfig)
    # Validate simulation parameters
    if config.simulation.total_time <= 0u"yr"
        @error "Total simulation time must be positive"
        return false
    end
    
    if config.simulation.time_step <= 0u"d"
        @error "Time step must be positive"
        return false
    end
    
    if !(0.0 <= config.simulation.ac_grid_availability <= 1.0)
        @error "AC grid availability must be between 0 and 1"
        return false
    end
    
    if !(0.0 <= config.simulation.battery_availability <= 1.0)
        @error "Battery availability must be between 0 and 1"
        return false
    end
    
    # Validate network configuration
    if isempty(config.network.feeders)
        @error "Network must have at least one feeder"
        return false
    end
    
    if isempty(config.network.protection_devices)
        @error "Network must have at least one protection device type"
        return false
    end
    
    # Validate paths
    if !isdir(config.paths.base_dir)
        @error "Base directory does not exist: $(config.paths.base_dir)"
        return false
    end
    
    @info "Configuration validation passed"
    return true
end

"""
    print_config(config::AnalysisConfig)

Print configuration summary.

# Arguments
- `config::AnalysisConfig`: Configuration to display
"""
function print_config(config::AnalysisConfig)
    println("=== Protection Device Analysis Configuration ===")
    println("Simulation Parameters:")
    println("  Total Time: $(config.simulation.total_time)")
    println("  Time Step: $(config.simulation.time_step)")
    println("  Battery Backup: $(config.simulation.battery_backup_time)")
    println("  AC Availability: $(config.simulation.ac_grid_availability)")
    println("  Battery Availability: $(config.simulation.battery_availability)")
    
    println("\nNetwork Configuration:")
    println("  Feeders: $(length(config.network.feeders))")
    for (name, range) in config.network.feeders
        println("    $name: $(first(range))-$(last(range))")
    end
    println("  Protection Devices: $(join(config.network.protection_devices, ", "))")
    
    println("\nData Paths:")
    println("  Base Directory: $(config.paths.base_dir)")
    println("  Data Directory: $(config.paths.data_dir)")
    println("  Results Directory: $(config.paths.results_dir)")
    println("===============================================")
end
