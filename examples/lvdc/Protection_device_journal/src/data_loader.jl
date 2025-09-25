################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
    DataLoader

Module for loading and processing pre-computed state transition diagram (STD) data
for protection device reliability analysis in LVDC networks.
"""

using MultiStateSystems
using Serialization
using Unitful

const _MSS = MultiStateSystems

"""
    load_std_data(file_path::String) -> Dict

Load serialized state transition diagram data from file.

# Arguments
- `file_path::String`: Path to the serialized data file

# Returns
- `Dict`: Dictionary containing the STD data

# Throws
- `ArgumentError`: If file does not exist
"""
function load_std_data(file_path::String)
    if !isfile(file_path)
        throw(ArgumentError("File does not exist: $file_path"))
    end
    
    try
        return deserialize(file_path)
    catch e
        throw(ArgumentError("Failed to deserialize file $file_path: $e"))
    end
end

"""
    create_time_vector(; tsim=25.0u"yr", dt=0.5u"d") -> Vector

Create a time vector for simulation.

# Keyword Arguments
- `tsim`: Total simulation time (default: 25.0 years)
- `dt`: Time step (default: 0.5 days)

# Returns
- `Vector`: Time vector in years
"""
function create_time_vector(; tsim=25.0u"yr", dt=0.5u"d")
    return collect(0.0u"yr":dt:tsim .|> u"yr")
end

"""
    FeederConfiguration

Structure to hold feeder configuration data.
"""
struct FeederConfiguration
    lengths::Dict{String, UnitRange}
    protection_devices::Vector{String}
    
    function FeederConfiguration(lengths::Dict{String, <:Any}, protection_devices::Vector{String})
        # Validate input
        if isempty(lengths)
            throw(ArgumentError("Feeder lengths cannot be empty"))
        end
        if isempty(protection_devices)
            throw(ArgumentError("Protection devices list cannot be empty"))
        end
        
        new(lengths, protection_devices)
    end
end

"""
    create_default_feeder_config() -> FeederConfiguration

Create default feeder configuration for the LVDC network.

# Returns
- `FeederConfiguration`: Default configuration with 5 feeders and 5 protection device types
"""
function create_default_feeder_config()
    # Define feeder lengths (more descriptive names and consistent structure)
    feeder_lengths = Dict(
        "C1" => 1u"m":1u"m":100u"m",  # Feeder 1: 1-100m
        "C2" => 1u"m":1u"m":200u"m",  # Feeder 2: 1-200m
        "C3" => 1u"m":1u"m":150u"m",  # Feeder 3: 1-150m
        "C4" => 1u"m":1u"m":50u"m",   # Feeder 4: 1-50m
        "C5" => 1u"m":1u"m":100u"m"   # Feeder 5: 1-100m
    )
    
    # Define protection device types
    protection_devices = ["SSCB", "MCCB", "HCB", "Fuse", "Fuse_MCCB"]
    
    return FeederConfiguration(feeder_lengths, protection_devices)
end

"""
    process_std_data(raw_data::Dict, time_vector::Vector, config::FeederConfiguration) -> Dict

Process raw STD data into organized solved STD objects.

# Arguments
- `raw_data::Dict`: Raw deserialized STD data
- `time_vector::Vector`: Time vector for the analysis
- `config::FeederConfiguration`: Feeder configuration

# Returns
- `Dict`: Processed STD data organized by protection device and feeder
"""
function process_std_data(raw_data::Dict, time_vector::Vector, config::FeederConfiguration)
    processed_data = Dict{String, Dict{String, Any}}()
    
    for (category, category_data) in raw_data
        processed_category = Dict{String, Any}()
        
        for (device_type, device_data) in category_data
            if !haskey(device_data, :prob) || !haskey(device_data, :power)
                @warn "Missing required keys (prob or power) for $category/$device_type"
                continue
            end
            
            try
                processed_category[device_type] = _MSS.solvedSTD(
                    prob = device_data[:prob],
                    time = time_vector,
                    power = device_data[:power]
                )
            catch e
                @error "Failed to create solvedSTD for $category/$device_type: $e"
                continue
            end
        end
        
        processed_data[category] = processed_category
    end
    
    return processed_data
end

"""
    validate_processed_data(data::Dict, config::FeederConfiguration) -> Bool

Validate that processed data contains expected structure.

# Arguments
- `data::Dict`: Processed STD data
- `config::FeederConfiguration`: Expected configuration

# Returns
- `Bool`: True if validation passes, false otherwise
"""
function validate_processed_data(data::Dict, config::FeederConfiguration)
    # Check if data is not empty
    if isempty(data)
        @error "Processed data is empty"
        return false
    end
    
    # Basic structure validation
    for category in keys(data)
        if !isa(data[category], Dict)
            @error "Category $category is not a dictionary"
            return false
        end
        
        for (device, std_obj) in data[category]
            if !isa(std_obj, _MSS.solvedSTD)
                @error "Object for $category/$device is not a solvedSTD"
                return false
            end
        end
    end
    
    @info "Data validation passed successfully"
    return true
end
