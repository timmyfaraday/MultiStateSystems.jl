using MultiStateSystems
using Serialization
using Unitful

const _MSS = MultiStateSystems

function read_std_s_data(file_path::String)
    if !isfile(file_path)
        error("File does not exist: $file_path")
    end
    return deserialize(file_path)
end

# =============================================================================
# NOTICE: This file has been improved and restructured!
# 
# The improved version is available in:
# ../Protection_device_journal/
#
# Key improvements:
# - Modular structure with separated concerns
# - Comprehensive error handling and validation
# - Configurable parameters instead of hard-coded values
# - Better documentation and type safety
# - Structured data types and workflows
#
# To use the improved version:
# 1. Navigate to ../Protection_device_journal/
# 2. Run: julia examples/improved_io_example.jl
# 3. Or run: julia src/main_analysis.jl (requires data migration)
#
# For data migration, run: julia data_migration.jl
# =============================================================================

# Original example usage (kept for backward compatibility)
input_file = joinpath(_MSS.BASE_DIR, "examples/lvdc/Short-circuit/data/std_s_data_data.dat")

# Improved version uses structured error handling
try
    std_s_data = read_std_s_data(input_file)
    println("Data successfully loaded from $input_file")
catch e
    @warn "Data loading failed: $e"
    @info "Consider using the improved version in ../Protection_device_journal/"
end

# Simulation parameters (consider using the improved configuration system)
tsim = 25.0u"yr"  # Simulation time
dt = 0.5u"d"     # Time step
time = 0.0u"yr":dt:tsim .|> u"yr"  # Time vector in years

# Feeder configuration (improved version uses FeederConfiguration struct)
L_c = Dict("C1" => 1u"m":1u"m":100u"m", 
           "C2" => 1u"m":1u"m":200u"m",
           "C3" => 1u"m":1u"m":150u"m",
           "C4" => 1u"m":1u"m":50u"m",
           "C5" => 1u"m":1u"m":100u"m")

# Protection device mapping (improved version uses NetworkConfig)
L_tot = Dict("SSCB" => L_c,
             "MCCB" => L_c,
             "HCB"  => L_c,
             "Fuse" => L_c,
             "Fuse_MCCB" => L_c)

# Data processing (improved version uses process_std_data() function)
std_s = Dict()
if @isdefined(std_s_data)
    for (key, value) in std_s_data
        std_s_i = Dict()
        for (cb, std_sol) in value
            try
                prob = std_sol[:prob]  # Extract the problem definition
                power = std_sol[:power]  # Extract the power data
                std_s_i[cb] = solvedSTD(prob = prob,
                                        time = collect(time),
                                        power = power)
            catch e
                @warn "Failed to process STD for $key/$cb: $e"
            end
        end
        std_s[key] = std_s_i
    end
    
    println("Processed $(length(std_s)) categories of STD data")
    @info "For improved data processing, see ../Protection_device_journal/src/data_loader.jl"
else
    @warn "STD data not available - using empty dictionary"
end