################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
    Simple Example: Data Loading and Processing

This example demonstrates the improved data loading and processing workflow
for protection device reliability analysis.

This replaces the original io.jl with a cleaner, more structured approach.
"""

# Add the source directory to the path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "config"))

using MultiStateSystems
using Unitful

# Load our improved modules
include("../config/config.jl")
include("../src/data_loader.jl")
include("../src/network_analysis.jl")

const _MSS = MultiStateSystems

function example_data_loading()
    println("=== Example: Improved Data Loading Workflow ===")
    
    # 1. Load configuration
    println("\n1. Loading configuration...")
    config = get_default_config()
    print_config(config)
    
    # 2. Create time vector using the improved function
    println("\n2. Creating time vector...")
    time_vector = create_time_vector(
        tsim = config.simulation.total_time,
        dt = config.simulation.time_step
    )
    println("Created time vector with $(length(time_vector)) points")
    
    # 3. Create feeder configuration
    println("\n3. Setting up feeder configuration...")
    feeder_config = create_default_feeder_config()
    println("Configured $(length(feeder_config.feeders)) feeders:")
    for (name, range) in feeder_config.feeders
        println("  $name: $(first(range)) to $(last(range))")
    end
    println("Protection devices: $(join(feeder_config.protection_devices, ", "))")
    
    # 4. Attempt to load STD data (with error handling)
    println("\n4. Loading STD data...")
    
    # Construct the path to the original data file
    original_data_path = joinpath(_MSS.BASE_DIR, "examples/lvdc/Short-circuit/data/std_s_data_data.dat")
    
    if isfile(original_data_path)
        println("Found original data file: $original_data_path")
        
        try
            # Load the data using our improved function
            raw_data = load_std_data(original_data_path)
            println("Successfully loaded STD data")
            
            # Process the data
            processed_data = process_std_data(raw_data, time_vector, feeder_config)
            
            # Validate the processed data
            if validate_processed_data(processed_data, feeder_config)
                println("Data processing and validation successful!")
                
                # Display summary of processed data
                println("\nProcessed data summary:")
                for (category, category_data) in processed_data
                    println("  Category: $category")
                    for (device, std_obj) in category_data
                        println("    Device: $device ($(typeof(std_obj)))")
                    end
                end
            else
                println("Data validation failed!")
            end
            
        catch e
            println("Error loading or processing data: $e")
            println("This is expected if the original data file structure differs from our assumptions")
        end
    else
        println("Original data file not found at: $original_data_path")
        println("This example shows the improved structure even without the actual data file")
    end
    
    # 5. Demonstrate the improved structure compared to original
    println("\n5. Test MultiStateSystems distribution usage:")
    
    # Test LogNormal distribution from MultiStateSystems
    try
        # Create a repair time distribution using MultiStateSystems
        repair_mean = log(14)u"d"
        repair_std = 0.1u"d"
        
        # Test the MultiStateSystems LogNormal distribution
        repair_dist = _MSS.LogNormal(ustrip(u"d", repair_mean), ustrip(u"d", repair_std))
        println("   ✓ Successfully created MultiStateSystems LogNormal distribution")
        
        # Test the ccdf function
        test_time = 10.0  # 10 days (unitless since LogNormal expects stripped units)
        ccdf_value = _MSS.ccdf(repair_dist, test_time)
        println("   ✓ Successfully computed CCDF: $(round(ccdf_value, digits=6))")
        
        # Test the battery analysis functionality
        network_availability = 0.999 * ones(100)  # Mock network availability
        backup_time = 10.0u"hr"
        time_step = 0.5u"d"
        
        battery_analysis = analyze_battery_system(
            network_availability,
            backup_time,
            repair_mean,
            repair_std,
            time_step
        )
        
        avg_battery_availability = mean(battery_analysis.availability_over_time)
        println("   ✓ Battery analysis completed. Average availability: $(round(avg_battery_availability, digits=6))")
        
    catch e
        println("   ✗ Error testing MultiStateSystems distributions: $e")
    end
    
    println("\n6. Improvements over original io.jl:")
    println("   ✓ Separated concerns into logical modules")
    println("   ✓ Added comprehensive error handling")
    println("   ✓ Created reusable configuration system") 
    println("   ✓ Added input validation and documentation")
    println("   ✓ Structured data types for better organization")
    println("   ✓ Eliminated hard-coded values")
    println("   ✓ Added logging and progress information")
    println("   ✓ Uses MultiStateSystems' native distributions (no external Distributions.jl dependency)")
    
    println("\n=== Example Complete ===")
end

# Demonstrate individual functions with MultiStateSystems distributions
function demonstrate_individual_functions()
    println("\n=== Demonstrating Individual Improved Functions ===")
    
    # Show the improved time vector creation
    println("\n1. Flexible time vector creation:")
    
    # Default parameters
    time1 = create_time_vector()
    println("  Default: $(length(time1)) points over 25 years")
    
    # Custom parameters
    time2 = create_time_vector(tsim=10.0u"yr", dt=1.0u"d")
    println("  Custom: $(length(time2)) points over 10 years with 1-day steps")
    
    # Show configuration validation
    println("\n2. Configuration validation:")
    valid_config = get_default_config()
    is_valid = validate_config(valid_config)
    println("  Default config validation: $is_valid")
    
    # Show feeder configuration
    println("\n3. Structured feeder configuration:")
    feeder_config = create_default_feeder_config()
    println("  Feeders: $(length(feeder_config.feeders))")
    println("  Protection devices: $(length(feeder_config.protection_devices))")
    
    # Test MultiStateSystems distributions
    println("\n4. MultiStateSystems distribution testing:")
    try
        # Test different MultiStateSystems distribution constructors
        exp_dist = _MSS.Exponential(100.0u"hr")
        weibull_dist = _MSS.Weibull(200.0u"hr", 2.0)
        lognormal_dist = _MSS.LogNormal(5.0, 1.0)
        
        println("  Exponential distribution created: $(typeof(exp_dist))")
        println("  Weibull distribution created: $(typeof(weibull_dist))")
        println("  LogNormal distribution created: $(typeof(lognormal_dist))")
        
        # Test distribution functions
        test_point = 50.0
        println("  LogNormal CCDF at t=50: $(round(_MSS.ccdf(lognormal_dist, test_point), digits=6))")
        println("  LogNormal PDF at t=50: $(round(_MSS.pdf(lognormal_dist, test_point), digits=6))")
        
    catch e
        println("  Error testing distributions: $e")
    end
    
    println("\n=== Function Demonstrations Complete ===")
end

# Run the examples
if abspath(PROGRAM_FILE) == @__FILE__
    example_data_loading()
    demonstrate_individual_functions()
end
