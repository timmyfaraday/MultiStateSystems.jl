################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
    Protection Device Reliability Analysis

Main execution script for analyzing protection device reliability in LVDC networks.
This script demonstrates the improved structure and workflow for the analysis.

Usage:
    julia main_analysis.jl

The script will:
1. Load configuration
2. Load and process STD data
3. Perform battery backup analysis
4. Calculate bus availability for different scenarios
5. Save results and generate summary
"""

# Load required modules
include("../config/config.jl")
include("data_loader.jl")
include("network_analysis.jl")

using MultiStateSystems
using Serialization
using Unitful
using Dates

const _MSS = MultiStateSystems

"""
    main()

Main execution function for the protection device analysis.
"""
function main()
    println("=== Protection Device Reliability Analysis ===")
    println("Starting analysis at $(now())")
    
    try
        # 1. Load and validate configuration
        println("\n1. Loading configuration...")
        config = get_default_config()
        
        if !validate_config(config)
            error("Configuration validation failed")
        end
        
        print_config(config)
        
        # 2. Create time vector
        println("\n2. Creating time vector...")
        time_vector = create_time_vector(
            tsim = config.simulation.total_time,
            dt = config.simulation.time_step
        )
        println("Time vector created: $(length(time_vector)) points over $(config.simulation.total_time)")
        
        # 3. Load STD data
        println("\n3. Loading STD data...")
        if !isfile(config.paths.std_data_file)
            error("STD data file not found: $(config.paths.std_data_file)")
        end
        
        raw_std_data = load_std_data(config.paths.std_data_file)
        println("Raw STD data loaded successfully")
        
        # 4. Process STD data
        println("\n4. Processing STD data...")
        feeder_config = create_default_feeder_config()
        processed_std_data = process_std_data(raw_std_data, time_vector, feeder_config)
        
        if !validate_processed_data(processed_std_data, feeder_config)
            error("Processed data validation failed")
        end
        
        println("STD data processed and validated successfully")
        
        # 5. Create AC network model
        println("\n5. Creating AC network model...")
        ac_network, ac_availability = create_ac_network(config, time_vector)
        println("AC network model created")
        
        # 6. Analyze battery backup system
        println("\n6. Analyzing battery backup system...")
        battery_analysis = analyze_battery_system(
            ac_availability,
            config.simulation.battery_backup_time,
            config.simulation.repair_mean,
            config.simulation.repair_std,
            config.simulation.time_step
        )
        
        avg_battery_availability = mean(battery_analysis.availability_over_time)
        println("Battery backup analysis complete. Average availability: $(round(avg_battery_availability, digits=6))")
        
        # 7. Calculate bus availability scenarios
        println("\n7. Calculating bus availability scenarios...")
        
        # Scenario 1: All feeders
        println("  - All feeders scenario...")
        bus_availability_all = calculate_bus_availability(
            processed_std_data, 
            time_vector, 
            config
        )
        
        # Scenario 2: Individual feeder analyses
        individual_analyses = Dict()
        for feeder in keys(config.network.feeders)
            println("  - Single feeder scenario: $feeder...")
            individual_analyses[feeder] = calculate_selective_bus_availability(
                processed_std_data,
                time_vector,
                config,
                [feeder]
            )
        end
        
        # Scenario 3: Excluding specific feeders
        excluded_scenarios = Dict()
        for feeder in keys(config.network.feeders)
            println("  - Excluding feeder: $feeder...")
            excluded_scenarios[feeder] = calculate_bus_availability(
                processed_std_data,
                time_vector,
                config,
                excluded_feeders = [feeder]
            )
        end
        
        # 8. Generate analysis summary
        println("\n8. Generating analysis summary...")
        
        println("\n--- All Feeders Analysis ---")
        print_analysis_summary(bus_availability_all)
        
        for (feeder, results) in individual_analyses
            println("\n--- Single Feeder Analysis: $feeder ---")
            print_analysis_summary(results)
        end
        
        # 9. Save results
        println("\n9. Saving results...")
        results_package = Dict(
            "config" => config,
            "time_vector" => time_vector,
            "battery_analysis" => battery_analysis,
            "bus_availability_all" => bus_availability_all,
            "individual_analyses" => individual_analyses,
            "excluded_scenarios" => excluded_scenarios,
            "timestamp" => now()
        )
        
        serialize(config.paths.results_file, results_package)
        println("Results saved to: $(config.paths.results_file)")
        
        # 10. Analysis complete
        println("\n=== Analysis Complete ===")
        println("Total analysis time: $(now() - start_time) seconds")
        println("Results available in: $(config.paths.results_dir)")
        
    catch e
        @error "Analysis failed: $e"
        rethrow(e)
    end
end

"""
    load_previous_results(config) -> Dict

Load previously saved analysis results.

# Arguments
- `config`: Configuration object

# Returns
- `Dict`: Previously saved results
"""
function load_previous_results(config)
    if !isfile(config.paths.results_file)
        error("Results file not found: $(config.paths.results_file)")
    end
    
    return deserialize(config.paths.results_file)
end

"""
    run_sensitivity_analysis(config, parameter_variations::Dict)

Run sensitivity analysis by varying key parameters.

# Arguments
- `config`: Base configuration
- `parameter_variations::Dict`: Parameters to vary and their ranges
"""
function run_sensitivity_analysis(config, parameter_variations::Dict)
    println("=== Sensitivity Analysis ===")
    
    # This is a placeholder for sensitivity analysis implementation
    # Would vary parameters like battery backup time, repair times, etc.
    
    @warn "Sensitivity analysis not yet implemented - placeholder function"
end

# Store start time for performance measurement
const start_time = now()

# Execute main analysis if script is run directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
