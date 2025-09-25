################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

using MultiStateSystems
using Unitful
using JLD2

"""
    run_single_bus_example(; save_stds::Bool = true, load_stds_from=nothing)

Run the single bus LVDC network example.

# Arguments
- `save_stds::Bool`: Whether to save solved STDs to file (default: true)
- `load_stds_from`: Path to pre-solved STD file to load (default: "MultiStateSystems.jl/examples/lvdc/ICDCM_2025/results/single_bus_stds.jld2", because the STDs are big to solve)

# Returns
- `Dict`: Dictionary containing the complete analysis results
"""
function run_single_bus_example(; save_stds::Bool = true, load_stds_from="MultiStateSystems.jl/examples/lvdc/ICDCM_2025/results/single_bus_stds.jld2")
    println("Running Single Bus LVDC Example...")
    println("="^50)
    
    # Get system configuration
    sys_params = get_default_system_parameters()
    int_params = get_default_integration_parameters()
    
    # Validate configuration
    validate_configuration(sys_params, int_params)
    
    # Run the example
    results = single_bus_analysis(sys_params, int_params; load_stds_from=load_stds_from)
    
    # Print summary
    print_summary(results)

    if save_stds && isnothing(load_stds_from)
        stds_path = "MultiStateSystems.jl/examples/lvdc/ICDCM_2025/results/single_bus_stds.jld2"
        save_solved_stds(results["solved_stds"], stds_path)
    end
    
    println("Single bus example completed successfully!")
    return results
end

"""
    run_double_bus_example(; save_stds::Bool = true, load_stds_from=nothing)

Run the double bus LVDC network example.

# Arguments  
- `save_stds::Bool`: Whether to save solved STDs to file (default: true)
- `load_stds_from`: Path to pre-solved STD file to load (default: "MultiStateSystems.jl/examples/lvdc/ICDCM_2025/results/double_bus_stds.jld2", because the STDs are big to solve)

# Returns
- `Dict`: Dictionary containing the complete analysis results
"""
function run_double_bus_example(; save_stds::Bool = true, load_stds_from="MultiStateSystems.jl/examples/lvdc/ICDCM_2025/results/double_bus_stds.jld2")
    println("Running Double Bus LVDC Example...")
    println("="^50)
    
    # Get system configuration
    sys_params = get_default_system_parameters()
    int_params = get_default_integration_parameters()
    
    # Validate configuration
    validate_configuration(sys_params, int_params)
    
    # Run the example
    results = double_bus_analysis(sys_params, int_params; load_stds_from=load_stds_from)

    results_no_bridge = double_bus_analysis(sys_params, int_params; include_bridge=false, load_stds_from=load_stds_from)
    
    if save_stds && isnothing(load_stds_from)
        stds_path = "MultiStateSystems.jl/examples/lvdc/ICDCM_2025/results/double_bus_stds.jld2"
        save_solved_stds(results["solved_stds"], stds_path)
    end

    println("Entire double bus example completed successfully!")
    return results, results_no_bridge
end

"""
    run_complete_analysis(; save_stds::Bool = true, load_stds_from=nothing)

Run both single bus and double bus examples with complete analysis.

# Arguments
- `save_stds::Bool`: Whether to save solved STDs to file (default: true)
- `load_stds_from`: Path to pre-solved STD file to load (default: nothing, solves from scratch)

# Returns
- `Dict`: Dictionary containing results from both examples
"""
function run_complete_analysis(; save_stds::Bool = true, load_single_bus_stds_from="MultiStateSystems.jl/examples/lvdc/ICDCM_2025/results/single_bus_stds.jld2", load_double_bus_stds_from="MultiStateSystems.jl/examples/lvdc/ICDCM_2025/results/double_bus_stds.jld2")
    println("Running Complete LVDC Analysis...")
    println("="^60)
    
    results = Dict()
    
    # Run single bus example
    results["single_bus"] = run_single_bus_example(save_stds=save_stds, load_stds_from=load_single_bus_stds_from)
    
    println("\n")

    # Run double bus example
    results["double_bus"], results["double_bus_no_bridge"] = run_double_bus_example(save_stds=save_stds, load_stds_from=load_double_bus_stds_from)

    # Print combined summary
    println("\n" * "="^60)
    println("COMPLETE ANALYSIS SUMMARY")
    println("="^60)
    println("Single Bus Analysis: ✓ Completed")
    println("Double Bus Analysis: ✓ Completed")
    
    println("\nComplete analysis finished successfully!")
    return results
end

"""
    generate_example_stds()

Generate example STD files for repository inclusion.
This creates lightweight files suitable for version control that contain
only the solved state transition diagrams.
"""
function generate_example_stds()
    println("Generating Example STD Files for Repository...")
    println("="^50)
    
    # Generate STD files for both configurations
    println("Generating single bus STDs...")
    run_single_bus_example(save_stds=true, load_stds_from=nothing)
    
    println("\nGenerating double bus STDs...")
    run_double_bus_example(save_stds=true, load_stds_from=nothing)
    
    println("\nGenerated example files:")
    println("• results/single_bus_stds.jld2")
    println("• results/double_bus_stds.jld2")
    println("\nThese files contain only the solved state transition diagrams")
    println("and are suitable for including in the repository.")
end

"""
    demo_load_from_stds()

Demonstrate loading pre-solved STDs and running network analysis.
This shows how to use the pre-solved STD files to skip the STD generation step.
"""
function demo_load_from_stds()
    println("Demonstration: Loading Pre-solved STDs")
    println("="^40)
    
    # First, ensure we have STD files
    single_stds_file = "MultiStateSystems.jl/examples/lvdc/ICDCM_2025/results/single_bus_stds.jld2"
    double_stds_file = "MultiStateSystems.jl/examples/lvdc/ICDCM_2025/results/double_bus_stds.jld2"
    
    if !isfile(single_stds_file) || !isfile(double_stds_file)
        println("STD files not found. Generating them first...")
        generate_example_stds()
    end
    
    # Demonstrate loading from STDs
    println("\nLoading single bus analysis from pre-solved STDs...")
    single_results = run_single_bus_example(save_stds=false, load_stds_from=single_stds_file)
    
    println("\nLoading double bus analysis from pre-solved STDs...")
    double_results = run_double_bus_example(save_stds=false, load_stds_from=double_stds_file)
    
    println("\n✓ Both analyses completed using pre-solved STDs")
    println("This approach is much faster as it skips STD generation and solving.")
    
    return Dict("single_bus" => single_results, "double_bus" => double_results)
end
