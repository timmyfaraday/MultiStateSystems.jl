################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
    Data Migration Utility

Utility script to copy data files from the original Short-circuit example
to the new Protection_device_journal structure.
"""

using MultiStateSystems
using Serialization

const _MSS = MultiStateSystems

"""
    copy_original_data()

Copy data files from the original Short-circuit example to the new structure.
"""
function copy_original_data()
    println("=== Data Migration Utility ===")
    
    # Define paths
    original_base = joinpath(_MSS.BASE_DIR, "examples/lvdc/Short-circuit")
    new_base = dirname(@__DIR__)
    
    # Create data directory if it doesn't exist
    new_data_dir = joinpath(new_base, "data")
    if !isdir(new_data_dir)
        mkpath(new_data_dir)
        println("Created data directory: $new_data_dir")
    end
    
    # Files to copy
    files_to_copy = [
        "data/std_s_data_data.dat",
        "data/source_data.dat", 
        "data/results.dat",
        "data/std_s_data.dat"
    ]
    
    copied_files = 0
    
    for file_path in files_to_copy
        original_file = joinpath(original_base, file_path)
        new_file = joinpath(new_base, file_path)
        
        if isfile(original_file)
            try
                cp(original_file, new_file, force=true)
                println("✓ Copied: $file_path")
                copied_files += 1
            catch e
                println("✗ Failed to copy $file_path: $e")
            end
        else
            println("✗ Original file not found: $original_file")
        end
    end
    
    println("\nCopied $copied_files out of $(length(files_to_copy)) files")
    
    # Also copy any Julia files that might contain data definitions
    julia_files = [
        "data/stds.jl"
    ]
    
    for file_path in julia_files
        original_file = joinpath(original_base, file_path)
        new_file = joinpath(new_base, file_path)
        
        if isfile(original_file)
            try
                cp(original_file, new_file, force=true)
                println("✓ Copied Julia file: $file_path")
            catch e
                println("✗ Failed to copy Julia file $file_path: $e")
            end
        end
    end
    
    println("\n=== Migration Complete ===")
    println("Data files are now available in: $new_data_dir")
    println("You can now run the improved analysis scripts.")
end

"""
    create_sample_data()

Create sample data files for testing when original data is not available.
"""
function create_sample_data()
    println("=== Creating Sample Data ===")
    
    new_base = dirname(@__DIR__)
    new_data_dir = joinpath(new_base, "data")
    
    if !isdir(new_data_dir)
        mkpath(new_data_dir)
    end
    
    # Create a simple sample data structure
    sample_std_data = Dict(
        "Source" => Dict(
            "SSCB" => Dict(
                "C1" => Dict(:prob => [0.999, 0.001], :power => [50.0, 0.0]),
                "C2" => Dict(:prob => [0.999, 0.001], :power => [50.0, 0.0]),
                "C3" => Dict(:prob => [0.999, 0.001], :power => [50.0, 0.0]),
                "C4" => Dict(:prob => [0.999, 0.001], :power => [50.0, 0.0]),
                "C5" => Dict(:prob => [0.999, 0.001], :power => [50.0, 0.0])
            ),
            "MCCB" => Dict(
                "C1" => Dict(:prob => [0.998, 0.002], :power => [50.0, 0.0]),
                "C2" => Dict(:prob => [0.998, 0.002], :power => [50.0, 0.0]),
                "C3" => Dict(:prob => [0.998, 0.002], :power => [50.0, 0.0]),
                "C4" => Dict(:prob => [0.998, 0.002], :power => [50.0, 0.0]),
                "C5" => Dict(:prob => [0.998, 0.002], :power => [50.0, 0.0])
            )
        )
    )
    
    # Save sample data
    sample_file = joinpath(new_data_dir, "sample_std_data.dat")
    serialize(sample_file, sample_std_data)
    
    println("✓ Created sample data file: $sample_file")
    println("This can be used for testing the improved structure")
    
    # Create a simple configuration file that points to the sample data
    config_override = """
# Configuration override for sample data
# Modify config.jl to use this file path:
# std_data_file = joinpath(data_dir, "sample_std_data.dat")
"""
    
    config_file = joinpath(new_data_dir, "sample_config_override.txt")
    open(config_file, "w") do f
        write(f, config_override)
    end
    
    println("✓ Created configuration notes: $config_file")
    println("\n=== Sample Data Creation Complete ===")
end

"""
    verify_data_availability()

Check what data files are available and provide guidance.
"""
function verify_data_availability()
    println("=== Data Availability Check ===")
    
    original_base = joinpath(_MSS.BASE_DIR, "examples/lvdc/Short-circuit")
    new_base = dirname(@__DIR__)
    
    println("Original Short-circuit directory: $original_base")
    println("New Protection_device_journal directory: $new_base")
    
    # Check original data
    println("\nOriginal data files:")
    original_data_files = [
        "data/std_s_data_data.dat",
        "data/source_data.dat",
        "data/results.dat"
    ]
    
    original_available = 0
    for file_path in original_data_files
        full_path = joinpath(original_base, file_path)
        if isfile(full_path)
            println("  ✓ $file_path")
            original_available += 1
        else
            println("  ✗ $file_path")
        end
    end
    
    # Check new data
    println("\nNew data directory:")
    new_data_dir = joinpath(new_base, "data")
    if isdir(new_data_dir)
        data_files = readdir(new_data_dir)
        if isempty(data_files)
            println("  (empty)")
        else
            for file in data_files
                println("  ✓ $file")
            end
        end
    else
        println("  (does not exist)")
    end
    
    # Provide recommendations
    println("\nRecommendations:")
    if original_available > 0
        println("  1. Run copy_original_data() to migrate existing data")
    else
        println("  1. Original data not available")
    end
    println("  2. Run create_sample_data() to create test data")
    println("  3. Use examples/improved_io_example.jl to test the structure")
    
    println("\n=== Check Complete ===")
end

# Run functions based on script arguments or interactively
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) > 0
        action = ARGS[1]
        if action == "copy"
            copy_original_data()
        elseif action == "sample"
            create_sample_data()
        elseif action == "verify"
            verify_data_availability()
        else
            println("Usage: julia data_migration.jl [copy|sample|verify]")
        end
    else
        # Interactive mode
        verify_data_availability()
        println("\nOptions:")
        println("  - Run with 'copy' to copy original data")
        println("  - Run with 'sample' to create sample data")  
        println("  - Run with 'verify' to check data availability")
    end
end
