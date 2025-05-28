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

# Example usage
input_file = joinpath(_MSS.BASE_DIR, "examples/lvdc/Short-circuit/data/std_s_data_data.dat")
std_s_data = read_std_s_data(input_file)
println("Data successfully loaded from $input_file")

tsim = 25.0u"yr"  # Simulation time
dt = 0.5u"d"     # Time step
time = 0.0u"yr":dt:tsim .|> u"yr"  # Time vector in years

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
             "Fuse" => L_c,
             "Fuse_MCCB" => L_c)

std_s = Dict()
for (key, value) in std_s_data
    std_s_i = Dict()
    for (cb, std_sol) in value
        prob = std_sol[:prob]  # Extract the problem definition
        power = std_sol[:power]  # Extract the power data
        std_s_i[cb] = solvedSTD(prob = prob,
                                time = collect(time),
                                power = power)
    end
    std_s[key] = std_s_i
end