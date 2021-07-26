################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

@testset "Dependence" begin
    
    @testset "Evaluation Dependence" begin
        # Example from: xxx
        ntwʰᵛ = include(joinpath(_MSS.BASE_DIR,"test/networks/high_voltage_system.jl"))
        solve!(ntwʰᵛ)

        #ntwˡᵛ = include(joinpath(_MSS.BASE_DIR,"test/networks/low_voltage_system"))
    end

end