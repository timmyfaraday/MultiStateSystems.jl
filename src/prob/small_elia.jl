################################################################################
#  Copyright 2021, Tom Van Acker, BASF                                         #
################################################################################
# The contents and works on these pages compiled by BASF are subject to        #
# copyright law. Copying, processing, distribution and any kind of use outside #
# the limits of copyright law require the written consent of BASF.             #
################################################################################

# load pkgs 
using Unitful 
using Measurements
using MultiStateSystems

### INPUT 

global DIR = "C:/Users/VANACKT5/OneDrive - BASF/Work/Availability/Code/basf_power_system/"

# availability of corresponding sources of the Elia network in 2004
global Aᵗˢᵒ = 0.999 
global Uᵗˢᵒ = 1.0 - Aᵗˢᵒ

# reliability and maintainability data
global λᵇʳ, μᵇʳ = 1.0u"1/yr", 1.0u"1/wk"
global λᵇᵇ, μᵇᵇ = 1.0u"1/yr", 1.0u"1/wk"
global λᵗʳ, μᵗʳ = 1.0u"1/yr", 1.0u"1/wk"
global λᵖʸ, μᵖʸ = 1.0u"1/yr", 1.0u"1/wk"
global λˢʷ, μˢʷ = 1.0u"1/yr", 1.0u"1/wk"

global λᶜᵇˡ, μᶜᵇˡ = 1.0u"1/(yr*km)", 1.0u"1/wk"

# include the necessary subnetworks
ntwⁿ  = include(joinpath(DIR,"sections/north_2004.jl"))
ntwᵉ  = include(joinpath(DIR,"sections/east_2004.jl"))
ntwʷ  = include(joinpath(DIR,"sections/west_2004.jl"))
ntwᶻᵛ = include(joinpath(DIR,"sections/zandvliet_at1_2004.jl"))
ntwᵛᵗ = include(joinpath(DIR,"sections/basf_vt_2004.jl"))
# include the necessary component state-transition diagrams
stdᵇᵇ = include(joinpath(DIR,"elements/busbar.jl"))
stdᵇʳ = include(joinpath(DIR,"elements/breaker.jl"))
stdˢʷ = include(joinpath(DIR,"elements/switch.jl"))

## MODEL

# initialize the network corresponding with the Elia network in 2004
ntwᵗˢᵒ = Network()

# add the sources to the network 

## east
add_source!(ntwᵗˢᵒ, node = 1, 
                    name = "Ekeren 150kV",
                    ugf  = UGF(:power, [0.0u"MW", 280.0u"MW"], [Uᵗˢᵒ, Aᵗˢᵒ]))
add_source!(ntwᵗˢᵒ, node = 2,
                    name = "Lillo 150kV",
                    ugf  = UGF(:power, [0.0u"MW", 280.0u"MW"], [Uᵗˢᵒ, Aᵗˢᵒ]))

# add the components to the network
## high voltage lines
add_components!(ntwᵗˢᵒ, edge = [(1,3), (1,4), (2,3), (2,3)],
                        name = ["150.116.a", "150.116.b", "150.117.a", "150.117.b"],
                        ntw  = [(ntwᵉ,1), (ntwᵉ,2), (ntwᵉ,3), (ntwᵉ,4)],
                        eval_dep = true)

## zandvliet 150kV
add_components!(ntwᵗˢᵒ, node = [3, 4],
                        name = ["Zandvliet 150kV bb.a", "Zandvliet 150kV bb.b"],
                        std  = stdᵇᵇ)
add_components!(ntwᵗˢᵒ, edge = [(3,5)],
                        name = ["150.111.a"],
                        std  = stdˢʷ)
add_components!(ntwᵗˢᵒ, edge = [(4,5)],
                        name = ["150.111.b"],
                        std  = stdˢʷ)
add_component!(ntwᵗˢᵒ,  edge = (3, 4),
                        name = "Zandvliet 150kV br.ab",
                        std  = stdᵇʳ)
## zandvliet 150kV - BASF 35kV
add_components!(ntwᵗˢᵒ, edge = [(5,6)],
                        name = ["150.111"],
                        ntw  = (ntwᵛᵗ,1))

# add the users to the network 
add_users!(ntwᵗˢᵒ, node = [6],
                   name = ["BASF vt.1"],
                   eval_dep = true)

solve!(ntwᵗˢᵒ)