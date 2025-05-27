################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models, often found in           #
# reliability engineering.                                                     #
# See https://github.com/timmyfaraday/MultiStateSystems.jl                     #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

# functions ####################################################################
## general
nu(ntw::AbstractNetwork) = length(ntw.usr)
nu(ntw::AbstractNetwork,u_node::Int) = length(ntw.ulib[u_node])
usr(ntw::AbstractNetwork) = ntw.usr
usr(ntw::AbstractNetwork,u_node::Int) = ntw.usr[ntw.ulib[u_node]]
usr_ids(ntw::AbstractNetwork) = 1:nu(ntw)
usr_ids(ntw::AbstractNetwork,u_node::Int) = ntw.ulib[u_node]
usr_nodes(ntw::AbstractNetwork) = keys(ntw.ulib)

## indices
"""
    EENS(usr::MultiStateSystems.PropDict)

Expected Energy Not Served (EENS) [MWh]

`EENS(usr)` gives the EENS when an user's demand equals the maximal output of 
the system for that user.
"""
EENS(usr::PropDict) =
    8760u"hr"*sum((maximum(usr[:ugf].val).-usr[:ugf].val).*usr[:ugf].prb) |> u"MWh"

"""
    GRA(usr::MultiStateSystems.PropDict)

Generation Ratio Availability (GRA) [-]

- `GRA(usr,GR)` gives the probability of transferring at least the GR to a
  user through the network.
- `GRA(usr)` gives a user's GRA for a GR ranging from 0.0 to 1.0.
- `GRO(usr,GR)` gives the output towards a user for a given GR
"""
GRA(usr::PropDict) = [GRA(usr,GR) for GR in 0.0:0.01:1.0]
GRA(usr::PropDict,GR::Float64) =
    sum(usr[:ugf].prb[findfirst(x -> x .>= GRO(usr,GR),usr[:ugf].val):end])
GRO(usr::PropDict,GR::Float64) = GR*maximum(usr[:ugf].val)

## user
"""
    add_user!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)

Adds a single user to the network `ntw` and fills their corresponding
`PropDict` with the named arguments `kwargs`.

# Example
```julia-repl
julia> ntwᵖʷʳ = ntw()
julia> add_source!(ntwᵖʷʳ, node = 1,
                           ind  = [:EENS])
```
"""
function add_user!(ntw::AbstractNetwork; kwargs...)
    haskey(kwargs,:node) || return false
    
    add_vertex!(ntw, kwargs[:node])

    info = UserInfo()

    push!(ntw.usr, PropDict(:info => info, kwargs...))
    update_lib!(:node, ntw.usr, ntw.ulib)
    return true
end
"""
    add_users!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)

Adds multiple users to the network `ntw` and fills their corresponding
`PropDict` with the named arguments `kwargs`. Either an uniform arguments is
given which holds for all components or an array is given whith specific
argument for each component.

# Example
```julia-repl
julia> ntwᵖʷʳ = ntw()
julia> add_sources!(ntwᵖʷʳ, node = [1,5,8],
                            ind  = [:EENS])
```
"""
function add_users!(ntw::AbstractNetwork; kwargs...)
    (test(kwargs) && haskey(kwargs,:node)) || return false

    node, Nn = kwargs[:node], length(kwargs[:node])
    
    eval_dep_ids = init_eval_dep_usr(ntw, kwargs, Nn)

    for ni in 1:Nn
        prop_dict = reduce(kwargs, ni, excl=[:node])
        set_eval_dep!(prop_dict, ni, eval_dep_ids)
        
        add_user!(ntw, node[ni], prop_dict...)
    end
    return true
end
