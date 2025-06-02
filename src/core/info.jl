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

# structs ######################################################################
## network
### cmp
""
mutable struct ComponentInfo{B<:Bool} <: AbstractInfo
    eval_dep::B
    eval_dep_id::B

    # default constructor 
    ComponentInfo() = new{Bool}(false,false)
end
### ntw
""
mutable struct NetworkInfo{B<:Bool} <: AbstractInfo
    solved::B
    dependent_sources::B
    eval_dep::B

    # default constructor
    NetworkInfo() = new{Bool}(false,false,false)
end
### src
""
mutable struct SourceInfo{B<:Bool} <: AbstractInfo
    eval_dep::B
    eval_dep_id::B

    # default constructor
    SourceInfo() = new{Bool}(false,false)
end
### usr
""
mutable struct UserInfo{B<:Bool} <: AbstractInfo 
    eval_dep::B
    eval_dep_id::B

    # default constructor
    UserInfo() = new{Bool}(false,false)
end

## state-transition diagram
### state
""
mutable struct StateInfo{B<:Bool} <: AbstractInfo
    renewal::B
    trapping::B

    # constructor
    StateInfo() = new{Bool}(true,false)
end
### std
""
mutable struct STDInfo{B<:Bool} <: AbstractInfo
    solved::B
    renewal::B
    markovian::B
    time_homogeneous::B

    # constructor
    STDInfo() = new{Bool}(false,true,true,true)
end
### trans
""
mutable struct TransInfo{B<:Bool} <: AbstractInfo
    renewal::B
    markovian::B
    time_homogeneous::B

    # constructor
    TransInfo() = new{Bool}(true,true,true)
end

# functions ####################################################################
## getters
### ntw
get_info(ntw::AbstractNetwork, info::Symbol) = getproperty(ntw.props[:info],info)
get_info(prt::PropDict, info::Symbol) = getproperty(prt[:info],info)
### std
get_info(std::AbstractSTD, info::Symbol) = getproperty(std.props[:info],info)
get_info(std::AbstractSTD, ns::Int, info::Symbol) =
    getproperty(std.sprops[ns][:info],info)
get_info(std::AbstractSTD, nt::Graphs.Edge, info::Symbol) =
    getproperty(std.tprops[nt][:info],info)

## dependence
set_eval_info!(info::AbstractInfo; kwargs...) =
    if haskey(kwargs, :eval_dep)
        for i_key in [:eval_dep,:eval_dep_id]
            set_info!(info, i_key, kwargs[i_key])
            # delete!(prop_dict, i_key)
        end
    end

## setters
### ntw
set_info!(ntw::AbstractNetwork, info::Symbol, value::Bool) =
    setproperty!(ntw.props[:info], info, value)
set_info!(prt::PropDict, info:: Symbol, value::Bool) =
    setproperty!(prt[:info], info, value)
set_info!(info::AbstractInfo, field::Symbol, value::Bool) =
    setproperty!(info, field, value)
### std
set_info!(std::AbstractSTD, info::Symbol, value::Bool) =
    setproperty!(std.props[:info],info,value)
set_info!(std::AbstractSTD, ns::Int, info::Symbol, value::Bool) =
    setproperty!(std.sprops[ns][:info],info,value)
set_info!(std::AbstractSTD, nt::Graphs.Edge, info::Symbol, value::Bool) =
    setproperty!(std.tprops[nt][:info],info,value)  