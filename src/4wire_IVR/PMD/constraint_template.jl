"""
    function constraint_mc_transformer_voltage(
        pm::_PMD.ExplicitNeutralModels,
        i::Int;
        nw::Int=nw_id_default,
        fix_taps::Bool=true
    )

For models with explicit neutrals,
links the voltage of the from-side and to-side transformer windings
"""
function constraint_mc_transformer_voltage(pm::_PMD.ExplicitNeutralModels, i::Int; nw::Int=nw_id_default, fix_taps::Bool=true)
    transformer = _PMD.ref(pm, nw, :transformer, i)
    f_bus = transformer["f_bus"]
    t_bus = transformer["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    configuration = transformer["configuration"]
    f_connections = transformer["f_connections"]
    t_connections = transformer["t_connections"]
    tm_set = transformer["tm_set"]
    tm_fixed = fix_taps ? ones(Bool, length(tm_set)) : transformer["tm_fix"]
    tm_scale = _PMD.calculate_tm_scale(transformer, _PMD.ref(pm, nw, :bus, f_bus), _PMD.ref(pm, nw, :bus, t_bus))

    #TODO change data model
    # there is redundancy in specifying polarity seperately on from and to side
    #TODO change this once migrated to new data model
    pol = transformer["polarity"]

    if configuration == _PMD.WYE
        constraint_mc_transformer_voltage_yy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == _PMD.DELTA
        constraint_mc_transformer_voltage_dy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == "zig-zag"
        error("Zig-zag not yet supported.")
    end
end


"""
	function constraint_mc_transformer_current(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
		i::Int;
		nw::Int=nw_id_default,
		fix_taps::Bool=true
	)

For IVR models with explicit neutrals,
links the current variables of the from-side and to-side transformer windings,
and creates expressions for the terminal current flows
"""
function constraint_mc_transformer_current(pm::_PMD.AbstractExplicitNeutralIVRModel, i::Int; nw::Int=nw_id_default, fix_taps::Bool=true)
    # if _PMD.ref(pm, nw_id_default, :conductors)!=3
    #     error("Transformers only work with networks with three conductors.")
    # end

    transformer = _PMD.ref(pm, nw, :transformer, i)
    f_bus = transformer["f_bus"]
    t_bus = transformer["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    configuration = transformer["configuration"]
    f_connections = transformer["f_connections"]
    t_connections = transformer["t_connections"]
    tm_set = transformer["tm_set"]
    tm_fixed = fix_taps ? ones(Bool, length(tm_set)) : transformer["tm_fix"]
    tm_scale = _PMD.calculate_tm_scale(transformer, _PMD.ref(pm, nw, :bus, f_bus), _PMD.ref(pm, nw, :bus, t_bus))

    #TODO change data model
    # there is redundancy in specifying polarity seperately on from and to side
    #TODO change this once migrated to new data model
    pol = transformer["polarity"]

    if configuration == _PMD.WYE
        constraint_mc_transformer_current_yy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == _PMD.DELTA
        constraint_mc_transformer_current_dy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == "zig-zag"
        error("Zig-zag not yet supported.")
    end
end


"""
	function constraint_mc_transformer_thermal_limit(
		pm::_PMD.ExplicitNeutralModels,
		id::Int;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

Imposes a bound on the total apparent at each transformer winding
"""
function constraint_mc_transformer_thermal_limit(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    trans = _PMD.ref(pm, nw, :transformer, id)
    f_bus = trans["f_bus"]
    t_bus = trans["t_bus"]
    f_idx = (id,f_bus,t_bus)
    t_idx = (id,t_bus,f_bus)
    f_conns = trans["f_connections"]
    t_conns = trans["t_connections"]
    config = trans["configuration"]
    sm_ub = trans["sm_ub"]

    constraint_mc_transformer_thermal_limit(pm, nw, id, f_idx, t_idx, f_bus, t_bus, f_conns, t_conns, config, sm_ub)
end



"""
    function constraint_mc_switch_current(
        pm::_PMD.ExplicitNeutralModels,
        id::Int;
        nw::Int=nw_id_default,
        report::Bool=true
    )

For models with explicit neutrals,
link the switch currents or create appropiate expressions for them.
"""
function constraint_mc_switch_current(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, report::Bool=true)
    switch = _PMD.ref(pm, nw, :switch, id)
    f_bus = switch["f_bus"]
    t_bus = switch["t_bus"]
    f_idx = (id, f_bus, t_bus)
    t_idx = (id, t_bus, f_bus)

    constraint_mc_switch_current(pm, nw, id, f_idx, t_idx, switch["f_connections"], switch["t_connections"])
end


"""
    constraint_mc_switch_current_limit(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for switch current limit constraints
"""
function constraint_mc_switch_current_limit(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    switch = _PMD.ref(pm, nw, :switch, i)

    if !haskey(_PMD.con(pm, nw), :mu_cm_switch)
        _PMD.con(pm, nw)[:mu_cm_switch] = Dict{Tuple{Int,Int,Int},Vector{JuMP.ConstraintRef}}()
    end

    if haskey(switch, "current_rating") && any(switch["current_rating"] .< Inf)
        f_idx = (i, switch["f_bus"], switch["t_bus"])
        constraint_mc_switch_current_limit(pm, nw, f_idx, switch["f_connections"], switch["current_rating"])
    end
    nothing
end


"""
    constraint_mc_switch_thermal_limit(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for switch thermal limit constraint
"""
function constraint_mc_switch_thermal_limit(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    switch = _PMD.ref(pm, nw, :switch, i)
    f_idx = (i, switch["f_bus"], switch["t_bus"])

    if !haskey(_PMD.con(pm, nw), :mu_sm_switch)
        _PMD.con(pm, nw)[:mu_sm_switch] = Dict{Tuple{Int,Int,Int},Vector{JuMP.ConstraintRef}}()
    end

    if haskey(switch, "thermal_rating") && any(switch["thermal_rating"] .< Inf)
        constraint_mc_switch_thermal_limit(pm, nw, f_idx, switch["f_connections"], switch["thermal_rating"])
    end
    nothing
end