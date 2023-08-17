"""
    function constraint_mc_voltage_absolute(
        pm::_PMD.RectangularVoltageExplicitNeutralModels,
        id::Int;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true,
    )

Imposes absolute voltage magnitude bounds for models with explicit neutrals
"""
function constraint_mc_voltage_absolute(pm::_PMD.RectangularVoltageExplicitNeutralModels, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    bus = _PMD.ref(pm, nw, :bus, id)

    constraint_mc_voltage_absolute(pm, nw, id, bus["terminals"], bus["grounded"], bus["vmin"], bus["vmax"])
end


"""
    function constraint_mc_voltage_pairwise(
        pm::_PMD.RectangularVoltageExplicitNeutralModels,
        id::Int;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true,
    )

Imposes pairwise voltage magnitude bounds, i.e. magnitude bounds on the voltage between to terminals, for models with explicit neutrals
"""
function constraint_mc_voltage_pairwise(pm::_PMD.RectangularVoltageExplicitNeutralModels, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    bus = _PMD.ref(pm, nw, :bus, id)

    vm_pair_lb = bus["vm_pair_lb"]
    vm_pair_ub = bus["vm_pair_ub"]

    constraint_mc_voltage_pairwise(pm, nw, id, bus["terminals"], bus["grounded"], vm_pair_lb, vm_pair_ub)
end


"""
    function constraint_mc_voltage_reference(
        pm::_PMD.ExplicitNeutralModels,
        id::Int;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true,
    )

Imposes suitable constraints for the voltage at the reference bus
"""
function constraint_mc_voltage_reference(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    bus = _PMD.ref(pm, nw, :bus, id)
    terminals = bus["terminals"]
    grounded = bus["grounded"]

    if haskey(bus, "va") && !haskey(bus, "vm")
        constraint_mc_theta_ref(pm, nw, id, bus["va"], terminals, grounded)
    elseif haskey(bus, "vm") && !haskey(bus, "va")
        constraint_mc_voltage_magnitude_fixed(pm, nw, id, bus["vm"], terminals, grounded)
    elseif haskey(bus, "vm") && haskey(bus, "va")
        constraint_mc_voltage_fixed(pm, nw, id, bus["vm"], bus["va"], terminals, grounded)
    end
end


"""
    function constraint_mc_generator_power(
        pm::_PMD.ExplicitNeutralModels,
        id::Int;
        nw::Int=nw_id_default,
        report::Bool=true
    )

Constrains generator power variables for models with explicit neutrals.
"""
function constraint_mc_generator_power(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, report::Bool=true)
    generator = _PMD.ref(pm, nw, :gen, id)
    bus = _PMD.ref(pm, nw,:bus, generator["gen_bus"])

    configuration = generator["configuration"]

    N = length(generator["connections"])
    pmin = get(generator, "pmin", fill(-Inf, N))
    pmax = get(generator, "pmax", fill( Inf, N))
    qmin = get(generator, "qmin", fill(-Inf, N))
    qmax = get(generator, "qmax", fill( Inf, N))

    if configuration==_PMD.WYE || length(pmin)==1
        constraint_mc_generator_power_wye(pm, nw, id, bus["index"], generator["connections"], pmin, pmax, qmin, qmax; report=report)
    else
        constraint_mc_generator_power_delta(pm, nw, id, bus["index"], generator["connections"], pmin, pmax, qmin, qmax; report=report)
    end
end


"""
	function constraint_mc_generator_current(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
		id::Int;
		nw::Int=nw_id_default,
		report::Bool=true,
		bounded::Bool=true
	)

For IVR models with explicit neutrals,
creates expressions for the terminal current flows `:crg_bus` and `:cig_bus`.
"""
function constraint_mc_generator_current(pm::_PMD.AbstractExplicitNeutralIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)
    generator = _PMD.ref(pm, nw, :gen, id)

    nphases = _PMD._infer_int_dim_unit(generator, false)
    # Note that one-dimensional delta generators are handled as wye-connected generators.
    # The distinction between one-dimensional wye and delta generators is purely semantic
    # when neutrals are modeled explicitly.
    if get(generator, "configuration", _PMD.WYE) == _PMD.WYE || nphases==1
        constraint_mc_generator_current_wye(pm, nw, id, generator["connections"]; report=report, bounded=bounded)
    else
        constraint_mc_generator_current_delta(pm, nw, id, generator["connections"]; report=report, bounded=bounded)
    end
end



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
    tm_scale = calculate_tm_scale(transformer, _PMD.ref(pm, nw, :bus, f_bus), _PMD.ref(pm, nw, :bus, t_bus))

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
    tm_scale = calculate_tm_scale(transformer, _PMD.ref(pm, nw, :bus, f_bus), _PMD.ref(pm, nw, :bus, t_bus))

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


# LOAD - Constraints

"""
	function constraint_mc_load_power(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
		id::Int;
		nw::Int=nw_id_default,
		report::Bool=true
	)

For IVR models with explicit neutrals,
the load power does not require any constraints.
"""
function constraint_mc_load_power(pm::_PMD.AbstractExplicitNeutralIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
    # nothing to do
end


"""
	function constraint_mc_load_power(
		pm::_PMD.ExplicitNeutralModels,
		id::Int;
		nw::Int=nw_id_default,
		report::Bool=true
	)

Constrains load power variables for models with explicit neutrals.
"""
function constraint_mc_load_power(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, report::Bool=true)
    load = _PMD.ref(pm, nw, :load, id)
    bus = _PMD.ref(pm, nw,:bus, load["load_bus"])

    configuration = load["configuration"]

    a, alpha, b, beta = _PMD._load_expmodel_params(load, bus)

    if configuration==_PMD.WYE || length(a)==1
        constraint_mc_load_power_wye(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    else
        constraint_mc_load_power_delta(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    end
end


"""
	function constraint_mc_load_current(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
		id::Int;
		nw::Int=nw_id_default,
		report::Bool=true
	)

For IVR models with explicit neutrals,
create non-linear expressions for the terminal current flows `:crd_bus` and `:cid_bus`
"""
function constraint_mc_load_current(pm::_PMD.AbstractExplicitNeutralIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
    load = _PMD.ref(pm, nw, :load, id)
    bus = _PMD.ref(pm, nw,:bus, load["load_bus"])

    configuration = load["configuration"]

    a, alpha, b, beta = _PMD._load_expmodel_params(load, bus)

    int_dim = _PMD._infer_int_dim_unit(load, false)
    if configuration==_PMD.WYE || int_dim==1
        constraint_mc_load_current_wye(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    else
        constraint_mc_load_current_delta(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    end
end



"""
    constraint_mc_current_from(pm::_PMD.AbstractUnbalancedIVRModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for current constraints on branches (from-side)
"""
function constraint_mc_current_from(pm::_PMD.AbstractUnbalancedIVRModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]

    constraint_mc_current_from(pm, nw, f_bus, f_idx, branch["f_connections"], g_fr, b_fr)
    nothing
end


"""
    constraint_mc_current_to(pm::_PMD.AbstractUnbalancedIVRModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for current constraints on branches (to-side)
"""
function constraint_mc_current_to(pm::_PMD.AbstractUnbalancedIVRModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g_to = branch["g_to"]
    b_to = branch["b_to"]

    constraint_mc_current_to(pm, nw, t_bus, f_idx, t_idx, branch["f_connections"], branch["t_connections"], g_to, b_to)
    nothing
end


"""
    constraint_mc_bus_voltage_drop(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for bus voltage drop constraints
"""
function constraint_mc_bus_voltage_drop(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    r = branch["br_r"]
    x = branch["br_x"]

    constraint_mc_bus_voltage_drop(pm, nw, i, f_bus, t_bus, f_idx, branch["f_connections"], branch["t_connections"], r, x)
    nothing
end


"""
	function constraint_mc_branch_current_limit(
		pm::_PMD.ExplicitNeutralModels,
		id::Int;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

For models with explicit neutrals,
imposes a bound on the current magnitude per conductor
at both ends of the branch (total current, i.e. including shunt contributions)
"""
function constraint_mc_branch_current_limit(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    branch = _PMD.ref(pm, nw, :branch, id)
    f_idx = (id,branch["f_bus"],branch["t_bus"])
    t_idx = (id,branch["t_bus"],branch["f_bus"])

    constraint_mc_branch_current_limit(pm, nw, f_idx, t_idx, branch["f_connections"], branch["t_connections"], branch["c_rating_a"])
end


"""
    constraint_mc_thermal_limit_from(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch thermal constraints (from-side)
"""
function constraint_mc_thermal_limit_from(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch, i)
    f_idx = (i, branch["f_bus"], branch["t_bus"])

    if !haskey(_PMD.con(pm, nw), :mu_sm_branch)
        _PMD.con(pm, nw)[:mu_sm_branch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    if haskey(branch, "rate_a") && any(branch["rate_a"] .< Inf)
        constraint_mc_thermal_limit_from(pm, nw, f_idx, branch["f_connections"], branch["rate_a"])
    end
    nothing
end


"""
    constraint_mc_thermal_limit_to(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch thermal constraints (to-side)
"""
function constraint_mc_thermal_limit_to(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    branch = _PMD.ref(pm, nw, :branch, i)
    t_idx = (i, branch["t_bus"], branch["f_bus"])

    if !haskey(_PMD.con(pm, nw), :mu_sm_branch)
        _PMD.con(pm, nw)[:mu_sm_branch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    if haskey(branch, "rate_a") && any(branch["rate_a"] .< Inf)
        constraint_mc_thermal_limit_to(pm, nw, t_idx, branch["t_connections"], branch["rate_a"])
    end
    nothing
end



"""
    constraint_mc_current_balance(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for KCL constraints in current-voltage variable space
"""
function constraint_mc_current_balance(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing
    bus = _PMD.ref(pm, nw, :bus, i)
    bus_arcs = _PMD.ref(pm, nw, :bus_arcs_conns_branch, i)
    bus_arcs_sw = _PMD.ref(pm, nw, :bus_arcs_conns_switch, i)
    bus_arcs_trans = _PMD.ref(pm, nw, :bus_arcs_conns_transformer, i)
    bus_gens = _PMD.ref(pm, nw, :bus_conns_gen, i)
    bus_storage = _PMD.ref(pm, nw, :bus_conns_storage, i)
    bus_loads = _PMD.ref(pm, nw, :bus_conns_load, i)
    bus_shunts = _PMD.ref(pm, nw, :bus_conns_shunt, i)

    constraint_mc_current_balance(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_shunts)
    nothing
end
