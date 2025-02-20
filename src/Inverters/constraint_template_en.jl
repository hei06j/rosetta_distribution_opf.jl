
# BRANCH - Constraints

"""
    constraint_mc_current_from(pm::AbstractUnbalancedIVRModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for current constraints on branches (from-side)
"""
function constraint_mc_current_from(pm::PMD.AbstractUnbalancedIVRModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    branch = PMD.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]

    constraint_mc_current_from(pm, nw, f_bus, f_idx, branch["f_connections"], g_fr, b_fr)
    nothing
end


"""
    constraint_mc_current_to(pm::AbstractUnbalancedIVRModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for current constraints on branches (to-side)
"""
function constraint_mc_current_to(pm::PMD.AbstractUnbalancedIVRModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    branch = PMD.ref(pm, nw, :branch, i)
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
	function constraint_mc_branch_current_limit(
		pm::ExplicitNeutralModels,
		id::Int;
		nw::Int=PMD.nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

For models with explicit neutrals,
imposes a bound on the current magnitude per conductor
at both ends of the branch (total current, i.e. including shunt contributions)
"""
function constraint_mc_branch_current_limit(pm::PMD.ExplicitNeutralModels, id::Int; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    branch = PMD.ref(pm, nw, :branch, id)
    f_idx = (id,branch["f_bus"],branch["t_bus"])
    t_idx = (id,branch["t_bus"],branch["f_bus"])

    constraint_mc_branch_current_limit(pm, nw, f_idx, t_idx, branch["f_connections"], branch["t_connections"], branch["c_rating_a"])
end



"""
    constraint_mc_thermal_limit_from(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch thermal constraints (from-side)
"""
function constraint_mc_thermal_limit_from(pm::PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    branch = PMD.ref(pm, nw, :branch, i)
    f_idx = (i, branch["f_bus"], branch["t_bus"])

    if !haskey(PMD.con(pm, nw), :mu_sm_branch)
        PMD.con(pm, nw)[:mu_sm_branch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    if haskey(branch, "rate_a") && any(branch["rate_a"] .< Inf)
        constraint_mc_thermal_limit_from(pm, nw, f_idx, branch["f_connections"], branch["rate_a"])
    end
    nothing
end



"""
    constraint_mc_thermal_limit_to(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch thermal constraints (to-side)
"""
function constraint_mc_thermal_limit_to(pm::PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    branch = PMD.ref(pm, nw, :branch, i)
    t_idx = (i, branch["t_bus"], branch["f_bus"])

    if !haskey(PMD.con(pm, nw), :mu_sm_branch)
        PMD.con(pm, nw)[:mu_sm_branch] = Dict{Tuple{Int,Int,Int}, Vector{JuMP.ConstraintRef}}()
    end

    if haskey(branch, "rate_a") && any(branch["rate_a"] .< Inf)
        constraint_mc_thermal_limit_to(pm, nw, t_idx, branch["t_connections"], branch["rate_a"])
    end
    nothing
end


# GENERATOR - Constraints

"""
    function constraint_mc_generator_power(
        pm::ExplicitNeutralModels,
        id::Int;
        nw::Int=PMD.nw_id_default,
        report::Bool=true
    )

Constrains generator power variables for models with explicit neutrals.
"""
function constraint_mc_generator_power(pm::PMD.ExplicitNeutralModels, id::Int; nw::Int=PMD.nw_id_default, report::Bool=true)
    generator = PMD.ref(pm, nw, :gen, id)
    bus = PMD.ref(pm, nw,:bus, generator["gen_bus"])

    configuration = generator["configuration"]

    N = length(generator["connections"])
    pmin = get(generator, "pmin", fill(-Inf, N))
    pmax = get(generator, "pmax", fill( Inf, N))
    qmin = get(generator, "qmin", fill(-Inf, N))
    qmax = get(generator, "qmax", fill( Inf, N))

    if configuration==PMD.WYE || length(pmin)==1
        constraint_mc_generator_power_wye(pm, nw, id, bus["index"], generator["connections"], pmin, pmax, qmin, qmax; report=report)
    else
        constraint_mc_generator_power_delta(pm, nw, id, bus["index"], generator["connections"], pmin, pmax, qmin, qmax; report=report)
    end
end


"""
	function constraint_mc_generator_current(
		pm::AbstractExplicitNeutralIVRModel,
		id::Int;
		nw::Int=nw_id_default,
		report::Bool=true,
		bounded::Bool=true
	)

For IVR models with explicit neutrals,
creates expressions for the terminal current flows `:crg_bus` and `:cig_bus`.
"""
function constraint_mc_generator_current(pm::PMD.AbstractExplicitNeutralIVRModel, id::Int; nw::Int=PMD.nw_id_default, report::Bool=true, bounded::Bool=true)
    generator = PMD.ref(pm, nw, :gen, id)

    nphases = PMD._infer_int_dim_unit(generator, false)
    # Note that one-dimensional delta generators are handled as wye-connected generators.
    # The distinction between one-dimensional wye and delta generators is purely semantic
    # when neutrals are modeled explicitly.
    if get(generator, "configuration", PMD.WYE) == PMD.WYE || nphases==1
        constraint_mc_generator_current_wye(pm, nw, id, generator["connections"]; report=report, bounded=bounded)
    else
        constraint_mc_generator_current_delta(pm, nw, id, generator["connections"]; report=report, bounded=bounded)
    end
end


function constraint_mc_generator_current_limit(pm::PMD.AbstractExplicitNeutralIVRModel, id::Int; nw::Int=PMD.nw_id_default, report::Bool=true, bounded::Bool=true)
    generator = PMD.ref(pm, nw, :gen, id)

    if haskey(generator, "c_rating") && any(generator["c_rating"] .< Inf)
        constraint_mc_generator_current_limit(pm, nw, id, generator["connections"], generator["c_rating"]; report=report, bounded=bounded)
    end
end


"""
	function constraint_mc_load_current(
		pm::AbstractExplicitNeutralIVRModel,
		id::Int;
		nw::Int=nw_id_default,
		report::Bool=true
	)

For IVR models with explicit neutrals,
create non-linear expressions for the terminal current flows `:crd_bus` and `:cid_bus`
"""
function constraint_mc_load_current(pm::PMD.AbstractExplicitNeutralIVRModel, id::Int; nw::Int=PMD.nw_id_default, report::Bool=true)
    load = PMD.ref(pm, nw, :load, id)
    bus = PMD.ref(pm, nw,:bus, load["load_bus"])

    configuration = load["configuration"]

    a, alpha, b, beta = PMD._load_expmodel_params(load, bus)

    int_dim = PMD._infer_int_dim_unit(load, false)
    if configuration==PMD.WYE || int_dim==1
        # constraint_mc_load_current_wye(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
        constraint_mc_load_current_wye(pm, nw, id, load["load_bus"], load["connections"], load["cm"], a, alpha, b, beta; report=report)
    else
        PMD.constraint_mc_load_current_delta(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    end
end



"""
    function constraint_mc_inverter_dc_link_ripple_power(
        pm::ExplicitNeutralModels,
        id::Int;
        nw::Int=PMD.nw_id_default,
        report::Bool=true
    )

Constrains generator power variables for models with explicit neutrals.
"""
function constraint_mc_inverter_dc_link_ripple_power(pm::PMD.ExplicitNeutralModels, id::Int; nw::Int=PMD.nw_id_default, report::Bool=true)
    generator = PMD.ref(pm, nw, :gen, id)
    bus = PMD.ref(pm, nw,:bus, generator["gen_bus"])

    pdcmin = get(generator, "pdcmin", -Inf)
    pdcmax = get(generator, "pdcmax", Inf)
    
    constraint_mc_inverter_dc_link_ripple_power(pm, nw, id, bus["index"], generator["connections"], pdcmin, pdcmax; report=report)
end


"""
    function constraint_mc_inverter_branch_dc_link_ripple_power(
        pm::ExplicitNeutralModels,
        id::Int;
        nw::Int=PMD.nw_id_default,
        report::Bool=true
    )

Constrains generator power variables for models with explicit neutrals.
"""
function constraint_mc_inverter_branch_dc_link_ripple_power(pm::PMD.ExplicitNeutralModels, id::Int; nw::Int=PMD.nw_id_default, report::Bool=true)
    branch = PMD.ref(pm, nw, :branch, id)
    f_idx = (id, branch["f_bus"], branch["t_bus"])

    pdcmin = get(branch, "pdcmin", -Inf)
    pdcmax = get(branch, "pdcmax", Inf)
    
    constraint_mc_inverter_branch_dc_link_ripple_power(pm, nw, id, f_idx, branch["f_connections"], pdcmin, pdcmax)
end
