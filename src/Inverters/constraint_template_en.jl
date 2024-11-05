
# BRANCH - Constraints

"""
	function constraint_mc_branch_current_limit(
		pm::ExplicitNeutralModels,
		id::Int;
		nw::Int=_PMD.nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

For models with explicit neutrals,
imposes a bound on the current magnitude per conductor
at both ends of the branch (total current, i.e. including shunt contributions)
"""
function constraint_mc_branch_current_limit(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=_PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    branch = _PMD.ref(pm, nw, :branch, id)
    f_idx = (id,branch["f_bus"],branch["t_bus"])
    t_idx = (id,branch["t_bus"],branch["f_bus"])

    constraint_mc_branch_current_limit(pm, nw, f_idx, t_idx, branch["f_connections"], branch["t_connections"], branch["c_rating_a"])
end



"""
    constraint_mc_thermal_limit_from(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch thermal constraints (from-side)
"""
function constraint_mc_thermal_limit_from(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=_PMD.nw_id_default)::Nothing
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
    constraint_mc_thermal_limit_to(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for branch thermal constraints (to-side)
"""
function constraint_mc_thermal_limit_to(pm::_PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=_PMD.nw_id_default)::Nothing
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


# GENERATOR - Constraints

"""
    function constraint_mc_generator_power(
        pm::ExplicitNeutralModels,
        id::Int;
        nw::Int=_PMD.nw_id_default,
        report::Bool=true
    )

Constrains generator power variables for models with explicit neutrals.
"""
function constraint_mc_generator_power(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=_PMD.nw_id_default, report::Bool=true)
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
    function constraint_mc_generator_power(
        pm::ExplicitNeutralModels,
        id::Int;
        nw::Int=_PMD.nw_id_default,
        report::Bool=true
    )

Constrains generator power variables for models with explicit neutrals.
"""
function constraint_mc_inverter_dc_link_ripple_power(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=_PMD.nw_id_default, report::Bool=true)
    generator = _PMD.ref(pm, nw, :gen, id)
    bus = _PMD.ref(pm, nw,:bus, generator["gen_bus"])

    pdcmin = get(generator, "pdcmin", -Inf)
    pdcmax = get(generator, "pdcmax", Inf)
    
    constraint_mc_inverter_dc_link_ripple_power(pm, nw, id, bus["index"], generator["connections"], pdcmin, pdcmax; report=report)
end
