"given a variable that is indexed by component ids, builds the standard solution structure"
function sol_component_value_comp(aim::_IM.AbstractInfrastructureModel, it::Symbol, n::Int, comp_name::Symbol, field_name::Symbol, comp_id, variable)
    @assert !haskey(_IM.sol(aim, it, n, comp_name, comp_id), field_name)
    _IM.sol(aim, it, n, comp_name, comp_id)[field_name] = variable
end

"given a variable that is indexed by component ids, builds the standard solution structure"
function sol_component_value(aim::_IM.AbstractInfrastructureModel, it::Symbol, n::Int, comp_name::Symbol, field_name::Symbol, comp_ids, variables)
    for i in comp_ids
        @assert !haskey(_IM.sol(aim, it, n, comp_name, i), field_name)
        _IM.sol(aim, it, n, comp_name, i)[field_name] = variables[i]
    end
    return 1
end


function variable_reconfigurable_inverter(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    pv_gen_ids = [i for (i, gen) in pm.data["gen"] if !occursin("source", gen["name"])]
    # connections = Dict(i => pm.data["gen"]["$i"]["connections"] for i in pv_gen_ids)
    connections = Dict(i => length(pm.data["gen"]["$i"]["connections"]) for i in pv_gen_ids)
    m_legs = Dict(i => pm.data["gen"]["$i"]["m_legs"] for i in pv_gen_ids)

    # bg = _PMD.var(pm, nw)[:bg] = Dict(parse(Int, i) => JuMP.@variable(pm.model, [connections[i], 1:m_legs[i]], base_name="bg_$i", Bin) for i in pv_gen_ids)
    ### report && sol_component_value(pm, _PMD.pmd_it_sym, nw, :gen, :bg, parse.(Int, pv_gen_ids)[1], bg)
    # for i in 1:length(pv_gen_ids)
    #     report && sol_component_value_comp(pm, _PMD.pmd_it_sym, nw, :gen, :bg, parse.(Int, i), bg[parse.(Int, i)])
    # end
    bg = _PMD.var(pm, nw)[:bg] = Dict(parse(Int, i) => JuMP.@variable(pm.model, [1:connections[i], 1:m_legs[i]], base_name="bg_$i", Bin) for i in pv_gen_ids)
    for i in pv_gen_ids
        report && sol_component_value_comp(pm, _PMD.pmd_it_sym, nw, :gen, :bg, parse.(Int, i), bg[parse.(Int, i)])
    end
end



"""
	function constraint_mc_branch_current_limit_mx(
		pm::ExplicitNeutralModels,
		id::Int;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

For models with explicit neutrals,
imposes a bound on the current magnitude per conductor
at both ends of the branch (total current, i.e. including shunt contributions)
"""
function constraint_mc_branch_current_limit_mx(pm::_PMD.ExplicitNeutralModels, id::Int, gen_id::Int; nw::Int=_PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    branch = _PMD.ref(pm, nw, :branch, id)
    f_idx = (id,branch["f_bus"],branch["t_bus"])
    t_idx = (id,branch["t_bus"],branch["f_bus"])
    constraint_mc_branch_current_limit_mx(pm, nw, f_idx, t_idx, branch["c_rating_a"], gen_id)
end

"""
	function constraint_mc_branch_current_limit_mx(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector,
		t_connections::Vector,
		c_rating::Vector{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
imposes a bound on the current magnitude per conductor
at both ends of the branch (total current, i.e. including shunt contributions).

```
cr_fr^2 + ci_fr^2 <= c_rating^2
cr_to^2 + ci_to^2 <= c_rating^2
```
"""
function constraint_mc_branch_current_limit_mx(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, c_rating::Vector{<:Real}, gen_id; report::Bool=true)
    cr_fr = _PMD.var(pm, nw, :cr, f_idx)
    ci_fr = _PMD.var(pm, nw, :ci, f_idx)
    cr_to = _PMD.var(pm, nw, :cr, t_idx)
    ci_to = _PMD.var(pm, nw, :ci, t_idx)
    bg = _PMD.var(pm, nw, :bg, gen_id)

    m_legs = pm.data["gen"]["$gen_id"]["m_legs"]
    alpha_g = 1/m_legs * ones(m_legs)

    JuMP.@constraint(pm.model, [k in 1:size(bg,2)], sum(bg[:,k]) == 1)

    c_rating = JuMP.@expression(pm.model,  sum(c_rating) * Array(bg) * alpha_g)
    _PMD.var(pm, nw, :c_rating)[f_idx[1]] = c_rating

    if report
        _PMD.sol(pm, nw, :branch, f_idx[1])[:c_rating] = c_rating
    end

    cnds_finite_rating = [c for (c,r) in enumerate(c_rating) if r!==Inf]
    JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_fr[c]^2+ci_fr[c]^2 <= c_rating[c]^2)
    JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_to[c]^2+ci_to[c]^2 <= c_rating[c]^2)
end


"""
	function build_mc_opf_mx(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for OPF in current-voltage variable space with explicit neutrals including reconfigurable inverters
"""
function build_mc_opf_mx(pm::_PMD.AbstractExplicitNeutralIVRModel)
    # Variables
    _PMD.variable_mc_bus_voltage(pm)
    _PMD.variable_mc_branch_current(pm)
    _PMD.variable_mc_load_current(pm)
    _PMD.variable_mc_load_power(pm)
    _PMD.variable_mc_generator_current(pm)
    _PMD.variable_mc_generator_power(pm)
    _PMD.variable_mc_transformer_current(pm)
    _PMD.variable_mc_transformer_power(pm)
    _PMD.variable_mc_switch_current(pm)

    inverter_branches = [branch["index"] for (i, branch) in pm.data["branch"] if occursin("inverter_branch", branch["name"])]
    gen_ids, gen_bus_ids, gen_branch_ids = get_pv_bus_branch(pm.ref[:it][:pmd][:nw][0])  # maybe output branch_ids and arcs seperately?
    branch_gens = Dict(branch_id[1] => gen_ids[i] for (i, branch_id) in enumerate(gen_branch_ids))

    if pm.setting["reconfigurable"]
        variable_reconfigurable_inverter(pm)
        _PMD.var(pm, 0)[:c_rating] = Dict{Int, Any}()
    end

    if pm.setting["dc_link"]
        _PMD.var(pm, 0)[:pdc_link] = Dict{Int, Any}()
    end

    # Constraints
    for i in _PMD.ids(pm, :bus)

        if i in _PMD.ids(pm, :ref_buses)
            _PMD.constraint_mc_voltage_reference(pm, i)
        end

        _PMD.constraint_mc_voltage_absolute(pm, i)
        _PMD.constraint_mc_voltage_pairwise(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in _PMD.ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
        # _PMD.constraint_mc_generator_power(pm, id)
        _PMD.constraint_mc_generator_current(pm, id)

        if id ∈ gen_ids && pm.setting["dc_link"]
            constraint_mc_inverter_dc_link_ripple_power(pm, id)
        end
    end

    for id in _PMD.ids(pm, :load)
        _PMD.constraint_mc_load_power(pm, id)
        _PMD.constraint_mc_load_current(pm, id)
    end

    for i in _PMD.ids(pm, :transformer)
        _PMD.constraint_mc_transformer_voltage(pm, i)
        _PMD.constraint_mc_transformer_current(pm, i)

        _PMD.constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in _PMD.ids(pm, :branch)
        _PMD.constraint_mc_current_from(pm, i)
        _PMD.constraint_mc_current_to(pm, i)
        _PMD.constraint_mc_bus_voltage_drop(pm, i)

        if (i ∈ inverter_branches) && pm.setting["reconfigurable"]
            constraint_mc_branch_current_limit_mx(pm, i, branch_gens[i])

        elseif (i ∈ inverter_branches) && pm.setting["ideal"]
            constraint_mc_branch_current_limit(pm, i)

        else # normal branch, or conventional inverter
            _PMD.constraint_mc_branch_current_limit(pm, i)
        end

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in _PMD.ids(pm, :switch)
        _PMD.constraint_mc_switch_current(pm, i)
        _PMD.constraint_mc_switch_state(pm, i)

        _PMD.constraint_mc_switch_current_limit(pm, i)
        _PMD.constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in _PMD.ids(pm, :bus)
        _PMD.constraint_mc_current_balance(pm, i)
    end

    # Objective
    _PMD.objective_mc_min_fuel_cost(pm)
    # objective_mc_min_IUF(pm)
end