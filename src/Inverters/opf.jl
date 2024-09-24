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

function variable_multiplexing_inverter(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    pv_gen_ids = [i for (i, gen) in pm.data["gen"] if !occursin("source", gen["name"])]
    connections = Dict(i => pm.data["gen"]["$i"]["connections"] for i in pv_gen_ids)
    m_legs = Dict(i => pm.data["gen"]["$i"]["m_legs"] for i in pv_gen_ids)
    
    bg = _PMD.var(pm, nw)[:bg] = Dict(parse(Int, i) => JuMP.@variable(pm.model, [connections[i], 1:m_legs[i]], base_name="bg_$i", Bin) for i in pv_gen_ids)
    # report && sol_component_value(pm, _PMD.pmd_it_sym, nw, :gen, :bg, parse.(Int, pv_gen_ids)[1], bg)
    for i in 1:length(pv_gen_ids)
        report && sol_component_value_comp(pm, _PMD.pmd_it_sym, nw, :gen, :bg, parse.(Int, i), bg[parse.(Int, i)])
    end
end



"""
	function constraint_mc_branch_current_limit(
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
function constraint_mc_branch_current_limit(pm::_PMD.ExplicitNeutralModels, id::Int, gen_id::Int; nw::Int=_PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    branch = _PMD.ref(pm, nw, :branch, id)
    f_idx = (id,branch["f_bus"],branch["t_bus"])
    t_idx = (id,branch["t_bus"],branch["f_bus"])
    constraint_mc_branch_current_limit(pm, nw, f_idx, t_idx, branch["c_rating_a"], gen_id)
end

"""
	function constraint_mc_branch_current_limit(
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
function constraint_mc_branch_current_limit(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, c_rating::Vector{<:Real}, gen_id; report::Bool=true)
    cr_fr = _PMD.var(pm, nw, :cr, f_idx)
    ci_fr = _PMD.var(pm, nw, :ci, f_idx)
    cr_to = _PMD.var(pm, nw, :cr, t_idx)
    ci_to = _PMD.var(pm, nw, :ci, t_idx)
    bg = _PMD.var(pm, nw, :bg, gen_id)

    m_legs = pm.data["gen"]["$gen_id"]["m_legs"]
    alpha_g = 1/m_legs * ones(m_legs)

    JuMP.@constraint(pm.model, [k in 1:size(bg,2)], sum(bg[:,k]) == 1)

    c_rating = JuMP.@expression(pm.model,  sum(c_rating) * Array(bg) * alpha_g)

    cnds_finite_rating = [c for (c,r) in enumerate(c_rating) if r!==Inf]
    JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_fr[c]^2+ci_fr[c]^2 <= c_rating[c]^2)
    JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_to[c]^2+ci_to[c]^2 <= c_rating[c]^2)
end


function objective_mc_min_IUF(pm::_PMD.AbstractUnbalancedPowerModel)
    alpha = exp(im*2/3*pi)
    T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
    Tre = real.(T)
    Tim = imag.(T)

    ref = pm.ref[:it][:pmd][:nw][0]   # TODO change 0 to nw, make this ref dependent
    _, _, arc, branch = get_ref_bus_branch(ref)
    nconds = Dict(l => length(ref[:branch][l]["f_connections"]) for l in [branch])
    conds = Dict(l => ref[:branch][l]["f_connections"] for l in [branch])
    
    n_ph = 4  # TODO, make this dependent on nconds, and set to zero the ones that are not equal to nconds[i]

    ### TODO these variables and constraints should be defiened outside of here, but only added to the model whithin this objective function, 
    ### so that there are only added to the model if this objective is called
    cr_012 = Dict((l,i,j) => JuMP.@variable(pm.model, [c in conds[l]], base_name="cr_012_$((l,i,j))") for (l,i,j) in [arc])
    cr_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in cr_012[(l,i,j)].axes[1] ? cr_012[(l,i,j)][c] : 0.0 for c in 1:n_ph, (l,i,j) in [arc]]), 1:n_ph, [arc])
    ci_012 = Dict((l,i,j) => JuMP.@variable(pm.model, [c in conds[l]], base_name="ci_012_$((l,i,j))") for (l,i,j) in [arc])
    ci_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in ci_012[(l,i,j)].axes[1] ? ci_012[(l,i,j)][c] : 0.0 for c in 1:n_ph, (l,i,j) in [arc]]), 1:n_ph, [arc])
    cm_012 = Dict((l,i,j) => JuMP.@variable(pm.model, [c in conds[l]], base_name="cm_012_$((l,i,j))", lower_bound=0) for (l,i,j) in [arc])
    cm_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in cm_012[(l,i,j)].axes[1] ? cm_012[(l,i,j)][c] : 0.0 for c in 1:n_ph, (l,i,j) in [arc]]), 1:n_ph, [arc])

    cr_bus = _PMD.var(pm, 0, :cr_bus, arc)  # TODO change 0 to nw
    ci_bus = _PMD.var(pm, 0, :ci_bus, arc)  # TODO change 0 to nw

    phases = 1:3
    JuMP.@constraint(pm.model, cr_012[phases,:] .== Tre * Array(cr_bus[phases]) .- Tim * Array(ci_bus[phases]))
    JuMP.@constraint(pm.model, ci_012[phases,:] .== Tre * Array(ci_bus[phases]) .+ Tim * Array(cr_bus[phases]))
    JuMP.@constraint(pm.model, cm_012[phases,:].^2 .== cr_012[phases,:].^2 .+ ci_012[phases,:].^2)

    # if objective == "IUF_inv"
    #     # JuMP.@objective(model, Min, cmg_012[3,1] / cmg_012[2,1])
    #     JuMP.@objective(model, Min, cm_012[3,arc] / cm_012[2,arc])

    # elseif objective == "IUF2_inv"
        # JuMP.@objective(model, Min, cmg_012[3,1])
        JuMP.@objective(pm.model, Min, cm_012[3,arc])
    # end

end



"""
	function build_mc_opf(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for OPF in current-voltage variable space with explicit neutrals
"""
function build_mc_opf(pm::_PMD.AbstractExplicitNeutralIVRModel)
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

    if pm.setting["multileg"] && pm.setting["multiplexing"]
        variable_multiplexing_inverter(pm)
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
        _PMD.constraint_mc_generator_power(pm, id)
        _PMD.constraint_mc_generator_current(pm, id)
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

    inverter_branches = [branch["index"] for (i, branch) in pm.data["branch"] if occursin("inverter_branch", branch["name"])]
    gen_ids, gen_bus_ids, gen_branch_ids = get_pv_bus_branch(pm.ref[:it][:pmd][:nw][0])  # maybe output branch_ids and arcs seperately?
    branch_gens = Dict(branch_id[1] => gen_ids[i] for (i, branch_id) in enumerate(gen_branch_ids))

    for i in _PMD.ids(pm, :branch)
        _PMD.constraint_mc_current_from(pm, i)
        _PMD.constraint_mc_current_to(pm, i)
        _PMD.constraint_mc_bus_voltage_drop(pm, i)

        if (i ∈ inverter_branches) && pm.setting["multileg"] && pm.setting["multiplexing"]
            constraint_mc_branch_current_limit(pm, i, branch_gens[i])
        else
            _PMD.constraint_mc_branch_current_limit(pm, i)
        end

        _PMD.constraint_mc_thermal_limit_from(pm, i)
        _PMD.constraint_mc_thermal_limit_to(pm, i)
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
    # _PMD.objective_mc_min_fuel_cost(pm)
    objective_mc_min_IUF(pm)
end
