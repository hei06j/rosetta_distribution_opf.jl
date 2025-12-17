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
function constraint_mc_branch_current_limit_mx(pm::PMD.ExplicitNeutralModels, id::Int, gen_id::Int; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    branch = PMD.ref(pm, nw, :branch, id)
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
function constraint_mc_branch_current_limit_mx(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, c_rating::Vector{<:Real}, gen_id; report::Bool=true)
    cr_fr = PMD.var(pm, nw, :cr, f_idx)
    ci_fr = PMD.var(pm, nw, :ci, f_idx)
    cr_to = PMD.var(pm, nw, :cr, t_idx)
    ci_to = PMD.var(pm, nw, :ci, t_idx)
    bg = PMD.var(pm, nw, :bg, gen_id)

    m_legs = pm.data["gen"]["$gen_id"]["m_legs"]
    alpha_g = 1/m_legs * ones(m_legs)

    JuMP.@constraint(pm.model, [k in 1:size(bg,2)], sum(bg[:,k]) == 1)

    c_rating = JuMP.@expression(pm.model,  sum(c_rating) * Array(bg) * alpha_g)
    PMD.var(pm, nw, :c_rating)[f_idx[1]] = c_rating

    if report
        PMD.sol(pm, nw, :branch, f_idx[1])[:c_rating] = c_rating
    end

    cnds_finite_rating = [c for (c,r) in enumerate(c_rating) if r!==Inf]
    JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_fr[c]^2+ci_fr[c]^2 <= c_rating[c]^2)
    JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_to[c]^2+ci_to[c]^2 <= c_rating[c]^2)
end


function constraint_inverter_branch_balance(pm::PMD.ExplicitNeutralModels, id::Int; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    branch = PMD.ref(pm, nw, :branch, id)
    f_idx = (id, branch["f_bus"], branch["t_bus"])
    t_idx = (id, branch["t_bus"], branch["f_bus"])
    # r = 0.015 / (230.94^2 / 1000)
    # x = 0
    # c_rating = 3 * 0.0033  # minimum(pmax1, pmax2)
    # m_legs = branch["m_legs"]
    constraint_inverter_branch_balance(pm, nw, id, f_idx, t_idx, branch["f_connections"], branch["t_connections"], branch["c_rating_a"])
end

function constraint_inverter_branch_balance(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, branch_id, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections, t_connections, c_rating_a::Vector{<:Real}; report::Bool=true)
    cr_fr = PMD.var(pm, nw, :cr, f_idx)
    ci_fr = PMD.var(pm, nw, :ci, f_idx)
    cr_to = PMD.var(pm, nw, :cr, t_idx)
    ci_to = PMD.var(pm, nw, :ci, t_idx)

    csr_fr = PMD.var(pm, nw, :csr, f_idx[1])
    csi_fr = PMD.var(pm, nw, :csi, f_idx[1])
    csr_to = PMD.var(pm, nw, :csr, t_idx[1])
    csi_to = PMD.var(pm, nw, :csi, t_idx[1])

    vr_fr = [PMD.var(pm, nw, :vr, f_idx[2])[t] for t in f_connections]
    vi_fr = [PMD.var(pm, nw, :vi, f_idx[2])[t] for t in f_connections]
    vr_to = [PMD.var(pm, nw, :vr, t_idx[2])[t] for t in t_connections]
    vi_to = [PMD.var(pm, nw, :vi, t_idx[2])[t] for t in t_connections]
    
    ### constraint sum(csr_fr) + sum(csr_to) = 0,   sum(csi_fr) + sum(csi_to) = 0
    JuMP.@constraint(pm.model, sum(cr_fr) == 0)
    JuMP.@constraint(pm.model, sum(ci_fr) == 0)
    # JuMP.@constraint(pm.model, sum(csr_fr) == 0)
    # JuMP.@constraint(pm.model, sum(csi_fr) == 0)
    # JuMP.@constraint(pm.model, sum(cr_fr) + sum(cr_to) == 0)
    # JuMP.@constraint(pm.model, sum(ci_fr) + sum(ci_to) == 0)
    # JuMP.@constraint(pm.model, sum(csr_fr) + sum(csr_to) == 0)
    # JuMP.@constraint(pm.model, sum(csi_fr) + sum(csi_to) == 0)



    PMD.var(pm, 0)[:pf_idx] = Dict{Int, Any}()
    PMD.var(pm, 0)[:pt_idx] = Dict{Int, Any}()

    pf_idx = JuMP.@expression(pm.model,  vr_fr .* cr_fr .+ vi_fr .* ci_fr)
    pt_idx = JuMP.@expression(pm.model,  vr_to .* cr_to .+ vi_to .* ci_to)
    JuMP.@constraint(pm.model, sum(pf_idx) == 0)
    # JuMP.@constraint(pm.model, sum(pf_idx) + sum(pt_idx) == 0)
    
    PMD.var(pm, nw, :pf_idx)[branch_id] = pf_idx
    PMD.var(pm, nw, :pt_idx)[branch_id] = pt_idx
    if report
        PMD.sol(pm, nw, :branch, branch_id)[:pf_idx] = pf_idx
        PMD.sol(pm, nw, :branch, branch_id)[:pt_idx] = pt_idx
    end



    ### constraint_mc_thermal_limit
    # if haskey(branch, "rate_a") && any(branch["rate_a"] .< Inf)
    #     ### constraint_mc_thermal_limit_from
    #     # pf_idx = JuMP.@expression(model,  vr_fr .* cr_fr .+ vi_fr .* ci_fr)
    #     qf_idx = JuMP.@expression(model, -vr_fr .* ci_fr .+ vi_fr .* cr_fr)
    #     JuMP.@constraint(model, pf_idx.^2 .+ qf_idx.^2 .<= branch["rate_a"].^2)

    #     ### constraint_mc_thermal_limit_to
    #     # pt_idx = JuMP.@expression(model,  vr_to .* cr_to .+ vi_to .* ci_to)
    #     qt_idx = JuMP.@expression(model, -vr_to .* ci_to .+ vi_to .* cr_to)
    #     JuMP.@constraint(model, pt_idx.^2 .+ qt_idx.^2 .<= branch["rate_a"].^2)
    # end

end


"""
	function build_mc_opf_mx(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for OPF in current-voltage variable space with explicit neutrals including reconfigurable inverters
"""
function build_mc_opf_mx(pm::PMD.AbstractExplicitNeutralIVRModel)
    inverter_branches = [branch["index"] for (i, branch) in pm.data["branch"] if occursin("inverter_branch", branch["name"])]
    gen_ids, gen_bus_ids, gen_branch_ids = get_pv_bus_branch(pm.ref[:it][:pmd][:nw][0])  # maybe output branch_ids and arcs seperately?
    branch_gens = Dict(branch_id[1] => gen_ids[i] for (i, branch_id) in enumerate(gen_branch_ids))


    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_load_current(pm)
    PMD.variable_mc_load_power(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_generator_power(pm)
    # variable_mc_generator_current(pm)
    # variable_mc_generator_power(pm)
    PMD.variable_mc_transformer_current(pm)
    PMD.variable_mc_transformer_power(pm)
    PMD.variable_mc_switch_current(pm)

    if pm.setting["reconfigurable"]
        variable_reconfigurable_inverter(pm)
        PMD.var(pm, 0)[:c_rating] = Dict{Int, Any}()
    end

    if pm.setting["dc_link"]
        PMD.var(pm, 0)[:pdc_link_sqr] = Dict{Int, Any}()
    end

    # Constraints
    for i in PMD.ids(pm, :bus)

        if i in PMD.ids(pm, :ref_buses)
            PMD.constraint_mc_voltage_reference(pm, i)
        end

        PMD.constraint_mc_voltage_absolute(pm, i)
        PMD.constraint_mc_voltage_pairwise(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in PMD.ids(pm, :gen)
        if id ∈ gen_ids  # Generators connected with inverter
            constraint_mc_generator_power(pm, id)
            # PMD.constraint_mc_generator_current(pm, id)
            constraint_mc_generator_current(pm, id)
            constraint_mc_generator_current_limit(pm, id)
            
            if pm.setting["dc_link"]
                constraint_mc_inverter_dc_link_ripple_power(pm, id)
            end

        else  # Other generators
            PMD.constraint_mc_generator_power(pm, id)
            PMD.constraint_mc_generator_current(pm, id)
        end
    end

    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
        PMD.constraint_mc_load_current(pm, id)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_voltage(pm, i)
        PMD.constraint_mc_transformer_current(pm, i)

        PMD.constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :branch)

        if i ∈ inverter_branches

            if  pm.setting["reconfigurable"]
                PMD.constraint_mc_current_from(pm, i)
                PMD.constraint_mc_current_to(pm, i)
                # constraint_mc_current_from(pm, i)
                # constraint_mc_current_to(pm, i)
                PMD.constraint_mc_bus_voltage_drop(pm, i)
                constraint_mc_branch_current_limit_mx(pm, i, branch_gens[i])
                constraint_inverter_branch_balance(pm, i)
                # PMD.constraint_mc_thermal_limit_from(pm, i)
                # PMD.constraint_mc_thermal_limit_to(pm, i)

            elseif pm.setting["ideal"]
                PMD.constraint_mc_current_from(pm, i)
                PMD.constraint_mc_current_to(pm, i)
                # constraint_mc_current_from(pm, i)
                # constraint_mc_current_to(pm, i)
                PMD.constraint_mc_bus_voltage_drop(pm, i)
                constraint_mc_branch_current_limit(pm, i)
                constraint_inverter_branch_balance(pm, i)
                # constraint_mc_thermal_limit_from(pm, i)
                # constraint_mc_thermal_limit_to(pm, i)

            elseif pm.setting["conventional"]
                PMD.constraint_mc_current_from(pm, i)
                PMD.constraint_mc_current_to(pm, i)
                # constraint_mc_current_from(pm, i)
                # constraint_mc_current_to(pm, i)
                PMD.constraint_mc_bus_voltage_drop(pm, i)
                PMD.constraint_mc_branch_current_limit(pm, i)
                # constraint_inverter_branch_balance(pm, i)
                # PMD.constraint_mc_thermal_limit_from(pm, i)
                # PMD.constraint_mc_thermal_limit_to(pm, i)
            end

        else  # normal branch
            PMD.constraint_mc_current_from(pm, i)
            PMD.constraint_mc_current_to(pm, i)
            PMD.constraint_mc_bus_voltage_drop(pm, i)
            PMD.constraint_mc_branch_current_limit(pm, i)
            # PMD.constraint_mc_thermal_limit_from(pm, i)
            # PMD.constraint_mc_thermal_limit_to(pm, i)
        end

    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_current(pm, i)
        PMD.constraint_mc_switch_state(pm, i)

        PMD.constraint_mc_switch_current_limit(pm, i)
        PMD.constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
    end

    # Objective
    # PMD.objective_mc_min_fuel_cost(pm)
    # objective_mc_min_IUF(pm)
    objective_mc_min_max_phase_current(pm)
    # objective_mc_min_ref_branch_loss(pm)
end



"""
	function build_mc_opf_mx_cost(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for OPF in current-voltage variable space with explicit neutrals including reconfigurable inverters
"""
function build_mc_opf_mx_cost(pm::PMD.AbstractExplicitNeutralIVRModel)
    inverter_branches = [branch["index"] for (i, branch) in pm.data["branch"] if occursin("inverter_branch", branch["name"])]
    gen_ids, gen_bus_ids, gen_branch_ids = get_pv_bus_branch(pm.ref[:it][:pmd][:nw][0])  # maybe output branch_ids and arcs seperately?
    branch_gens = Dict(branch_id[1] => gen_ids[i] for (i, branch_id) in enumerate(gen_branch_ids))


    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_load_current(pm)
    PMD.variable_mc_load_power(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_generator_power(pm)
    # variable_mc_generator_current(pm)
    # variable_mc_generator_power(pm)
    PMD.variable_mc_transformer_current(pm)
    PMD.variable_mc_transformer_power(pm)
    PMD.variable_mc_switch_current(pm)

    if pm.setting["reconfigurable"]
        variable_reconfigurable_inverter(pm)
        PMD.var(pm, 0)[:c_rating] = Dict{Int, Any}()
    end

    if pm.setting["dc_link"]
        PMD.var(pm, 0)[:pdc_link_sqr] = Dict{Int, Any}()
    end

    # Constraints
    for i in PMD.ids(pm, :bus)

        if i in PMD.ids(pm, :ref_buses)
            PMD.constraint_mc_voltage_reference(pm, i)
        end

        PMD.constraint_mc_voltage_absolute(pm, i)
        PMD.constraint_mc_voltage_pairwise(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in PMD.ids(pm, :gen)
        if id ∈ gen_ids  # Generators connected with inverter
            constraint_mc_generator_power(pm, id)
            PMD.constraint_mc_generator_current(pm, id)
            # constraint_mc_generator_current(pm, id)
            constraint_mc_generator_current_limit(pm, id)
            
            if pm.setting["dc_link"]
                constraint_mc_inverter_dc_link_ripple_power(pm, id)
            end

        else  # Other generators
            PMD.constraint_mc_generator_power(pm, id)
            PMD.constraint_mc_generator_current(pm, id)
        end
    end

    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
        PMD.constraint_mc_load_current(pm, id)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_voltage(pm, i)
        PMD.constraint_mc_transformer_current(pm, i)

        PMD.constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :branch)

        if i ∈ inverter_branches

            if  pm.setting["reconfigurable"]
                PMD.constraint_mc_current_from(pm, i)
                PMD.constraint_mc_current_to(pm, i)
                # constraint_mc_current_from(pm, i)
                # constraint_mc_current_to(pm, i)
                PMD.constraint_mc_bus_voltage_drop(pm, i)
                constraint_mc_branch_current_limit_mx(pm, i, branch_gens[i])
                constraint_inverter_branch_balance(pm, i)
                # PMD.constraint_mc_thermal_limit_from(pm, i)
                # PMD.constraint_mc_thermal_limit_to(pm, i)

            elseif pm.setting["ideal"]
                PMD.constraint_mc_current_from(pm, i)
                PMD.constraint_mc_current_to(pm, i)
                # constraint_mc_current_from(pm, i)
                # constraint_mc_current_to(pm, i)
                PMD.constraint_mc_bus_voltage_drop(pm, i)
                constraint_mc_branch_current_limit(pm, i)
                constraint_inverter_branch_balance(pm, i)
                # constraint_mc_thermal_limit_from(pm, i)
                # constraint_mc_thermal_limit_to(pm, i)

            elseif pm.setting["conventional"]
                PMD.constraint_mc_current_from(pm, i)
                PMD.constraint_mc_current_to(pm, i)
                # constraint_mc_current_from(pm, i)
                # constraint_mc_current_to(pm, i)
                PMD.constraint_mc_bus_voltage_drop(pm, i)
                PMD.constraint_mc_branch_current_limit(pm, i)
                # constraint_inverter_branch_balance(pm, i)
                # PMD.constraint_mc_thermal_limit_from(pm, i)
                # PMD.constraint_mc_thermal_limit_to(pm, i)
            end

        else  # normal branch
            PMD.constraint_mc_current_from(pm, i)
            PMD.constraint_mc_current_to(pm, i)
            PMD.constraint_mc_bus_voltage_drop(pm, i)
            PMD.constraint_mc_branch_current_limit(pm, i)
            # PMD.constraint_mc_thermal_limit_from(pm, i)
            # PMD.constraint_mc_thermal_limit_to(pm, i)
        end

    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_current(pm, i)
        PMD.constraint_mc_switch_state(pm, i)

        PMD.constraint_mc_switch_current_limit(pm, i)
        PMD.constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
    end

    # Objective
    PMD.objective_mc_min_fuel_cost(pm)
end