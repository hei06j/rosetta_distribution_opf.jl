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

function variable_multiplexing_sop_inverter(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    sop_branches = [i for (i, branch) in pm.data["branch"] if occursin("SOP_branch", branch["name"])]
    # f_connections = Dict(i => pm.data["branch"]["$i"]["f_connections"] for i in sop_branches)
    # t_connections = Dict(i => pm.data["branch"]["$i"]["f_connections"] for i in sop_branches)
    connections = Dict(i => length(pm.data["branch"]["$i"]["f_connections"]) + length(pm.data["branch"]["$i"]["t_connections"]) for i in sop_branches)
    m_legs = Dict(i => pm.data["branch"]["$i"]["m_legs"] for i in sop_branches) 
    # TODO add assert to make sure that m_legs >= connections for each sop, unless the neutral is not linked to legs and directly connected to each other
    
    # bg = _PMD.var(pm, nw)[:bg] = Dict(parse(Int, i) => JuMP.@variable(pm.model, [connections[i], 1:m_legs[i]], base_name="bg_$i", Bin) for i in sop_branches)
    bg = _PMD.var(pm, nw)[:bg] = Dict(parse(Int, i) => JuMP.@variable(pm.model, [1:connections[i], 1:m_legs[i]], base_name="bg_$i", Bin) for i in sop_branches)
    for i in sop_branches
        report && sol_component_value_comp(pm, _PMD.pmd_it_sym, nw, :branch, :bg, parse.(Int, i), bg[parse.(Int, i)])
    end

    # bg = _PMD.var(pm, nw)[:bg] = Dict(parse(Int, i) => JuMP.@variable(pm.model, [connections[i], 1:m_legs[i]], base_name="bg_$i", Bin) for i in sop_branches)
    cr_max = _PMD.var(pm, nw)[:cr_max] = Dict(parse(Int, i) => JuMP.@variable(pm.model, [1:connections[i]], base_name="cr_max_$i", lower_bound=0) for i in sop_branches)
    for i in sop_branches
        report && sol_component_value_comp(pm, _PMD.pmd_it_sym, nw, :branch, :cr_max, parse.(Int, i), cr_max[parse.(Int, i)])
    end


    
end


function constraint_sop_branch(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=_PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    branch = _PMD.ref(pm, nw, :branch, id)
    f_idx = (id, branch["f_bus"], branch["t_bus"])
    t_idx = (id, branch["t_bus"], branch["f_bus"])
    # r = 0.015 / (230.94^2 / 1000)
    # x = 0
    # c_rating = 3 * 0.0033  # minimum(pmax1, pmax2)
    # m_legs = branch["m_legs"]
    constraint_sop_branch(pm, nw, id, f_idx, t_idx, branch["g_fr"], branch["g_to"], branch["b_fr"], branch["b_to"], branch["br_r"], branch["br_x"], branch["f_connections"], branch["t_connections"], branch["c_rating_a"], branch["m_legs"])
end


function constraint_sop_branch(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, branch_id, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, g_fr, g_to, b_fr, b_to, r, x, f_connections, t_connections, c_rating::Vector{<:Real}, m_legs; report::Bool=true)
    # pg1_sum = JuMP.@NLexpression(pm.model, sum(pgi for pgi in _PMD.var(pm, 0, :pg, 1)))
    # pg2_sum = JuMP.@NLexpression(pm.model, sum(pgi for pgi in _PMD.var(pm, 0, :pg, 2)))
    # Rin_pu = 0.015 / (230.94^2 / 1000)
    # Vdc_pu = 600 / 1000
    # JuMP.@NLconstraint(pm.model, pg1_sum + pg2_sum == 0 )#+ (Rin_pu / Vdc_pu^2)*(pg1_sum + pg2_sum)^2 == 0)

    f_bus = f_idx[2]
    t_bus = t_idx[2]

    vr_fr = [_PMD.var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr = [_PMD.var(pm, nw, :vi, f_bus)[c] for c in f_connections]
    vr_to = [_PMD.var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to = [_PMD.var(pm, nw, :vi, t_bus)[c] for c in t_connections]

    cr_fr = _PMD.var(pm, nw, :cr, f_idx)
    ci_fr = _PMD.var(pm, nw, :ci, f_idx)
    cr_to = _PMD.var(pm, nw, :cr, t_idx)
    ci_to = _PMD.var(pm, nw, :ci, t_idx)

    csr_fr = _PMD.var(pm, nw, :csr, f_idx[1])
    csi_fr = _PMD.var(pm, nw, :csi, f_idx[1])
    csr_to = _PMD.var(pm, nw, :csr, t_idx[1])
    csi_to = _PMD.var(pm, nw, :csi, t_idx[1])

    # ### constraint_mc_current_from
    # JuMP.@constraint(pm.model, cr_fr .== csr_fr + g_fr*vr_fr - b_fr*vi_fr)
    # JuMP.@constraint(pm.model, ci_fr .== csi_fr + g_fr*vi_fr + b_fr*vr_fr)
    # _PMD.var(pm, nw, :cr_bus)[f_idx] = cr_bus_fr = _PMD._merge_bus_flows(pm, cr_fr, f_connections)
    # _PMD.var(pm, nw, :ci_bus)[f_idx] = ci_bus_fr = _PMD._merge_bus_flows(pm, ci_fr, f_connections)
    # if report
    #     _PMD.sol(pm, nw, :branch, f_idx[1])[:cr_fr] = csr_fr + g_fr*vr_fr - b_fr*vi_fr
    #     _PMD.sol(pm, nw, :branch, f_idx[1])[:ci_fr] = csi_fr + g_fr*vi_fr + b_fr*vr_fr
    #     _PMD.sol(pm, nw, :branch, f_idx[1])[:pf] =  cr_fr.*vr_fr .+ ci_fr.*vi_fr
    #     _PMD.sol(pm, nw, :branch, f_idx[1])[:qf] = -cr_fr.*vi_fr .+ ci_fr.*vr_fr
    # end

    # ### constraint_mc_current_to
    # JuMP.@constraint(pm.model, cr_to .== csr_to .+ g_to * vr_to .- b_to * vi_to)
    # JuMP.@constraint(pm.model, ci_to .== csi_to .+ g_to * vi_to .+ b_to * vr_to)
    # _PMD.var(pm, nw, :cr_bus)[t_idx] = cr_bus_to = _PMD._merge_bus_flows(pm, cr_to, t_connections)
    # _PMD.var(pm, nw, :ci_bus)[t_idx] = ci_bus_to = _PMD._merge_bus_flows(pm, ci_to, t_connections)
    # if report
    #     _PMD.sol(pm, nw, :branch, f_idx[1])[:cr_to] = csr_to .+ g_to * vr_to .- b_to * vi_to
    #     _PMD.sol(pm, nw, :branch, f_idx[1])[:ci_to] = csi_to .+ g_to * vi_to .+ b_to * vr_to
    #     _PMD.sol(pm, nw, :branch, t_idx[1])[:pt] =  cr_to.*vr_to .+ ci_to.*vi_to
    #     _PMD.sol(pm, nw, :branch, t_idx[1])[:qt] = -cr_to.*vi_to .+ ci_to.*vr_to
    # end

    
    # ### constraint_mc_bus_voltage_drop
    # # JuMP.@constraint(pm.model, vr_to .== vr_fr .- r*csr_fr .+ x*csi_fr)
    # # JuMP.@constraint(pm.model, vi_to .== vi_fr .- r*csi_fr .- x*csr_fr)
    # JuMP.@constraint(pm.model, vr_to .== vr_fr - r*csr_fr + x*csi_fr)
    # JuMP.@constraint(pm.model, vi_to .== vi_fr - r*csi_fr - x*csr_fr)


    
    ### constraint_mc_branch_current_limit
    bg = _PMD.var(pm, nw, :bg, branch_id)
    alpha_g = 1/m_legs * ones(m_legs)
    @assert size(bg,1) == length(f_connections) + length(t_connections)
    @assert size(bg,1) <= size(bg,2)

    cr_max = _PMD.var(pm, nw, :cr_max, branch_id)
    @assert size(cr_max,1) == length(f_connections) + length(t_connections)

    JuMP.@constraint(pm.model, [k in 1:size(bg,2)], sum(bg[:,k]) == 1)

    # c_rating = JuMP.@expression(pm.model, sum(c_rating) .* Array(bg) * alpha_g)
    # c_rating_fr = c_rating[1:length(f_connections)]
    # c_rating_to = c_rating[length(f_connections)+1:length(c_rating)]
    # JuMP.@constraint(pm.model, [c in 1:length(f_connections)], cr_fr[c]^2+ci_fr[c]^2 <= c_rating_fr[c]^2)
    # JuMP.@constraint(pm.model, [c in 1:length(t_connections)], cr_to[c]^2+ci_to[c]^2 <= c_rating_to[c]^2)
    
    mult = 1E10
    JuMP.@constraint(pm.model, cr_max .== sum(c_rating) * Array(bg) * alpha_g)
    JuMP.@NLconstraint(pm.model, [c in 1:length(f_connections)], mult * (cr_fr[c]^2+ci_fr[c]^2) <= mult * cr_max[c]^2)
    JuMP.@NLconstraint(pm.model, [c in 1:length(t_connections)], mult * (cr_to[c]^2+ci_to[c]^2) <= mult * cr_max[c+length(f_connections)]^2)

    ### constraint sum(csr_fr) + sum(csr_to) = 0,   sum(csi_fr) + sum(csi_to) = 0
    JuMP.@constraint(pm.model, sum(csr_fr) + sum(csr_to) == 0)
    JuMP.@constraint(pm.model, sum(csi_fr) + sum(csi_to) == 0)

    # ### constraint_mc_thermal_limit
    # if haskey(branch, "rate_a") && any(branch["rate_a"] .< Inf)
    #     ### constraint_mc_thermal_limit_from
    #     pf_idx = JuMP.@expression(model,  vr_fr .* cr_fr .+ vi_fr .* ci_fr)
    #     qf_idx = JuMP.@expression(model, -vr_fr .* ci_fr .+ vi_fr .* cr_fr)
    #     JuMP.@constraint(model, pf_idx.^2 .+ qf_idx.^2 .<= branch["rate_a"].^2)

    #     ### constraint_mc_thermal_limit_to
    #     pt_idx = JuMP.@expression(model,  vr_to .* cr_to .+ vi_to .* ci_to)
    #     qt_idx = JuMP.@expression(model, -vr_to .* ci_to .+ vi_to .* cr_to)
    #     JuMP.@constraint(model, pt_idx.^2 .+ qt_idx.^2 .<= branch["rate_a"].^2)
    # end

end



function objective_mc_min_IUF(pm::_PMD.AbstractUnbalancedPowerModel)
    alpha = exp(im*2/3*pi)
    T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
    Tre = real.(T)
    Tim = imag.(T)

    ref = pm.ref[:it][:pmd][:nw][0]   # TODO change 0 to nw, make this ref dependent
    _, _, arc, branch = get_ref_bus_branch(ref)
    # arc = (3, 1, 2)
    # branch = 3
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



function objective_mc_min_losses_active(pm::_PMD.AbstractUnbalancedPowerModel)
    # vr_fr = [_PMD.var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    # vi_fr = [_PMD.var(pm, nw, :vi, f_bus)[c] for c in f_connections]
    # vr_to = [_PMD.var(pm, nw, :vr, f_bus)[c] for c in t_connections]
    # vi_to = [_PMD.var(pm, nw, :vi, f_bus)[c] for c in t_connections]

    # cr_fr = _PMD.var(pm, nw, :cr, f_idx)
    # ci_fr = _PMD.var(pm, nw, :ci, f_idx)
    # cr_to = _PMD.var(pm, nw, :cr, f_idx)
    # ci_to = _PMD.var(pm, nw, :ci, f_idx)
    
    # pf_idx = JuMP.@expression(model,  vr_fr .* cr_fr .+ vi_fr .* ci_fr)
    # pt_idx = JuMP.@expression(model,  vr_to .* cr_to .+ vi_to .* ci_to)
    
    obj = sum(
            sum( [_PMD.var(pm, n, :vr, branch["f_bus"])[c] for c in branch["f_connections"]] .* _PMD.var(pm, n, :cr, (i, branch["f_bus"], branch["t_bus"]))
               .+[_PMD.var(pm, n, :vi, branch["f_bus"])[c] for c in branch["f_connections"]] .* _PMD.var(pm, n, :ci, (i, branch["f_bus"], branch["t_bus"]))
               .+[_PMD.var(pm, n, :vr, branch["t_bus"])[c] for c in branch["t_connections"]] .* _PMD.var(pm, n, :cr, (i, branch["t_bus"], branch["f_bus"]))
               .+[_PMD.var(pm, n, :vi, branch["t_bus"])[c] for c in branch["t_connections"]] .* _PMD.var(pm, n, :ci, (i, branch["t_bus"], branch["f_bus"]))
            for (i,branch) in nw_ref[:branch])
        for (n, nw_ref) in _PMD.nws(pm))
    # obj = sum(
    #         sum( _PMD.var(pm, n, :pf, i) + _PMD.var(pm, n, :pt, i) for (i,branch) in nw_ref[:branch])
    #     for (n, nw_ref) in _PMD.nws(pm))

    return JuMP.@objective(pm.model, Min, obj)
end


"""
	function build_mc_opf(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for OPF in current-voltage variable space with explicit neutrals
"""
function build_mc_opf_sop(pm::_PMD.AbstractExplicitNeutralIVRModel)
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
        # variable_multiplexing_inverter(pm)
        variable_multiplexing_sop_inverter(pm)
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

    # inverter_branches = [branch["index"] for (i, branch) in pm.data["branch"] if occursin("inverter_branch", branch["name"])]
    # gen_ids, gen_bus_ids, gen_branch_ids = get_pv_bus_branch(pm.ref[:it][:pmd][:nw][0])  # maybe output branch_ids and arcs seperately?
    # branch_gens = Dict(branch_id[1] => gen_ids[i] for (i, branch_id) in enumerate(gen_branch_ids))

    sop_branches = [parse(Int, i) for (i, branch) in pm.data["branch"] if occursin("SOP_branch", branch["name"])]

    for i in _PMD.ids(pm, :branch)
        _PMD.constraint_mc_current_from(pm, i)
        _PMD.constraint_mc_current_to(pm, i)
        _PMD.constraint_mc_bus_voltage_drop(pm, i)

        if i ∈ sop_branches && pm.setting["multileg"] && pm.setting["multiplexing"]
            constraint_sop_branch(pm, i)
        else
        # if (i ∈ inverter_branches) && pm.setting["multileg"] && pm.setting["multiplexing"]
        #     constraint_mc_branch_current_limit(pm, i, branch_gens[i])
        # else
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
    _PMD.objective_mc_min_fuel_cost(pm)
    # objective_mc_min_IUF(pm)
    # objective_mc_min_losses_active(pm)
end