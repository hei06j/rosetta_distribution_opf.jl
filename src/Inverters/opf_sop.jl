function variable_reconfigurable_sop_inverter(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    sop_branches = [i for (i, branch) in pm.data["branch"] if occursin("SOP_branch", branch["name"])]
    # f_connections = Dict(i => pm.data["branch"]["$i"]["f_connections"] for i in sop_branches)
    # t_connections = Dict(i => pm.data["branch"]["$i"]["f_connections"] for i in sop_branches)
    connections = Dict(i => length(pm.data["branch"]["$i"]["f_connections"]) + length(pm.data["branch"]["$i"]["t_connections"]) for i in sop_branches)
    m_legs = Dict(i => pm.data["branch"]["$i"]["m_legs"] for i in sop_branches) 
    # TODO add assert to make sure that m_legs >= connections for each sop, unless the neutral is not linked to legs and directly connected to each other
    
    bg = PMD.var(pm, nw)[:bg] = Dict(parse(Int, i) => JuMP.@variable(pm.model, [1:connections[i], 1:m_legs[i]], base_name="bg_$i", Bin) for i in sop_branches)
    for i in sop_branches
        report && sol_component_value_comp(pm, PMD.pmd_it_sym, nw, :branch, :bg, parse.(Int, i), bg[parse.(Int, i)])
    end

    PMD.var(pm, 0)[:c_rating] = Dict{Int, Any}()
end


function constraint_sop_branch_current_limit(pm::PMD.ExplicitNeutralModels, id::Int; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    branch = PMD.ref(pm, nw, :branch, id)
    f_idx = (id, branch["f_bus"], branch["t_bus"])
    t_idx = (id, branch["t_bus"], branch["f_bus"])
    constraint_sop_branch_current_limit(pm, nw, id, f_idx, t_idx, branch["f_connections"], branch["t_connections"], branch["c_rating_a"], branch["m_legs"])
end


function constraint_sop_branch_current_limit(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, branch_id, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections, t_connections, c_rating_a::Vector{<:Real}, m_legs; report::Bool=true)
    ### constraint_mc_branch_current_limit
    bg = PMD.var(pm, nw, :bg, branch_id)
    alpha_g = 1/m_legs * ones(m_legs)
    @assert size(bg,1) == length(f_connections) + length(t_connections)
    @assert size(bg,1) <= size(bg,2)

    cr_fr = PMD.var(pm, nw, :cr, f_idx)
    ci_fr = PMD.var(pm, nw, :ci, f_idx)
    cr_to = PMD.var(pm, nw, :cr, t_idx)
    ci_to = PMD.var(pm, nw, :ci, t_idx)

    # c_rating = PMD.var(pm, nw, :c_rating, branch_id)
    # @assert size(c_rating,1) == length(f_connections) + length(t_connections)

    JuMP.@constraint(pm.model, [k in 1:size(bg,2)], sum(bg[:,k]) == 1)
    # JuMP.@constraint(pm.model, [l in 1:size(bg,1)], sum(bg[l,:]) >= 1)

    c_rating = Vector{JuMP.AffExpr}([])

    for idx in collect(1:size(bg,1))
        push!(c_rating, JuMP.@expression(pm.model,  sum(c_rating_a) * sum(bg[idx,k]*alpha_g[k] for k in collect(1:length(alpha_g))) ))
    end

    JuMP.@constraint(pm.model, [c in 1:length(f_connections)], cr_fr[c]^2+ci_fr[c]^2 <= c_rating[c]^2)
    JuMP.@constraint(pm.model, [c in 1:length(t_connections)], cr_to[c]^2+ci_to[c]^2 <= c_rating[c+length(f_connections)]^2)

    PMD.var(pm, nw, :c_rating)[branch_id] = c_rating
    if report
        PMD.sol(pm, nw, :branch, branch_id)[:c_rating] = c_rating
    end

    # c_rating = JuMP.@expression(pm.model, sum(c_rating) .* Array(bg) * alpha_g)
    # c_rating_fr = c_rating[1:length(f_connections)]
    # c_rating_to = c_rating[length(f_connections)+1:length(c_rating)]
    # JuMP.@constraint(pm.model, [c in 1:length(f_connections)], cr_fr[c]^2+ci_fr[c]^2 <= c_rating_fr[c]^2)
    # JuMP.@constraint(pm.model, [c in 1:length(t_connections)], cr_to[c]^2+ci_to[c]^2 <= c_rating_to[c]^2)
end


function constraint_sop_branch_balance(pm::PMD.ExplicitNeutralModels, id::Int; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    branch = PMD.ref(pm, nw, :branch, id)
    f_idx = (id, branch["f_bus"], branch["t_bus"])
    t_idx = (id, branch["t_bus"], branch["f_bus"])
    # r = 0.015 / (230.94^2 / 1000)
    # x = 0
    # c_rating = 3 * 0.0033  # minimum(pmax1, pmax2)
    # m_legs = branch["m_legs"]
    constraint_sop_branch_balance(pm, nw, id, f_idx, t_idx, branch["f_connections"], branch["t_connections"], branch["c_rating_a"], branch["m_legs"])
end

function constraint_sop_branch_balance(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, branch_id, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections, t_connections, c_rating_a::Vector{<:Real}, m_legs; report::Bool=true)
    # pg1_sum = JuMP.@expression(pm.model, sum(pgi for pgi in PMD.var(pm, 0, :pg, 1)))
    # pg2_sum = JuMP.@expression(pm.model, sum(pgi for pgi in PMD.var(pm, 0, :pg, 2)))
    # Rin_pu = 0.015 / (230.94^2 / 1000)
    # Vdc_pu = 600 / 1000
    # JuMP.@constraint(pm.model, pg1_sum + pg2_sum == 0 )#+ (Rin_pu / Vdc_pu^2)*(pg1_sum + pg2_sum)^2 == 0)

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
    JuMP.@constraint(pm.model, sum(cr_fr) + sum(cr_to) == 0)
    JuMP.@constraint(pm.model, sum(ci_fr) + sum(ci_to) == 0)
    # JuMP.@constraint(pm.model, sum(csr_fr) + sum(csr_to) == 0)
    # JuMP.@constraint(pm.model, sum(csi_fr) + sum(csi_to) == 0)

    PMD.var(pm, 0)[:pf_idx] = Dict{Int, Any}()
    PMD.var(pm, 0)[:pt_idx] = Dict{Int, Any}()

    pf_idx = JuMP.@expression(pm.model,  vr_fr .* cr_fr .+ vi_fr .* ci_fr)
    pt_idx = JuMP.@expression(pm.model,  vr_to .* cr_to .+ vi_to .* ci_to)
    JuMP.@constraint(pm.model, sum(pf_idx) + sum(pt_idx) == 0)
    
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


function constraint_mc_sop_dc_link_ripple_power(pm::PMD.ExplicitNeutralModels, id::Int; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    branch = PMD.ref(pm, nw, :branch, id)
    f_idx = (id, branch["f_bus"], branch["t_bus"])
    t_idx = (id, branch["t_bus"], branch["f_bus"])
    # r = 0.015 / (230.94^2 / 1000)
    # x = 0
    # c_rating = 3 * 0.0033  # minimum(pmax1, pmax2)
    # m_legs = branch["m_legs"]

    pdcmin = get(branch, "pdcmin", 0)
    pdcmax = get(branch, "pdcmax", 100)

    constraint_mc_sop_dc_link_ripple_power(pm, nw, id, f_idx, t_idx, branch["f_connections"], branch["t_connections"], pdcmin, pdcmax)
end


function constraint_mc_sop_dc_link_ripple_power(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, branch_id, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections, t_connections, pdcmin::Real, pdcmax::Real; report::Bool=true)
    cr_fr = PMD.var(pm, nw, :cr, f_idx)
    ci_fr = PMD.var(pm, nw, :ci, f_idx)
    cr_to = PMD.var(pm, nw, :cr, t_idx)
    ci_to = PMD.var(pm, nw, :ci, t_idx)

    vr_fr = [PMD.var(pm, nw, :vr, f_idx[2])[t] for t in f_connections]
    vi_fr = [PMD.var(pm, nw, :vi, f_idx[2])[t] for t in f_connections]
    vr_to = [PMD.var(pm, nw, :vr, t_idx[2])[t] for t in t_connections]
    vi_to = [PMD.var(pm, nw, :vi, t_idx[2])[t] for t in t_connections]
    
    # if pdcmax > 0
    #     if pdcmax < Inf
    #         JuMP.@constraint(pm.model, pdcmax^2 >= 
    #                         sum( (vr_fr[p]*cr_fr[idx] - vi_fr[p]*ci_fr[idx])^2 + (vr_fr[p]*ci_fr[idx] + vi_fr[p]*cr_fr[idx])^2 
    #                             for (idx, p) in enumerate(f_connections))
    #                             +
    #                         sum( (vr_to[p]*cr_to[idx] - vi_to[p]*ci_to[idx])^2 + (vr_to[p]*ci_to[idx] + vi_to[p]*cr_to[idx])^2 
    #                             for (idx, p) in enumerate(t_connections))
    #                         )
    #     end
    # else
    #     JuMP.@constraint(pm.model, sum( (vr_fr[p]*cr_fr[idx] - vi_fr[p]*ci_fr[idx])^2 + (vr_fr[p]*ci_fr[idx] + vi_fr[p]*cr_fr[idx])^2 
    #         for (idx, p) in enumerate(f_connections)) == 0)
            
    #     JuMP.@constraint(pm.model, sum( (vr_to[p]*cr_to[idx] - vi_to[p]*ci_to[idx])^2 + (vr_to[p]*ci_to[idx] + vi_to[p]*cr_to[idx])^2 
    #         for (idx, p) in enumerate(t_connections)) == 0)
    # end

    pdc_link_sqr = JuMP.@expression(pm.model,  
        sum( (vr_fr[p]*cr_fr[idx] - vi_fr[p]*ci_fr[idx])^2 + (vr_fr[p]*ci_fr[idx] + vi_fr[p]*cr_fr[idx])^2 
                for (idx, p) in enumerate(f_connections))
            +
        sum( (vr_to[p]*cr_to[idx] - vi_to[p]*ci_to[idx])^2 + (vr_to[p]*ci_to[idx] + vi_to[p]*cr_to[idx])^2 
            for (idx, p) in enumerate(t_connections))
        )

    PMD.var(pm, nw, :pdc_link_sqr)[branch_id] = pdc_link_sqr
    JuMP.@constraint(pm.model, pdc_link_sqr <= pdcmax^2)
    JuMP.@constraint(pm.model, pdc_link_sqr >= 0)
    if report
        PMD.sol(pm, nw, :branch, branch_id)[:pdc_link_sqr] = pdc_link_sqr
    end
end


function induction_motor_derating(pm)
    a0 = 0.033125
    a1 = 2.75
    a2 = 56.25

    # vmnegsqr = PMD.var(pm, 0, :vmnegsqr)[bus_id]
    # vmneg    = PMD.var(pm, 0, :vmneg)[bus_id]
    Db_curve(vmneg, vmnegsqr) = vmneg <= 0.01 ? 100 : 
            ((vmneg >= 0.01 && vmneg <= 0.05) ? 100 - a2 * vmnegsqr + a1 * vmneg + a0 : 
            0)
    # JuMP.@operator(pm.model, Db, 2, Db_curve)
    JuMP.add_nonlinear_operator(pm.model, 2, Db_curve; name = :Db)

    JuMP.register(pm.model, :Db, 2, Db_curve, autodiff=true)

    # _ = JuMP.add_nonlinear_operator(pm.model, 2, Db_curve; name = :Db)

    # Db(x,y) = JuMP.NonlinearExpr(:Db, Any[x], Any[y])

end


function objective_utilize_ripple(pm)
    a0 = 0.033125 * 1
    a1 = 2.75 * 1
    a2 = 56.25 * 1
    Db_curve(vmneg, vmnegsqr) = vmneg <= 0.01 ? 1.0 : 
            ((vmneg >= 0.01 && vmneg <= 0.05) ? 1.0 - (a2 * vmnegsqr + a1 * vmneg - a0) : 
            0.0)
    JuMP.@operator(pm.model, Db, 2, Db_curve)
    # JuMP.add_nonlinear_operator(pm.model, 2, Db_curve; name = :Db)
    
    IM_load_ids = [(load["load_bus"], sum(sqrt.(load["pd"].^2 .+ load["qd"].^2))) 
        for (i, load) in PMD.ref(pm, 0, :load) if startswith(load["name"], "IM")]
    # @show IM_load_ids
    # for (bus_id, sd) in IM_load_ids
    #     @show PMD.var(pm, 0, :vmneg)[bus_id]
    #     @show PMD.var(pm, 0, :vmnegsqr)[bus_id]
    #     @show Db(PMD.var(pm, 0, :vmneg)[bus_id], PMD.var(pm, 0, :vmnegsqr)[bus_id])
    #     @show sd
    # end
    
    induction_obj = JuMP.@expression(pm.model,   
        sum(sd * (1 - Db(PMD.var(pm, 0, :vmneg)[bus_id], PMD.var(pm, 0, :vmnegsqr)[bus_id])) 
            for (bus_id, sd) in IM_load_ids)
        )
    
    ripple_obj =  JuMP.@expression(pm.model, sum(pdclinksqr for (id, pdclinksqr) in PMD.var(pm, 0)[:pdc_link_sqr]))
    # ripple_obj =  JuMP.@expression(pm.model, PMD.var(pm, 0)[:pdc_link_sqr][110])
    # ripple_obj =  JuMP.@expression(pm.model, PMD.var(pm, 0)[:pdc_link_sqr][192])
    alpha = 0.0001

    beta = 1
    
    JuMP.@objective(pm.model, Min, beta * induction_obj + alpha * ripple_obj)
    # JuMP.@objective(pm.model, Min, induction_obj)
end


"""
	function build_mc_opf(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for OPF in current-voltage variable space with explicit neutrals
"""
function build_mc_opf_mx_sop(pm::PMD.AbstractExplicitNeutralIVRModel)
    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_load_current(pm)
    PMD.variable_mc_load_power(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_generator_power(pm)
    PMD.variable_mc_transformer_current(pm)
    PMD.variable_mc_transformer_power(pm)
    PMD.variable_mc_switch_current(pm)

    if pm.setting["reconfigurable"]
        # variable_reconfigurable_inverter(pm)
        variable_reconfigurable_sop_inverter(pm)
        # PMD.var(pm, 0)[:c_rating] = Dict{Int, Any}()
    end

    if pm.setting["dc_link"]
        PMD.var(pm, 0)[:pdc_link_sqr] = Dict{Int, Any}()
    end

    if pm.setting["induction_motor"]
        PMD.var(pm, 0)[:vmneg] = Dict{Int, Any}()
        PMD.var(pm, 0)[:vmnegsqr] = Dict{Int, Any}()
        # induction_motor_derating(pm)
    end

    # Constraints
    for i in PMD.ids(pm, :bus)

        if i in PMD.ids(pm, :ref_buses)
            PMD.constraint_mc_voltage_reference(pm, i)
        end

        PMD.constraint_mc_voltage_absolute(pm, i)
        PMD.constraint_mc_voltage_pairwise(pm, i)

        if pm.setting["induction_motor"]
            constraint_mc_bus_voltage_balance(pm, i)
        end
    end
    
    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
        PMD.constraint_mc_generator_current(pm, id)
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

    # inverter_branches = [branch["index"] for (i, branch) in pm.data["branch"] if occursin("inverter_branch", branch["name"])]
    # gen_ids, gen_bus_ids, gen_branch_ids = get_pv_bus_branch(pm.ref[:it][:pmd][:nw][0])  # maybe output branch_ids and arcs seperately?
    # branch_gens = Dict(branch_id[1] => gen_ids[i] for (i, branch_id) in enumerate(gen_branch_ids))

    sop_branches = [parse(Int, i) for (i, branch) in pm.data["branch"] if occursin("SOP_branch", branch["name"])]
    
    for i in PMD.ids(pm, :branch)
        # PMD.constraint_mc_current_from(pm, i)
        # PMD.constraint_mc_current_to(pm, i)
        # PMD.constraint_mc_bus_voltage_drop(pm, i)
        
        if i ∈ sop_branches

            if pm.setting["reconfigurable"]
                constraint_mc_current_from(pm, i)
                constraint_mc_current_to(pm, i)
                # PMD.constraint_mc_bus_voltage_drop(pm, i)
                constraint_sop_branch_current_limit(pm, i)
                constraint_sop_branch_balance(pm, i)
                ## TODO why not having PMD.constraint_mc_thermal_limit_from and PMD.constraint_mc_thermal_limit_to

            elseif pm.setting["ideal"]
                constraint_mc_current_from(pm, i)
                constraint_mc_current_to(pm, i)
                # PMD.constraint_mc_bus_voltage_drop(pm, i)
                constraint_mc_branch_current_limit(pm, i)
                constraint_sop_branch_balance(pm, i)
                PMD.constraint_mc_thermal_limit_from(pm, i)  # TODO check if it works for ideal inverters
                PMD.constraint_mc_thermal_limit_to(pm, i)    # TODO check if it works for ideal inverters

            elseif pm.setting["conventional"]
                constraint_mc_current_from(pm, i)
                constraint_mc_current_to(pm, i)
                # PMD.constraint_mc_bus_voltage_drop(pm, i)
                PMD.constraint_mc_branch_current_limit(pm, i)
                constraint_sop_branch_balance(pm, i)
                # PMD.constraint_mc_thermal_limit_from(pm, i)
                # PMD.constraint_mc_thermal_limit_to(pm, i)
            end

            if pm.setting["dc_link"]
                constraint_mc_sop_dc_link_ripple_power(pm, i)
            end

        else  # normal branch
            PMD.constraint_mc_current_from(pm, i)
            PMD.constraint_mc_current_to(pm, i)
            PMD.constraint_mc_bus_voltage_drop(pm, i)
            PMD.constraint_mc_branch_current_limit(pm, i)
            PMD.constraint_mc_thermal_limit_from(pm, i)
            PMD.constraint_mc_thermal_limit_to(pm, i)
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
    if pm.setting["induction_motor"]
        objective_utilize_ripple(pm)
    else
        PMD.objective_mc_min_fuel_cost(pm)
    end
    # objective_mc_min_IUF(pm)
    # objective_mc_min_losses_active(pm)
end