"""
objective_mc_min_losses_active(pm::PMD.AbstractUnbalancedPowerModel)
    Minimising negative current sequence on the main feeder branch
"""
function objective_mc_min_IUF(pm::PMD.AbstractUnbalancedPowerModel)
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
    
    cr_bus = PMD.var(pm, 0, :cr_bus, arc)  # TODO change 0 to nw
    ci_bus = PMD.var(pm, 0, :ci_bus, arc)  # TODO change 0 to nw

    phases = 1:3
    JuMP.@constraint(pm.model, cr_012[phases,:] .== Tre * Array(cr_bus[phases]) .- Tim * Array(ci_bus[phases]))
    JuMP.@constraint(pm.model, ci_012[phases,:] .== Tre * Array(ci_bus[phases]) .+ Tim * Array(cr_bus[phases]))
    JuMP.@constraint(pm.model, cm_012[phases,:].^2 .== cr_012[phases,:].^2 .+ ci_012[phases,:].^2)

    # JuMP.@objective(model, Min, cm_012[3,arc] / cm_012[2,arc])  # minimising ratio of negative to positive sequence
    JuMP.@objective(pm.model, Min, cm_012[3,arc])  # minimising negative sequence
end


"""
objective_mc_min_losses_active(pm::PMD.AbstractUnbalancedPowerModel)
    Minimising totale active losses in the network
"""
function objective_mc_min_losses_active(pm::PMD.AbstractUnbalancedPowerModel)
    obj = sum(
            sum( [PMD.var(pm, n, :vr, branch["f_bus"])[c] for c in branch["f_connections"]] .* PMD.var(pm, n, :cr, (i, branch["f_bus"], branch["t_bus"]))
               .+[PMD.var(pm, n, :vi, branch["f_bus"])[c] for c in branch["f_connections"]] .* PMD.var(pm, n, :ci, (i, branch["f_bus"], branch["t_bus"]))
               .+[PMD.var(pm, n, :vr, branch["t_bus"])[c] for c in branch["t_connections"]] .* PMD.var(pm, n, :cr, (i, branch["t_bus"], branch["f_bus"]))
               .+[PMD.var(pm, n, :vi, branch["t_bus"])[c] for c in branch["t_connections"]] .* PMD.var(pm, n, :ci, (i, branch["t_bus"], branch["f_bus"]))
            for (i,branch) in nw_ref[:branch])
        for (n, nw_ref) in PMD.nws(pm))

    return JuMP.@objective(pm.model, Min, obj)
end



"""
objective_mc_min_max_phase_current(pm::PMD.AbstractUnbalancedPowerModel)
    Minimising maximum current magnitude on the main feeder branch
"""
function objective_mc_min_max_phase_current(pm::PMD.AbstractUnbalancedPowerModel)
    ref = pm.ref[:it][:pmd][:nw][0]       # TODO change 0 to nw, make this ref dependent
    _, _, arc, branch = get_ref_bus_branch(ref)
    nconds = Dict(l => length(ref[:branch][l]["f_connections"]) for l in [branch])
    conds = Dict(l => ref[:branch][l]["f_connections"] for l in [branch])
    
    cr_bus = PMD.var(pm, 0, :cr_bus, arc)  # TODO change 0 to nw
    ci_bus = PMD.var(pm, 0, :ci_bus, arc)  # TODO change 0 to nw
    n_ph = 4                               # TODO, make this dependent on nconds, and set to zero the ones that are not equal to nconds[i]
    conds = 1:n_ph

    ### TODO these variables and constraints should be defiened outside of here, but only added to the model whithin this objective function, 
    ### so that there are only added to the model if this objective is called
    cm_max = JuMP.@variable(pm.model, base_name="cm_max", lower_bound=0)
    
    JuMP.@constraint(pm.model, [i=1:n_ph], cm_max .>= cr_bus[i].^2 .+ ci_bus[i].^2)

    return JuMP.@objective(pm.model, Min, cm_max)
end


"""
objective_mc_min_ref_branch_loss(pm::PMD.AbstractUnbalancedPowerModel)
    Minimising active power losses on the main feeder branch
"""
function objective_mc_min_ref_branch_loss(pm::PMD.AbstractUnbalancedPowerModel)
    ref = pm.ref[:it][:pmd][:nw][0]   # TODO change 0 to nw, make this ref dependent
    _, _, arc, branch_id = get_ref_bus_branch(ref)
    branch = ref[:branch][branch_id]

    ref_branch_loss = sum(
        sum(
                [PMD.var(pm, n, :vr, branch["f_bus"])[c] for c in branch["f_connections"]] .* PMD.var(pm, n, :cr, arc)
                .+[PMD.var(pm, n, :vi, branch["f_bus"])[c] for c in branch["f_connections"]] .* PMD.var(pm, n, :ci, arc)
                .+[PMD.var(pm, n, :vr, branch["t_bus"])[c] for c in branch["t_connections"]] .* PMD.var(pm, n, :cr, arc)
                .+[PMD.var(pm, n, :vi, branch["t_bus"])[c] for c in branch["t_connections"]] .* PMD.var(pm, n, :ci, arc)
            for (n, nw_ref) in PMD.nws(pm)))

    JuMP.@objective(pm.model, Min, ref_branch_loss)
end