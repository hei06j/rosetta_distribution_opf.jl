function solve_opf_acp(data, optimizer; verbose=true)
    time_data_start = time()

    ref = _IM.build_ref(data, _PMD.ref_add_core!, _PMD._pmd_global_keys, _PMD.pmd_it_name)[:it][:pmd][:nw][0]
    data_load_time = time() - time_data_start
    time_model_start = time()

    model = JuMP.Model(optimizer)
    JuMP.set_optimizer_attribute(model, "print_level", 0)
    n_ph = 3
        # JuMP.@variable(model, ref[:bus][i]["vmin"][ph] <= vm[ph in 1:n_ph, i in keys(ref[:bus])] <= ref[:bus][i]["vmax"][ph], start=1.0)
    # JuMP.@variable(model, -pi/2 <= va[ph in 1:n_ph, i in keys(ref[:bus])] <= pi/2 , start=0.0)

    JuMP.@variable(model, vm[ph in 1:n_ph, i in keys(ref[:bus])], start=1.0)  
    JuMP.@variable(model, va[ph in 1:n_ph, i in keys(ref[:bus])]) 

    va_start = collect(0:-1:-2)*2/3*pi
    for i in 1:3
        for j in keys(ref[:bus])
            JuMP.set_start_value(va[i,j], va_start[i])
        end
    end

    JuMP.@variable(model, ref[:gen][i]["pmin"][ph] <= pg[ph in 1:n_ph, i in keys(ref[:gen])] <= ref[:gen][i]["pmax"][ph])
    JuMP.@variable(model, ref[:gen][i]["qmin"][ph] <= qg[ph in 1:n_ph, i in keys(ref[:gen])] <= ref[:gen][i]["qmax"][ph])

    JuMP.@variable(model, -ref[:branch][l]["rate_a"][ph] <= p[ph in 1:n_ph, (l,i,j) in ref[:arcs_branch]] <= ref[:branch][l]["rate_a"][ph])
    JuMP.@variable(model, -ref[:branch][l]["rate_a"][ph] <= q[ph in 1:n_ph, (l,i,j) in ref[:arcs_branch]] <= ref[:branch][l]["rate_a"][ph])

    # JuMP.@objective(model, Min, sum(gen["cost"][1]*sum(pg[:,i].^2) + gen["cost"][2]*sum(pg[:,i]) for (i,gen) in ref[:gen]))
    # JuMP.@objective(model, Min, sum(gen["cost"][1]*sum(pg[:,i]) for (i,gen) in ref[:gen]))
    JuMP.@objective(model, Min, sum(gen["cost"][1]*sum(pg[:,i]) + gen["cost"][2] for (i,gen) in ref[:gen]))

    for (i,bus) in ref[:ref_buses]
        JuMP.@constraint(model, vm[:,i] .== bus["vm"])
        JuMP.@constraint(model, va[:,i] .== bus["va"])
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        JuMP.@constraint(model,
            sum(p[:,a] for a in ref[:bus_arcs_branch][i]) .==
            sum([pg[:,g] for g in ref[:bus_gens][i]], init=[0;0;0]) .-
            sum([load["pd"] for load in bus_loads], init=[0;0;0]) .-
            sum([shunt["gs"] for shunt in bus_shunts], init=[0;0;0]).*(vm[:,i].^2)
        )

        JuMP.@constraint(model,
            sum(q[:,a] for a in ref[:bus_arcs_branch][i]) .==
            sum([qg[:,g] for g in ref[:bus_gens][i]], init=[0;0;0]) .-
            sum([load["qd"] for load in bus_loads], init=[0;0;0]) .+
            sum([shunt["bs"] for shunt in bus_shunts], init=[0;0;0]).*(vm[:,i].^2)
        )
        if bus["bus_type"] != 3
            JuMP.@constraint(model, bus["vmin"] .<= vm[:, i]  .<= bus["vmax"])
        end
    end

    # Branch power flow physics and limit constraints
    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[:,f_idx]
        q_fr = q[:,f_idx]
        p_to = p[:,t_idx]
        q_to = q[:,t_idx]

        vm_fr = vm[:,branch["f_bus"]]
        vm_to = vm[:,branch["t_bus"]]
        va_fr = va[:,branch["f_bus"]]
        va_to = va[:,branch["t_bus"]]

        G, B = _PMD.calc_branch_y(branch)
        G_fr = branch["g_fr"]
        B_fr = branch["b_fr"]
        G_to = branch["g_to"]
        B_to = branch["b_to"]

        for pp in 1:3      #p
            JuMP.@NLconstraint(model, p_fr[pp] == 
                    (G[pp,pp]+G_fr[pp,pp])*vm_fr[pp]^2
                    +sum( 
                         (G[pp,qq]+G_fr[pp,qq])*vm_fr[pp]*vm_fr[qq]*cos(va_fr[pp]-va_fr[qq])
                        +(B[pp,qq]+B_fr[pp,qq])*vm_fr[pp]*vm_fr[qq]*sin(va_fr[pp]-va_fr[qq])
                        for qq in 1:3 if pp != qq)
                    +sum( 
                        -G[pp,qq]*vm_fr[pp]*vm_to[qq]*cos(va_fr[pp]-va_to[qq])
                        -B[pp,qq]*vm_fr[pp]*vm_to[qq]*sin(va_fr[pp]-va_to[qq])
                        for qq in 1:3)
                    )
            
        
            JuMP.@NLconstraint(model, q_fr[pp] == 
                    -(B[pp,pp]+B_fr[pp,pp])*vm_fr[pp]^2
                    -sum( 
                         (B[pp,qq]+B_fr[pp,qq])*vm_fr[pp]*vm_fr[qq]*cos(va_fr[pp]-va_fr[qq])
                        -(G[pp,qq]+G_fr[pp,qq])*vm_fr[pp]*vm_fr[qq]*sin(va_fr[pp]-va_fr[qq])
                        for qq in 1:3 if pp != qq)
                    -sum(
                        - B[pp,qq]*vm_fr[pp]*vm_to[qq]*cos(va_fr[pp]-va_to[qq])
                        + G[pp,qq]*vm_fr[pp]*vm_to[qq]*sin(va_fr[pp]-va_to[qq])
                        for qq in 1:3)
                    )
            JuMP.@NLconstraint(model, p_to[pp] == 
                    (G[pp,pp]+G_to[pp,pp])*vm_to[pp]^2
                    +sum( 
                         (G[pp,qq]+G_to[pp,qq])*vm_to[pp]*vm_to[qq]*cos(va_to[pp]-va_to[qq])
                        +(B[pp,qq]+B_to[pp,qq])*vm_to[pp]*vm_to[qq]*sin(va_to[pp]-va_to[qq])
                        for qq in 1:3 if pp != qq)
                    +sum( 
                        -G[pp,qq]*vm_to[pp]*vm_fr[qq]*cos(va_to[pp]-va_fr[qq])
                        -B[pp,qq]*vm_to[pp]*vm_fr[qq]*sin(va_to[pp]-va_fr[qq])
                        for qq in 1:3)
                    )
            
        
            JuMP.@NLconstraint(model, q_to[pp] == 
                    -(B[pp,pp]+B_to[pp,pp])*vm_to[pp]^2
                    -sum( 
                         (B[pp,qq]+B_to[pp,qq])*vm_to[pp]*vm_to[qq]*cos(va_to[pp]-va_to[qq])
                        -(G[pp,qq]+G_to[pp,qq])*vm_to[pp]*vm_to[qq]*sin(va_to[pp]-va_to[qq])
                        for qq in 1:3 if pp != qq)
                    -sum(
                        - B[pp,qq]*vm_to[pp]*vm_fr[qq]*cos(va_to[pp]-va_fr[qq])
                        + G[pp,qq]*vm_to[pp]*vm_fr[qq]*sin(va_to[pp]-va_fr[qq])
                        for qq in 1:3)
                    )
        end
     
        # Voltage angle difference limit
        JuMP.@constraint(model, branch["angmin"] .<= va_fr .- va_to .<= branch["angmax"])

        if haskey(branch, "rate_a") && any(branch["rate_a"] .< Inf)
        # Apparent power limit, from side and to side
            JuMP.@constraint(model, p_fr.^2 .+ q_fr.^2 .<= branch["rate_a"].^2)
            JuMP.@constraint(model, p_to.^2 .+ q_to.^2 .<= branch["rate_a"].^2)
        end
    end

    model_variables = JuMP.num_variables(model)

    # for consistency with other solvers, skip the variable bounds in the constraint count
    non_nl_constraints = sum(JuMP.num_constraints(model, ft, st) for (ft, st) in JuMP.list_of_constraint_types(model) if ft != JuMP.VariableRef)
    model_constraints = JuMP.num_nonlinear_constraints(model) + non_nl_constraints

    model_build_time = time() - time_model_start


    time_solve_start = time()

    JuMP.optimize!(model)
    cost = JuMP.objective_value(model)
    feasible = (JuMP.termination_status(model) == JuMP.LOCALLY_SOLVED)

    solve_time = time() - time_solve_start
    total_time = time() - time_data_start

    nlp_block = JuMP.MOI.get(model, JuMP.MOI.NLPBlock())
    total_callback_time =
        nlp_block.evaluator.eval_objective_timer +
        nlp_block.evaluator.eval_objective_gradient_timer +
        nlp_block.evaluator.eval_constraint_timer +
        nlp_block.evaluator.eval_constraint_jacobian_timer +
        nlp_block.evaluator.eval_hessian_lagrangian_timer

    if verbose
        println("")
        println("\033[1mSummary\033[0m")
        println("   case........: $(data["name"])")
        println("   variables...: $(model_variables)")
        println("   constraints.: $(model_constraints)")
        println("   feasible....: $(feasible)")
        println("   cost........: $(round(Int, cost))")
        println("   total time..: $(total_time)")
        println("     data time.: $(data_load_time)")
        println("     build time: $(model_build_time)")
        println("     solve time: $(solve_time)")
        println("      callbacks: $(total_callback_time)")
        println("")
        println("   callbacks time:")
        println("   * obj.....: $(nlp_block.evaluator.eval_objective_timer)")
        println("   * grad....: $(nlp_block.evaluator.eval_objective_gradient_timer)")
        println("   * cons....: $(nlp_block.evaluator.eval_constraint_timer)")
        println("   * jac.....: $(nlp_block.evaluator.eval_constraint_jacobian_timer)")
        println("   * hesslag.: $(nlp_block.evaluator.eval_hessian_lagrangian_timer)")
        println("")
    end

    return Dict(
        "case" => data["name"],
        "variables" => model_variables,
        "constraints" => model_constraints,
        "feasible" => feasible,
        "cost" => cost,
        "time_total" => total_time,
        "time_data" => data_load_time,
        "time_build" => model_build_time,
        "time_solve" => solve_time,
        "time_callbacks" => total_callback_time,
        "solution" => Dict(
            "bus" => Dict(
                "$i" => Dict( "vm" => [JuMP.value.(vm[j,i]) for j in 1:3],
                              "va" => [JuMP.value.(va[j,i]) for j in 1:3])
                for (i, bus) in ref[:bus]
            )
        )
    )
end