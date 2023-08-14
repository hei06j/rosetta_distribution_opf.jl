#!/usr/bin/env julia
###### AC-OPF using JuMP ######
#
# implementation reference: https://github.com/lanl-ansi/PowerModelsAnnex.jl/blob/master/src/model/ac-opf.jl
# only the built-in AD library is supported
#

import PowerModelsDistribution
import Ipopt
import JuMP
import InfrastructureModels
import LinearAlgebra: diag
const PMD = PowerModelsDistribution
const IM = InfrastructureModels

file_name = "./data/case3_unbalanced.dss"

# function solve_opf(file_name)
    time_data_start = time()

    data = PMD.parse_file(file_name)
    data["settings"]["sbase_default"] = 1.0

    #remove voltage source internal impedance branch
    data["voltage_source"]["source"]["rs"]*=0 
    data["voltage_source"]["source"]["xs"]*=0 

    if PMD.iseng(data)
        data = PMD.transform_data_model(data)
    end
    data["gen"]["2"] = deepcopy( data["gen"]["1"])
    data["gen"]["2"]["gen_bus"] = data["load"]["1"]["load_bus"]
    data["gen"]["2"]["cost"]*=2.0
    for (i, bus) in data["bus"]
        bus["vmin"] = 0.90*ones(3)
        bus["vmax"] = 1.10*ones(3)
    end
    for (g,gen) in data["gen"]
        gen["pmin"] =   0*ones(3);
        gen["pmax"] =  20*ones(3);
        gen["qmin"] = -20*ones(3);
        gen["qmax"] =  20*ones(3);
    end
    for (b,branch) in data["branch"]
        branch["rate_a"] = 12*ones(3)
    end
    ref = IM.build_ref(data, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]
    data_load_time = time() - time_data_start
    time_model_start = time()

    model = JuMP.Model(Ipopt.Optimizer)
    #JuMP.set_optimizer_attribute(model, "print_level", 0)
    n_ph = 3
    JuMP.@variable(model, vr[ph in 1:n_ph, i in keys(ref[:bus])])  
    JuMP.@variable(model, vi[ph in 1:n_ph, i in keys(ref[:bus])]) 

    v_start = exp.(im.*collect(0:-1:-2)*2/3*pi)
    for i in 1:3
        for j in keys(ref[:bus])
            JuMP.set_start_value(vr[i,j], real(v_start[i]))
            JuMP.set_start_value(vi[i,j], imag(v_start[i]))
        end
    end

    # JuMP.@variable(model, ref[:bus][i]["vmin"][ph] <= vm[ph in 1:n_ph, i in keys(ref[:bus])] <= ref[:bus][i]["vmax"][ph], start=1.0)
    # JuMP.@variable(model, -pi/2 <= va[ph in 1:n_ph, i in keys(ref[:bus])] <= pi/2 , start=0.0)

    JuMP.@variable(model, ref[:gen][i]["pmin"][ph] <= pg[ph in 1:n_ph, i in keys(ref[:gen])] <= ref[:gen][i]["pmax"][ph])
    JuMP.@variable(model, ref[:gen][i]["qmin"][ph] <= qg[ph in 1:n_ph, i in keys(ref[:gen])] <= ref[:gen][i]["qmax"][ph])

    JuMP.@variable(model, -ref[:branch][l]["rate_a"][ph] <= p[ph in 1:n_ph, (l,i,j) in ref[:arcs_branch]] <= ref[:branch][l]["rate_a"][ph])
    JuMP.@variable(model, -ref[:branch][l]["rate_a"][ph] <= q[ph in 1:n_ph, (l,i,j) in ref[:arcs_branch]] <= ref[:branch][l]["rate_a"][ph])

    JuMP.@objective(model, Min, sum(gen["cost"][1]*sum(pg[:,i].^2) + gen["cost"][2]*sum(pg[:,i]) for (i,gen) in ref[:gen]))




    for (i,bus) in ref[:ref_buses]
        vref = bus["vm"] .*exp.(im*bus["va"])
        vrefre = real.(vref)
        vrefim = imag.(vref)
        JuMP.@constraint(model, vr[:,i] .== vrefre)
        JuMP.@constraint(model, vr[:,i] .== vrefim)
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        JuMP.@constraint(model,
            sum(p[:,a] for a in ref[:bus_arcs_branch][i]) .==
            sum([pg[:,g] for g in ref[:bus_gens][i]], init=[0;0;0]) .-
            sum([load["pd"] for load in bus_loads], init=[0;0;0]) .-
            sum([shunt["gs"] for shunt in bus_shunts], init=[0;0;0]).*(vr[:,i].^2 + vi[:,i].^2)
        )

        JuMP.@constraint(model,
            sum(q[:,a] for a in ref[:bus_arcs_branch][i]) .==
            sum([qg[:,g] for g in ref[:bus_gens][i]], init=[0;0;0]) .-
            sum([load["qd"] for load in bus_loads], init=[0;0;0]) .+
            sum([shunt["bs"] for shunt in bus_shunts], init=[0;0;0]).*(vr[:,i].^2 + vi[:,i].^2)
        )
        if bus["bus_type"] != 3
            JuMP.@constraint(model, bus["vmin"].^2 .<= vr[:, i].^2 + vi[:, i].^2  .<= bus["vmax"].^2)
        end
    end

    # Branch power flow physics and limit constraints
    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = [p[:,f_idx]...]
        q_fr = [q[:,f_idx]...]
        p_to = [p[:,t_idx]...]
        q_to = [q[:,t_idx]...]

        vm_fr = [vm[:,branch["f_bus"]]...]
        vm_to = [vm[:,branch["t_bus"]]...]
        va_fr = [va[:,branch["f_bus"]]...]
        va_to = [va[:,branch["t_bus"]]...]

        vr_fr = [vr[:,branch["f_bus"]]...]
        vi_fr = [vi[:,branch["f_bus"]]...]
        vr_to = [vr[:,branch["t_bus"]]...]
        vi_to = [vi[:,branch["t_bus"]]...]


        g, b = PMD.calc_branch_y(branch)
        # tr, ti = PowerModels.calc_branch_t(branch)
        # ttm = tr^2 + ti^2
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]

        # From side of the branch flow
        JuMP.@constraint(model, p_fr .== diag(
            (vr_fr*vr_fr' + vi_fr*vi_fr')*g_fr'
            + 
            (vi_fr*vr_fr' - vr_fr*vi_fr')*b_fr'
            +
            (vr_fr*(vr_fr-vr_to)' + vi_fr*(vi_fr-vi_to)')*g'
            +
            (vi_fr*(vr_fr-vr_to)' - vr_fr*(vi_fr-vi_to)')*b'
            ))

        JuMP.@constraint(model, q_fr .== diag(
            -(vr_fr*vr_fr' + vi_fr*vi_fr')*b_fr'
            + 
            (vi_fr*vr_fr' - vr_fr*vi_fr')*g_fr'
            -
            (vr_fr*(vr_fr-vr_to)' + vi_fr*(vi_fr-vi_to)')*b'
            +
            (vi_fr*(vr_fr-vr_to)' - vr_fr*(vi_fr-vi_to)')*g'
            ))
        
            JuMP.@constraint(model, p_to .== diag(
            (vr_to*vr_to' + vi_to*vi_to')*g_to'
            + 
            (vi_to*vr_to' - vr_to*vi_to')*b_to'
            +
            (vr_to*(vr_to-vr_fr)' + vi_to*(vi_to-vi_fr)')*g'
            +
            (vi_to*(vr_to-vr_fr)' - vr_to*(vi_to-vi_fr)')*b'
            ))

        JuMP.@constraint(model, q_to .== diag(
            -(vr_to*vr_to' + vi_to*vi_to')*b_to'
            + 
            (vi_to*vr_to' - vr_to*vi_to')*g_to'
            -
            (vr_to*(vr_to-vr_fr)' + vi_to*(vi_to-vi_fr)')*b'
            +
            (vi_to*(vr_to-vr_fr)' - vr_to*(vi_to-vi_fr)')*g'
            ))

        # Voltage angle difference limit
        # JuMP.@constraint(model, branch["angmin"] .<= va_fr .- va_to .<= branch["angmax"])

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

    println("")
    println("\033[1mSummary\033[0m")
    println("   case........: $(file_name)")
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

    return Dict(
        "case" => file_name,
        "variables" => model_variables,
        "constraints" => model_constraints,
        "feasible" => feasible,
        "cost" => cost,
        "time_total" => total_time,
        "time_data" => data_load_time,
        "time_build" => model_build_time,
        "time_solve" => solve_time,
        "time_callbacks" => total_callback_time,
    )
# end

# # if isinteractive() == false
#     result = solve_opf(file_name)
# # end