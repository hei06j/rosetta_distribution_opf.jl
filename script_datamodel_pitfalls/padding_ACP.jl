#!/usr/bin/env julia
###### AC-OPF using JuMP ######
#
# implementation reference: https://github.com/lanl-ansi/PowerModelsAnnex.jl/blob/master/src/model/ac-opf.jl
# only the built-in AD library is supported

import PowerModelsDistribution
import Ipopt
import JuMP
import InfrastructureModels
using rosetta_distribution_opf
import LinearAlgebra: diag
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

# file_name = "./data/case0_load_2ph_wye_cp.dss"
# data = PMD.parse_file(file_name)
# n_ph = 2

# file_name = "./data/case1_load_2ph_wye_cp.dss"
# data = PMD.parse_file(file_name)
# n_ph = 3


file_name = "./data/Network1_Feeder1_2w/Master.dss"
data = PMD.parse_file(file_name)
n_ph = 2


# file_name = "./data/Network1_Feeder1_2w_padded/Master.dss"
# data = PMD.parse_file(file_name)
# n_ph = 3


##
data["settings"]["sbase_default"] = 1.0
data["voltage_source"]["source"]["rs"] *= 0
data["voltage_source"]["source"]["xs"] *= 0

if PMD.iseng(data)
    data = PMD.transform_data_model(data)
end
data["gen"]["2"] = deepcopy( data["gen"]["1"])
data["gen"]["2"]["gen_bus"] = data["load"]["1"]["load_bus"]
data["gen"]["2"]["cost"]*=2.0
for (i, bus) in data["bus"]
    # bus["vmin"] = 0.90*ones(length(bus["vmin"]))
    # bus["vmax"] = 1.10*ones(length(bus["vmax"]))
    bus["vmin"] = 0.90*ones(n_ph)
    bus["vmax"] = 1.10*ones(n_ph)
end
for (g,gen) in data["gen"]
    gen["pmin"] =   0*ones(n_ph);
    gen["pmax"] =  20*ones(n_ph);
    gen["qmin"] = -20*ones(n_ph);
    gen["qmax"] =  20*ones(n_ph);
    gen["cost"] *= 1000
end
for (b,branch) in data["branch"]
    branch["rate_a"] = 12*ones(n_ph)
end
for (l,load) in data["load"]
    pd = zeros(n_ph); pd[load["connections"]] .= load["pd"]
    load["pd"] = pd

    qd = zeros(n_ph); qd[load["connections"]] .= load["qd"]
    load["qd"] = qd
end


time_data_start = time()

ref = IM.build_ref(data, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]
data_load_time = time() - time_data_start
time_model_start = time()

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "sb"=>"yes","warm_start_init_point"=>"yes", "max_iter"=>1000)
model = JuMP.Model(ipopt_solver)

JuMP.@variable(model, vm[ph in 1:n_ph, i in keys(ref[:bus])], start=1.0)  
JuMP.@variable(model, va[ph in 1:n_ph, i in keys(ref[:bus])]) 

JuMP.@variable(model, ref[:gen][i]["pmin"][ph] <= pg[ph in 1:n_ph, i in keys(ref[:gen])] <= ref[:gen][i]["pmax"][ph])
JuMP.@variable(model, ref[:gen][i]["qmin"][ph] <= qg[ph in 1:n_ph, i in keys(ref[:gen])] <= ref[:gen][i]["qmax"][ph])

JuMP.@variable(model, -ref[:branch][l]["rate_a"][ph] <= p[ph in 1:n_ph, (l,i,j) in ref[:arcs_branch]] <= ref[:branch][l]["rate_a"][ph])
JuMP.@variable(model, -ref[:branch][l]["rate_a"][ph] <= q[ph in 1:n_ph, (l,i,j) in ref[:arcs_branch]] <= ref[:branch][l]["rate_a"][ph])

JuMP.@objective(model, Min, sum(gen["cost"][1]*sum(pg[:,i]) + gen["cost"][2] for (i,gen) in ref[:gen]))


for (i,bus) in ref[:ref_buses]
    JuMP.@constraint(model, vm[:,i] .== bus["vm"][1:n_ph])
    JuMP.@constraint(model, va[:,i] .== bus["va"][1:n_ph])
end

for (i,bus) in ref[:bus]
    bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
    bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

    JuMP.@constraint(model,
        Array(sum(p[:,a] for a in ref[:bus_arcs_branch][i])) .==
        Array(sum([pg[:,g] for g in ref[:bus_gens][i]], init=zeros(n_ph))) .-
        Array(sum([load["pd"] for load in bus_loads], init=zeros(n_ph))) .-
        Array(sum([shunt["gs"] for shunt in bus_shunts], init=zeros(n_ph)).*(vm[1:n_ph,i].^2))
    )

    JuMP.@constraint(model,
        Array(sum(q[:,a] for a in ref[:bus_arcs_branch][i])) .==
        Array(sum([qg[:,g] for g in ref[:bus_gens][i]], init=zeros(n_ph))) .-
        Array(sum([load["qd"] for load in bus_loads], init=zeros(n_ph))) .+
        Array(sum([shunt["bs"] for shunt in bus_shunts], init=zeros(n_ph)).*(vm[1:n_ph,i].^2))
    )
    if bus["bus_type"] != 3
        JuMP.@constraint(model, bus["vmin"] .<= vm[1:n_ph, i]  .<= bus["vmax"])
    end
end

for (i,branch) in ref[:branch]
    f_idx = (i, branch["f_bus"], branch["t_bus"])
    t_idx = (i, branch["t_bus"], branch["f_bus"])

    p_fr = p[:,f_idx]
    q_fr = q[:,f_idx]
    p_to = p[:,t_idx]
    q_to = q[:,t_idx]

    vm_fr = vm[:,branch["f_bus"]][branch["f_connections"]]
    vm_to = vm[:,branch["t_bus"]][branch["t_connections"]]
    va_fr = va[:,branch["f_bus"]][branch["f_connections"]]
    va_to = va[:,branch["t_bus"]][branch["t_connections"]]

    G, B = PMD.calc_branch_y(branch)
    G_fr = branch["g_fr"]
    B_fr = branch["b_fr"]
    G_to = branch["g_to"]
    B_to = branch["b_to"]

    for pp in 1:n_ph      #p
        JuMP.@NLconstraint(model, p_fr[pp] == 
                (G[pp,pp]+G_fr[pp,pp])*vm_fr[pp]^2
                +sum( 
                        (G[pp,qq]+G_fr[pp,qq])*vm_fr[pp]*vm_fr[qq]*cos(va_fr[pp]-va_fr[qq])
                    +(B[pp,qq]+B_fr[pp,qq])*vm_fr[pp]*vm_fr[qq]*sin(va_fr[pp]-va_fr[qq])
                    for qq in 1:n_ph if pp != qq)
                +sum( 
                    -G[pp,qq]*vm_fr[pp]*vm_to[qq]*cos(va_fr[pp]-va_to[qq])
                    -B[pp,qq]*vm_fr[pp]*vm_to[qq]*sin(va_fr[pp]-va_to[qq])
                    for qq in 1:n_ph)
                )
        
    
        JuMP.@NLconstraint(model, q_fr[pp] == 
                -(B[pp,pp]+B_fr[pp,pp])*vm_fr[pp]^2
                -sum( 
                        (B[pp,qq]+B_fr[pp,qq])*vm_fr[pp]*vm_fr[qq]*cos(va_fr[pp]-va_fr[qq])
                    -(G[pp,qq]+G_fr[pp,qq])*vm_fr[pp]*vm_fr[qq]*sin(va_fr[pp]-va_fr[qq])
                    for qq in 1:n_ph if pp != qq)
                -sum(
                    - B[pp,qq]*vm_fr[pp]*vm_to[qq]*cos(va_fr[pp]-va_to[qq])
                    + G[pp,qq]*vm_fr[pp]*vm_to[qq]*sin(va_fr[pp]-va_to[qq])
                    for qq in 1:n_ph)
                )
        JuMP.@NLconstraint(model, p_to[pp] == 
                (G[pp,pp]+G_to[pp,pp])*vm_to[pp]^2
                +sum( 
                        (G[pp,qq]+G_to[pp,qq])*vm_to[pp]*vm_to[qq]*cos(va_to[pp]-va_to[qq])
                    +(B[pp,qq]+B_to[pp,qq])*vm_to[pp]*vm_to[qq]*sin(va_to[pp]-va_to[qq])
                    for qq in 1:n_ph if pp != qq)
                +sum( 
                    -G[pp,qq]*vm_to[pp]*vm_fr[qq]*cos(va_to[pp]-va_fr[qq])
                    -B[pp,qq]*vm_to[pp]*vm_fr[qq]*sin(va_to[pp]-va_fr[qq])
                    for qq in 1:n_ph)
                )
        
    
        JuMP.@NLconstraint(model, q_to[pp] == 
                -(B[pp,pp]+B_to[pp,pp])*vm_to[pp]^2
                -sum( 
                        (B[pp,qq]+B_to[pp,qq])*vm_to[pp]*vm_to[qq]*cos(va_to[pp]-va_to[qq])
                    -(G[pp,qq]+G_to[pp,qq])*vm_to[pp]*vm_to[qq]*sin(va_to[pp]-va_to[qq])
                    for qq in 1:n_ph if pp != qq)
                -sum(
                    - B[pp,qq]*vm_to[pp]*vm_fr[qq]*cos(va_to[pp]-va_fr[qq])
                    + G[pp,qq]*vm_to[pp]*vm_fr[qq]*sin(va_to[pp]-va_fr[qq])
                    for qq in 1:n_ph)
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

res = JuMP.optimize!(model)
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

    
verbose = true
if verbose
    println("")
    println("\033[1mSummary\033[0m")
    println("   case........: $(data["name"])")
    println("   variables...: $(model_variables)")
    println("   constraints.: $(model_constraints)")
    println("   feasible....: $(feasible)")
    # println("   cost........: %$(round(Int, cost))")
    println("   cost........: $cost")
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
            "$i" => Dict( "vm" => [JuMP.value.(vm[j,i]) for j in 1:n_ph],
                            "va" => [JuMP.value.(va[j,i]) for j in 1:n_ph])
            for (i, bus) in ref[:bus]
        )
    )
)

