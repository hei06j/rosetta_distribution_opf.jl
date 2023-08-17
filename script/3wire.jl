#!/usr/bin/env julia
###### AC-OPF using JuMP ######
#
# implementation reference: https://github.com/lanl-ansi/PowerModelsAnnex.jl/blob/master/src/model/ac-opf.jl
# only the built-in AD library is supported
#

import PowerModelsDistribution
import Ipopt
import JuMP
using jump_pmd_ivr
const PMD = PowerModelsDistribution
const RPMD = jump_pmd_ivr

file_name = "./data/case3_unbalanced.dss"

function get_data_math(file_name)
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
        gen["cost"] *= 1000
    end
    for (b,branch) in data["branch"]
        branch["rate_a"] = 12*ones(3)
    end
    return data
end

function get_vm_max_error(solution1, solution2, data)
    vm_gap_acr = [abs.(solution1["bus"][i]["vm"] - solution2["bus"][i]["vm"]) for (i, bus) in data["bus"]]
    return maximum(maximum(vm_gap_acr))
end


data = get_data_math(file_name)
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")

## run rosetta PMD ACR
results_rosetta_acr = RPMD.solve_opf_acr(data, ipopt_solver)

pm_acr = PMD.instantiate_mc_model(data, PMD.ACRUPowerModel, PMD.build_mc_opf)
for (i, bus) in PMD.ref(pm_acr, 0, :ref_buses)
    vref = bus["vm"] .* exp.(im*bus["va"])
    vrefre = real.(vref)
    vrefim = imag.(vref)
    JuMP.@constraint(pm_acr.model, PMD.var(pm_acr, 0, :vr, i) .== vrefre)
    JuMP.@constraint(pm_acr.model, PMD.var(pm_acr, 0, :vi, i) .== vrefim)
end
results_pmd_acr = PMD.optimize_model!(pm_acr, optimizer=ipopt_solver; solution_processors=[PMD.sol_data_model!])

get_vm_max_error(results_rosetta_acr["solution"], results_pmd_acr["solution"], data)

## run rosetta PMD ACP
results_rosetta_acp = RPMD.solve_opf_acp(data, ipopt_solver)

pm_acp = PMD.instantiate_mc_model(data, PMD.ACPUPowerModel, PMD.build_mc_opf)
for (i, bus) in PMD.ref(pm_acp, 0, :ref_buses)
    JuMP.@constraint(pm_acp.model, PMD.var(pm_acp, 0, :vm, i) .== bus["vm"])
    # JuMP.@constraint(pm_acp.model, PMD.var(pm_acp, 0, :va, i) .== bus["va"])
end
results_pmd_acp = PMD.optimize_model!(pm_acp, optimizer=ipopt_solver; solution_processors=[PMD.sol_data_model!])

get_vm_max_error(results_rosetta_acp["solution"], results_pmd_acp["solution"], data)
