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

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")

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

rosetta_acr = Dict()
rosetta_acp = Dict()
data_dir = "./data/ENWL/enwl_phase_to_neutral_embedded"
for network in readdir(data_dir)
    if !startswith(network, ".")
        for feeder in readdir(data_dir*"/"*network)
            if !startswith(feeder, ".")
                filename = data_dir*"/"*network*"/"*feeder*"/Master.dss"
                data = get_data_math(file_name)
                results_rosetta_acr = RPMD.solve_opf_acr(data, ipopt_solver; verbose=false);
                results_rosetta_acp = RPMD.solve_opf_acp(data, ipopt_solver; verbose=false);
                if results_rosetta_acr["feasible"]
                    rosetta_acr[filename] = results_rosetta_acr["solution"];
                else
                    rosetta_acr[filename] = results_rosetta_acr["feasible"];
                end
                if results_rosetta_acp["feasible"]
                    rosetta_acp[filename] = results_rosetta_acp["solution"];
                else
                    rosetta_acp[filename] = results_rosetta_acp["feasible"];
                end
            end
        end
    end
end

data = get_data_math(file_name)

## run rosetta PMD ACP and ACR
# if isinteractive() == false
results_rosetta_acr = RPMD.solve_opf_acr(data, ipopt_solver)
# # end
vm_rosetta_acr = results_rosetta_acr["solution"]["vm"]
va_rosetta_acr = results_rosetta_acr["solution"]["va"]
sg_rosetta_acr = results_rosetta_acr["solution"]["sg"]
s_rosetta_acr = results_rosetta_acr["solution"]["s"]


# # if isinteractive() == false
resulta_rosetta_acp = RPMD.solve_opf_acp(data, ipopt_solver)
# # end
vm_rosetta_acp = resulta_rosetta_acp["solution"]["vm"]
va_rosetta_acp = resulta_rosetta_acp["solution"]["va"]
sg_rosetta_acp = resulta_rosetta_acp["solution"]["sg"]
s_rosetta_acp = resulta_rosetta_acp["solution"]["s"]


## run PMD ACP and ACR
# pm_pmd_acr = PMD.instantiate_mc_model(data, PMD.ACRUPowerModel, PMD.build_mc_opf)
# println(pm_pmd_acr.model)
results_pmd_acr = PMD.solve_mc_opf(data, PMD.ACRUPowerModel, ipopt_solver; solution_processors=[PMD.sol_data_model!])
vm_pmd_acr = [bus["vm"] for (i, bus) in results_pmd_acr["solution"]["bus"]]
va_pmd_acr = [bus["va"] for (i, bus) in results_pmd_acr["solution"]["bus"]]
sg_pmd_acr = [gen["pg"]+im*gen["qg"] for (i, gen) in results_pmd_acr["solution"]["gen"]]
s_pmd_acr = [branch["pf"]+im*branch["qf"] for (i, branch) in results_pmd_acr["solution"]["branch"]]
# s_pmd_acr_to = [branch["pt"]+im*branch["qt"] for (i, branch) in results_pmd_acr["solution"]["branch"]]


# pm_pmd_acp = PMD.instantiate_mc_model(data, PMD.ACPUPowerModel, PMD.build_mc_opf)
# println(pm_pmd_acp.model)
results_pmd_acp = PMD.solve_mc_opf(data, PMD.ACPUPowerModel, ipopt_solver; solution_processors=[PMD.sol_data_model!])
vm_pmd_acp = [bus["vm"] for (i, bus) in results_pmd_acp["solution"]["bus"]]
va_pmd_acp = [bus["va"] for (i, bus) in results_pmd_acp["solution"]["bus"]]
sg_pmd_acp = [gen["pg"]+im*gen["qg"] for (i, gen) in results_pmd_acp["solution"]["gen"]]
s_pmd_acp = [branch["pf"]+im*branch["qf"] for (i, branch) in results_pmd_acp["solution"]["branch"]]
