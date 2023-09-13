import PowerModelsDistribution
import Ipopt
import JuMP
import InfrastructureModels
using rosetta_distribution_opf
import LinearAlgebra: diag
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "sb"=>"yes","warm_start_init_point"=>"yes", "max_iter"=>1000)

file_name = "./data/Network1_Feeder1_Impedance_sections_Normal/Master.dss"

function solve_4w_IVR(file_name; Zsmall_section=false)
    data = PMD.parse_file(file_name, transformations=[PMD.remove_all_bounds!, PMD.transform_loops!])
    RPMD.pv1_correction!(data)
    data["settings"]["sbase_default"] = 1.0
    data["voltage_source"]["source"]["rs"] *= 0
    data["voltage_source"]["source"]["xs"] *= 0
    if PMD.iseng(data)
        data = PMD.transform_data_model(data, multinetwork=false, kron_reduce=false, phase_project=false)
    end
    PMD.add_start_vrvi!(data)

    for i in 2:10
        data["gen"]["$i"] = deepcopy( data["gen"]["1"])
        data["gen"]["$i"]["gen_bus"] = data["load"]["$i"]["load_bus"]
    end
    for (g,gen) in data["gen"]
        gen["pmin"] =   0*ones(3);
        gen["pmax"] =  20*ones(3);
        gen["qmin"] = -20*ones(3);
        gen["qmax"] =  20*ones(3);
        gen["cost"] *= 1000
    end

    if Zsmall_section
        bus1 = data["gen"]["2"]["gen_bus"]
        bus2 = data["gen"]["3"]["gen_bus"]
        branches = []
        bus1 = [bus1]
        iter = 1
        while bus2 ∉ bus1 && iter < 400
            branch_ids = [b for (b, branch) in data["branch"] if (branch["t_bus"] in bus1) || (branch["f_bus"] in bus1)]
            bus1 = []
            for id in branch_ids
                append!(bus1, data["branch"][id]["f_bus"]∈ bus1 ? data["branch"][id]["t_bus"] : data["branch"][id]["f_bus"])
            end
            append!(branches, branch_ids)
            iter += 1
        end
        @show length(unique(branches))

        for i in unique(branches)
            data["branch"][i]["br_r"] *= 1E-6
            data["branch"][i]["br_x"] *= 1E-6
        end
    end

    return PMD.solve_mc_opf(data, PMD.IVRENPowerModel, ipopt_solver)
end

##

file_name = "./data/Network1_Feeder1_Impedance_sections_Normal/Master.dss"
results_Zsmall = solve_4w_IVR(file_name; Zsmall_section=true)

results = solve_4w_IVR(file_name)


# file_name_small = "./data/Network1_Feeder1_Impedance_sections_Small/Master.dss"
# results_small = solve_4w_IVR(file_name_small)
