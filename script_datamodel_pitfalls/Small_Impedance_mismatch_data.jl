import PowerModelsDistribution
import Ipopt
import JuMP
import InfrastructureModels
using rosetta_distribution_opf
import LinearAlgebra: diag, diagm, I
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "sb"=>"yes","warm_start_init_point"=>"yes", "max_iter"=>1000)


function solve_4w_IVR(file_name; symmetry=false)
    data = PMD.parse_file(file_name, transformations=[PMD.remove_all_bounds!, PMD.transform_loops!, PMD.join_lines!])
    RPMD.pv1_correction!(data)
    data["settings"]["sbase_default"] = 1.0
    data["voltage_source"]["source"]["rs"] *= 0
    data["voltage_source"]["source"]["xs"] *= 0
    # for (i, bus) in data["bus"]
    #     line["length"] = 1
    # end
    if PMD.iseng(data)
        data = PMD.transform_data_model(data, multinetwork=false, kron_reduce=false, phase_project=false)
    end
    PMD.add_start_vrvi!(data)

    for i in 2:10
        data["gen"]["$i"] = deepcopy( data["gen"]["1"])
        data["gen"]["$i"]["gen_bus"] = data["load"]["$i"]["load_bus"]
    end
    for (g,gen) in data["gen"]
        if g !== "1"
            gen["pmin"] =  0*ones(3);
            gen["pmax"] =  2*ones(3);
            gen["qmin"] = -0*ones(3);
            gen["qmax"] =  2*ones(3);
            gen["cost"] *= 10
        else
            gen["cost"] *= 1000
        end
    end

    if symmetry
        for (b, branch) in data["branch"]
            diags = diag(branch["br_r"])
            off_diags = branch["br_r"] .- diagm(diag(branch["br_r"]))
            branch["br_r"] = sum(diags)/4 * I(4) + zeros(4,4) #.+ sum(off_diags)/12 * (ones(4).*ones(4)' - I(4))

            diags = diag(branch["br_x"])
            off_diags = branch["br_x"] .- diagm(diag(branch["br_x"]))
            branch["br_x"] = sum(diags)/4 * I(4)  + zeros(4,4) #.+ 0 sum(off_diags)/12 * (ones(4).*ones(4)' - I(4))
        end
    end
    
    for (i, load) in data["load"]
        load["pd"] .*= 0.01
        load["qd"] .*= 0.01
    end

    result = PMD.solve_mc_opf(data, PMD.IVRENPowerModel, ipopt_solver)

    return data, result
end

##
file_name_normal = "./data/Network1_Feeder1_Impedance_sections_Normal/Master.dss"


# file_name_normal = "./data/Network1_Feeder1_Impedance_sections_Normal_mismatch/Master.dss"

data, results_normal = solve_4w_IVR(file_name_normal)

data, results_mismatch = solve_4w_IVR(file_name_normal; symmetry=true)

#
sum([load["pd"] for (i,load) in data["load"]])

[data["load"]["$i"]["load_bus"] for i in collect(1:11)]



