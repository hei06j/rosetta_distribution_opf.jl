using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP  # bl/array_nl
using JSON
using Statistics
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

data_path = "./data/test_gen_3ph_wye.dss"
# data_path = "./data/test_gen_3ph_delta.dss"
# data_path = "./data/test_gen_1ph_wye.dss"
# data_path = "./data/test_gen_1ph_delta.dss"
# data_path = "./data/test_load_1ph_delta_cp.dss"
# data_path = "./data/test_load_1ph_wye_cp.dss"
# data_path = "./data/test_load_3ph_delta_cp.dss"
# data_path = "./data/test_load_3ph_wye_cp.dss"
# data_path = "./data/test_load_3ph_delta_ci.dss"
# data_path = "./data/test_load_3ph_wye_ci.dss"
# data_path = "./data/test_load_3ph_delta_cz.dss"
# data_path = "./data/test_load_3ph_wye_cz.dss"

function clean_export_JSON(data_path)
    de = PMD.parse_file(data_path)

    delete!(de, "data_path")
    delete!(de, "settings")
    delete!(de, "files")
    delete!(de, "data_model")
    delete!(de, "conductor_ids")
    delete!(de, "solar")

    for (i,bus) in de["bus"]
        delete!(bus,"grounded")
        delete!(bus,"rg")
        delete!(bus,"xg")
        delete!(bus,"status")
        bus["vmin"] = [0.9*ones(3); 0]
        bus["vmax"] = [1.1*ones(3); 0.1]
        bus["vpnmin"] = 0.9*ones(3)
        bus["vpnmax"] = 1.1*ones(3)
    end

    for (l,line) in de["line"]
        delete!(line,"status")
        delete!(line,"source_id")
        delete!(line,"f_connections")
        delete!(line,"t_connections")
    end

    for (s, source) in de["voltage_source"]
        delete!(source,"rs")
        delete!(source,"xs")
        delete!(source,"status")
        delete!(source,"configuration")
        delete!(source,"source_id")
        delete!(source,"connections")
    end 

    for (d, load) in de["load"]
        delete!(load,"source_id")
        delete!(load,"status")
        delete!(load,"model")
        delete!(load,"connections")
        delete!(load,"vm_nom")
        delete!(load,"dispatchable")
        delete!(load,"configuration")
    end 


    # for (d, load) in de["gen"]
    #     delete!(load,"source_id")
    #     delete!(load,"status")
    #     delete!(load,"model")
    #     delete!(load,"connections")
    #     delete!(load,"vm_nom")
    #     delete!(load,"dispatchable")
    #     delete!(load,"configuration")
    # end 

    PMD.print_file(data_path[1:end-4]*".json", de)
    return de
end
de = clean_export_JSON(data_path)


function parse_reduced_JSON_to_PMD(red_eng)
    de = deepcopy(red_eng)
    de["settings"] = Dict()
    de["settings"]["sbase_default"] = 1000
    de["settings"]["voltage_scale_factor"] = 1000
    de["settings"]["power_scale_factor"] = 1000
    de["settings"]["base_frequency"] = 60
    de["settings"]["vbases_default"] = Dict()
    
    de["data_model"] = PMD.ENGINEERING
    de["conductor_ids"] = [1;2;3;4]
    de["solar"] = Dict()

    for (i,bus) in de["bus"]
        bus["grounded"] = []
        bus["rg"] = []
        bus["xg"] =[]
        bus["status"] = PMD.ENABLED
        bus["terminals"] = [1;2;3;4]
    end

    for (l,line) in de["line"]
        line["status"] = PMD.ENABLED
        line["f_connections"] = [1;2;3;4]
        line["t_connections"] = [1;2;3;4]
    end

    for (s, source) in de["voltage_source"]
        source["rs"] = zeros(Float64, 4, 4)
        source["xs"] = zeros(Float64, 4, 4)
        source["status"] = PMD.ENABLED
        source["configuration"] = PMD.WYE
        source["connections"] = [1;2;3;4]
        de["settings"]["vbases_default"][source["bus"]] = mean(source["vm"])
    end 

    for (d, load) in de["load"]
        load["status"]  = PMD.ENABLED
        load["model"] = PMD.POWER
        load["connections"] = [1;2;3;4]
        load["vm_nom"] = 1
        load["dispatchable"] = PMD.NO
        load["configuration"] = PMD.WYE
    end 
    return de
end

de2 = parse_reduced_JSON_to_PMD(de)



data_math = PMD.transform_data_model(de2, multinetwork=false, kron_reduce=false, phase_project=false)
PMD.add_start_vrvi!(data_math)
pm = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, PMD.build_mc_opf);
res = PMD.optimize_model!(pm, optimizer=ipopt_solver)



