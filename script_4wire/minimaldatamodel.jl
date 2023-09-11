using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP  # bl/array_nl
using JSON
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
    PMD.print_file(data_path[1:end-4]*".json", de)
    return de
end

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!])
de = clean_export_JSON(data_path)

RPMD.pv1_correction!(data_eng)
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)
