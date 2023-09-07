using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

# data_path = "./data/test_gen_3ph_wye.dss"
# data_path = "./data/test_gen_3ph_delta.dss"
# data_path = "./data/test_gen_1ph_wye.dss"
# data_path = "./data/test_gen_1ph_delta.dss"
# data_path = "./data/test_load_1ph_delta_cp.dss"
data_path = "./data/test_load_1ph_wye_cp.dss"
# data_path = "./data/test_load_3ph_delta_cp.dss"
# data_path = "./data/test_load_3ph_wye_cp.dss"
# data_path = "./data/test_load_3ph_delta_ci.dss"
# data_path = "./data/test_load_3ph_wye_ci.dss"
# data_path = "./data/test_load_3ph_delta_cz.dss"
# data_path = "./data/test_load_3ph_wye_cz.dss"
# data_path = "./data/test_trans_dy.dss"  # grounding reactor is removed + initialization should be updated
# data_path = "./data/test_trans_yy.dss"  # grounding reactor is removed + initialization should be updated

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!])

RPMD.pv1_correction!(data_eng)

data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)
# PMD.add_start_vrvi!(data_math)

ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]


##
model = JuMP.Model(Ipopt.Optimizer)
RPMD.IVR_EN(model, ref)

###
JuMP.optimize!(model)
feasible = (JuMP.termination_status(model) == JuMP.LOCALLY_SOLVED)
cost = JuMP.objective_value(model)

##
PMD.add_start_vrvi!(data_math)
pm = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, PMD.build_mc_opf);
res = PMD.optimize_model!(pm, optimizer=ipopt_solver)
println(pm.model)