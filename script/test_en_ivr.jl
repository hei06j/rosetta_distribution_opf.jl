using Pkg
Pkg.activate("./")
using jump_pmd_ivr
import PowerModelsDistribution
using Ipopt
using JuMP
const _PMD = PowerModelsDistribution
const _RPMD = jump_pmd_ivr

data_path = "./data/test_gen_3ph_wye.dss"

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
data_eng = _PMD.parse_file(data_path, transformations=[_PMD.remove_all_bounds!])

"The solar parsing is a bit off, this method corrects for that."
function pv1_correction!(data_eng)
    pv1 = data_eng["solar"]["pv1"]
    pv1["pg_lb"] = pv1["pg_ub"] = pv1["pg"]
    pv1["qg_lb"] = pv1["qg_ub"] = pv1["qg"]
end

pv1_correction!(data_eng)

data_math = _PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)
_PMD.add_start_vrvi!(data_math)


##
pm = _PMD.instantiate_mc_model(data_math, _PMD.IVRENPowerModel, _RPMD.build_mc_opf);
res = _PMD.optimize_model!(pm, optimizer=ipopt_solver)