using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using Juniper
using HiGHS
using JuMP  # bl/array_nl
using LinearAlgebra
import LinearAlgebra: diag, diagm
using LaTeXStrings

const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
highs_solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
juniper_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt_solver, "mip_solver" => highs_solver)

""" TODO
- assymetrical alpha's 
- change the testcase
- timeseries analysis
"""

include("./core/solution.jl")

# data_path = "./data/inverter_4w_wye_unbalanced_loads.dss"
data_path = "./data/case5_gen_3ph_wye.dss"


##
data_eng = PMD.parse_file(data_path)
data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)

### add inverter lossy branches
pv_gen_ids = [i for (i, gen) in data_math["gen"] if !occursin("source", gen["name"])]
for gen_id in pv_gen_ids
    RPMD.add_inverter_losses!(data_math, gen_id; multiplexing=true)
end


## 4-leg inverters, set multiplexing true or false, set using binary variables true or false
setting = Dict("multileg" => true, "multiplexing" => true, "multiplexing_binary" => true)
# setting = Dict("multileg" => true, "multiplexing" => false, "multiplexing_binary" => false)

model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, RPMD.build_mc_opf; setting=setting)
PMD.add_start_vrvi!(data_math)

# if setting["multiplexing_binary"]
result = PMD.optimize_model!(model, optimizer=juniper_solver)
# else
# result = PMD.optimize_model!(model, optimizer=ipopt_solver)
# end


##

round.(value.(result["solution"]["gen"]["1"]["bg"]), digits=2)
round.(value.(result["solution"]["gen"]["2"]["bg"]), digits=2)

results_GFL_4w_mx = get_solutions(model, result)

### load, inverters, and source currents
results_GFL_4w_mx["c_load"]
results_GFL_4w_mx["c_inverter_branch"]
results_GFL_4w_mx["c_source_branch"]

### load, inverters, and source currents sequence values in complex
results_GFL_4w_mx["c_load_012"]
results_GFL_4w_mx["c_inverter_branch_012"]
results_GFL_4w_mx["c_source_branch_012"]

### load, inverters, and source currents sequence values in magnitude
round.(abs.(results_GFL_4w_mx["c_load_012"]), digits=6)
round.(abs.(results_GFL_4w_mx["c_inverter_branch_012"]), digits=6)
round.(abs.(results_GFL_4w_mx["c_source_branch_012"]), digits=6)