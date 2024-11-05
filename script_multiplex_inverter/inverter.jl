using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP
import LinearAlgebra: diag, diagm
using LaTeXStrings

const _PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")

##
# data_path = "./data/inverter_3w_wye_unbalanced_loads.dss"
# data_eng = _PMD.parse_file(data_path, transformations=[_PMD.transform_loops!])
# data_eng["settings"]["sbase_default"] = 1
# data_eng["voltage_source"]["source"]["rs"] *= 0
# data_eng["voltage_source"]["source"]["xs"] *= 0
# data_math = _PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)
# # RPMD.set_gen_max_powers!(data_math, "1", 10)
# _PMD.add_start_vrvi!(data_math)
# result = _PMD.solve_mc_opf(data_math, _PMD.IVRENPowerModel, ipopt_solver)

# ref_branch = 1
# c_ref = result["solution"]["branch"]["$ref_branch"]["cr_fr"] .+ im * result["solution"]["branch"]["$ref_branch"]["ci_fr"]
# -sum(c_ref[1:3])

# c_ref_gen = result["solution"]["gen"]["2"]["crg"] .+ im * result["solution"]["gen"]["2"]["cig"]

# c_inv = result["solution"]["gen"]["1"]["crg"] .+ im * result["solution"]["gen"]["1"]["cig"]
# c_INV = [c_inv[1]-c_inv[3] ; c_inv[2]-c_inv[1] ; c_inv[3]-c_inv[2]]
# sum(c_INV)

# c_ref_to = result["solution"]["branch"]["$ref_branch"]["cr_to"] .+ im * result["solution"]["branch"]["$ref_branch"]["ci_to"]



## 3-leg inverters
data_path = "./data/inverter_3w_wye_unbalanced_loads.dss"

### 4-leg inverters, set multiplexing true or false, set using binary variables true or false
setting = Dict("multileg" => true, "multiplexing" => false, "multiplexing_binary" => false, "ideal" => false, "dc_link" => true)

### parse data
data_eng = _PMD.parse_file(data_path, transformations=[_PMD.transform_loops!])
data_eng["settings"]["sbase_default"] = 1
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = _PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)

### add inverter lossy branches
pv_gen_ids = [i for (i, gen) in data_math["gen"] if !occursin("source", gen["name"])]
for gen_id in pv_gen_ids
    RPMD.add_inverter_losses!(data_math, gen_id; three_wire=true, multiplexing=setting["multiplexing"], dc_link=setting["dc_link"])
end

### solve opf
_PMD.add_start_vrvi!(data_math)
model = _PMD.instantiate_mc_model(data_math, _PMD.IVRENPowerModel, RPMD.build_mc_opf_mx; setting=setting)
result = _PMD.optimize_model!(model, optimizer=ipopt_solver)
results_GFL_3w = RPMD.get_solutions(model, result)

### load, inverters, and source currents
c_load_3w = length(results_GFL_3w["c_load"])==4 ? results_GFL_3w["c_load"] : results_GFL_3w["c_load"][1:4].+results_GFL_3w["c_load"][5:8]
results_GFL_3w["c_inverter_branch"]
results_GFL_3w["c_source_branch"]

[abs.(c_load_3w) angle.(c_load_3w).*180/pi]
[abs.(results_GFL_3w["c_inverter_branch"]) angle.(results_GFL_3w["c_inverter_branch"]).*180/pi]
[abs.(results_GFL_3w["c_source_branch"]) angle.(results_GFL_3w["c_source_branch"]).*180/pi]

### load, inverters, and source currents sequence values in complex
c_load_3w_012 = length(results_GFL_3w["c_load_012"])==3 ? results_GFL_3w["c_load_012"] : results_GFL_3w["c_load_012"][1:3].+results_GFL_3w["c_load_012"][4:6]
results_GFL_3w["c_inverter_branch_012"]
results_GFL_3w["c_source_branch_012"]

### load, inverters, and source currents sequence values in magnitude
round.(abs.(results_GFL_3w["c_load_012"]), digits=6)
round.(abs.(results_GFL_3w["c_inverter_branch_012"]), digits=6)
round.(abs.(results_GFL_3w["c_source_branch_012"]), digits=6)

result["solution"]["gen"]["1"]

# ref_branch = 1
# c_ref = result["solution"]["branch"]["$ref_branch"]["cr_fr"] .+ im * result["solution"]["branch"]["$ref_branch"]["ci_fr"]
# -sum(c_ref[1:3])

# c_ref_gen = result["solution"]["gen"]["2"]["crg"] .+ im * result["solution"]["gen"]["2"]["cig"]
# sum(c_ref_gen)

# c_inv = result["solution"]["gen"]["1"]["crg"] .+ im * result["solution"]["gen"]["1"]["cig"]
# c_INV = [c_inv[1]-c_inv[3] ; c_inv[2]-c_inv[1] ; c_inv[3]-c_inv[2]]
# sum(c_INV)

# c_ref_to = result["solution"]["branch"]["$ref_branch"]["cr_to"] .+ im * result["solution"]["branch"]["$ref_branch"]["ci_to"]




## 4-leg inverters
data_path = "./data/inverter_4w_wye_unbalanced_loads.dss"

### 4-leg inverters, set multiplexing true or false, set using binary variables true or false
setting = Dict("multileg" => true, "multiplexing" => false, "multiplexing_binary" => false, "ideal" => false, "dc_link" => true)

### parse data
data_eng = _PMD.parse_file(data_path, transformations=[_PMD.transform_loops!])
data_eng["settings"]["sbase_default"] = 1
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = _PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)

### add inverter lossy branches
pv_gen_ids = [i for (i, gen) in data_math["gen"] if !occursin("source", gen["name"])]
for gen_id in pv_gen_ids
    RPMD.add_inverter_losses!(data_math, gen_id; three_wire=false, multiplexing=setting["multiplexing"])
end

### solve opf
_PMD.add_start_vrvi!(data_math)
model = _PMD.instantiate_mc_model(data_math, _PMD.IVRENPowerModel, RPMD.build_mc_opf_mx; setting=setting)
result = _PMD.optimize_model!(model, optimizer=ipopt_solver)
results_GFL_4w = RPMD.get_solutions(model, result)

### load, inverters, and source currents
c_load_4w = length(results_GFL_4w["c_load"])==4 ? results_GFL_4w["c_load"] : results_GFL_4w["c_load"][1:4].+results_GFL_4w["c_load"][5:8]
results_GFL_4w["c_inverter_branch"]
results_GFL_4w["c_source_branch"]

[abs.(c_load_4w) angle.(c_load_4w).*180/pi]
[abs.(results_GFL_4w["c_inverter_branch"]) angle.(results_GFL_4w["c_inverter_branch"]).*180/pi]
[abs.(results_GFL_4w["c_source_branch"]) angle.(results_GFL_4w["c_source_branch"]).*180/pi]

### load, inverters, and source currents sequence values in complex
c_load_4w_012 = length(results_GFL_4w["c_load_012"])==3 ? results_GFL_4w["c_load_012"] : results_GFL_4w["c_load_012"][1:3].+results_GFL_4w["c_load_012"][4:6]
results_GFL_4w["c_inverter_branch_012"]
results_GFL_4w["c_source_branch_012"]

### load, inverters, and source currents sequence values in magnitude
round.(abs.(results_GFL_4w["c_load_012"]), digits=6)
round.(abs.(results_GFL_4w["c_inverter_branch_012"]), digits=6)
round.(abs.(results_GFL_4w["c_source_branch_012"]), digits=6)


##
using Plots
mkpath("./Figures")

##
# Imax_3w = maximum([abs.(results_GFL_3w["c_source_branch"])  abs.(c_load_3w) abs.(results_GFL_3w["c_inverter_branch"])])
# Imax_4w = maximum([abs.(results_GFL_4w["c_source_branch"])  abs.(c_load_4w) abs.(results_GFL_4w["c_inverter_branch"])])
# Imax = maximum([Imax_3w, Imax_4w])
Imax = maximum([abs.(c_load_3w) abs.(c_load_4w)] )

I0_3w = round(abs(results_GFL_3w["c_source_branch_012"][1]), digits=2)
I2_3w = round(abs(results_GFL_3w["c_source_branch_012"][3]), digits=2)
GFL_3w_vuf_c = RPMD.plot_phasors(results_GFL_3w["c_source_branch"], Imax; I2=I2_3w, I0=I0_3w)
Plots.savefig(GFL_3w_vuf_c, "./Figures/GFL_3w_vuf_c.pdf")
GFL_3w_vuf_cg = RPMD.plot_phasors(results_GFL_3w["c_inverter_branch"], Imax)
Plots.savefig(GFL_3w_vuf_cg, "./Figures/GFL_3w_vuf_cg.pdf")
GFL_3w_vuf_cd = RPMD.plot_phasors(c_load_3w, Imax)
Plots.savefig(GFL_3w_vuf_cd, "./Figures/GFL_3w_vuf_cd.pdf")
GFL_3w_vuf = Plots.plot(GFL_3w_vuf_c, GFL_3w_vuf_cg, GFL_3w_vuf_cd, layout=(1,3), size=(2200,700))
Plots.savefig(GFL_3w_vuf, "./Figures/GFL_3w_vuf.pdf")
Plots.savefig(GFL_3w_vuf, "./Figures/GFL_3w_vuf.png")


I0_4w = round(abs(results_GFL_4w["c_source_branch_012"][1]), digits=2)
I2_4w = round(abs(results_GFL_4w["c_source_branch_012"][3]), digits=2)
GFL_4w_vuf_c = plot_phasors(results_GFL_4w["c_source_branch"], Imax; I2=I2_4w, I0=I0_4w)
Plots.savefig(GFL_4w_vuf_c, "./Figures/GFL_4w_vuf_c.pdf")
GFL_4w_vuf_cg = plot_phasors(results_GFL_4w["c_inverter_branch"], Imax)
Plots.savefig(GFL_4w_vuf_cg, "./Figures/GFL_4w_vuf_cg.pdf")
GFL_4w_vuf_cd = plot_phasors(c_load_4w, Imax, labeled=true)
Plots.savefig(GFL_4w_vuf_cd, "./Figures/GFL_4w_vuf_cd.pdf")
GFL_4w_vuf = Plots.plot(GFL_4w_vuf_c, GFL_4w_vuf_cg, GFL_4w_vuf_cd, layout=(1,3), size=(2200,700))
Plots.savefig(GFL_4w_vuf, "./Figures/GFL_4w_vuf.pdf")
Plots.savefig(GFL_4w_vuf, "./Figures/GFL_4w_vuf.png")


GFL_4w_3w_vuf = plot(GFL_3w_vuf, GFL_4w_vuf, layout=(2,1), size=(2200,1400))
Plots.savefig(GFL_4w_3w_vuf, "./Figures/GFL_4w_3w_vuf.pdf")

