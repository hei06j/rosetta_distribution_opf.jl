using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
import OpenDSSDirect
using Ipopt
using JuMP  # bl/array_nl
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels
const ODSS = OpenDSSDirect


##
# ODSS.dss("Redirect ./data/test_gen_4w_wye.dss")
# vm_noControl = ODSS.Circuit.AllBusVolts()
# cm_noControl = ODSS.PDElements.AllCurrentsAllCurrents()[5:8]
# vm1_noControl = abs.(vm_noControl)[1:3]
# vm2_noControl = abs.(vm_noControl)[4:7]
# pg_noControl = ODSS.PVsystems.kW()
# qg_noControl = ODSS.PVsystems.kvar()

ODSS.dss("Redirect ./data/test_gen_4w_wye_InvControl_p2g.dss")
vm_Control_p2g = ODSS.Circuit.AllBusVolts()
cm_Control_p2g = ODSS.PDElements.AllCurrentsAllCurrents()[5:8]
vm1_Control_p2g = abs.(vm_Control_p2g)[1:3]
vm2_Control_p2g = abs.(vm_Control_p2g)[4:7]./230.94
pg_Control_p2g = ODSS.PVsystems.kW()
qg_Control_p2g = ODSS.PVsystems.kvar()

# s_load = [9+im*4.358898943540673 ; 4.5+im*2.1794494717703365 ; 0]*1E3
# v_load = vm_Control_p2g[8:10]
# c_load_transpose = s_load./v_load
# s_line21 = (vm_Control_p2n[4:6].-vm_Control_p2n[7]) .* transpose(cm_Control_p2n[5:7])'
# s_g2 = (vm_Control_p2n[4:6].-vm_Control_p2n[7]) .* (transpose(cm_Control_p2n[5:7])' .+ c_load_transpose)
# pg3 = real.(s_g2) ./ 1E3
# qg3 = imag.(s_g2) ./ 1E3
# sum(pg3)
# sum(qg3)
# round.(vm2_Control_p2g, digits=4)



##
data_path = "./data/test_gen_4w_wye_p2g.dss"

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!, PMD.transform_loops!])
RPMD.pv1_correction!(data_eng)
data_eng["settings"]["sbase_default"] = 1
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)

data_math["bus"]["1"]["vmin"] = [0.9 * ones(3) ; 0 ]
data_math["bus"]["1"]["vmax"] = [1.1 * ones(3) ; Inf]

gen_id = 1
gen = data_math["gen"]["$gen_id"]
smax = 40
pmax = 35
gen["pmax"] = [pmax/3 * ones(3) ; 0]
gen["pmin"] = [pmax/3 * ones(3) ; 0]
gen["pg"] = [pmax/3 * ones(3) ; 0]
gen["qmax"] = sqrt.(smax^2 - pmax^2)/3 * ones(4)
gen["qmin"] = -gen["qmax"]

gen["cost"] = [0 0]
data_math["gen"]["2"]["cost"] = [1000 0]
###
ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]

model = JuMP.Model(Ipopt.Optimizer)

objective = "cost"
include("./core/variables.jl")
include("./core/constraints.jl")
include("./objectives.jl")


### Inverter control Volt-var
gen_id = 1
gen = data_math["gen"]["$gen_id"]
bus_id = gen["gen_bus"]
vmin = 0.9; vmax = 1.1;
phases = gen["connections"][1:end-1]
n = gen["connections"][end]

terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref[:bus])
vm = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm", lower_bound = 0 ) for i in keys(ref[:bus]))
vm = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm[i].axes[1] ? vm[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
for (i, bus) in ref[:bus]
    JuMP.@constraint(model, [p in ref[:bus][i]["terminals"]], vm[p,i]^2 == vr[p,i]^2 + vi[p,i]^2)
end

vmn = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vmn", lower_bound = 0 ) for i in keys(ref[:bus]))
vmn = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vmn[i].axes[1] ? vmn[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
for (i, bus) in ref[:bus]
    JuMP.@constraint(model, [p in ref[:bus][i]["terminals"]], vmn[p,i]^2 == (vr[p,i]-vr[4,i])^2 + (vi[p,i]-vi[4,i])^2)
end

# vm_square = JuMP.@expression(model, vr[phases,bus_id].^2 .+ vi[phases,bus_id].^2)
# # vm_pv = JuMP.@expression(model, sqrt.((vr[phases,bus_id].-vr[n,bus_id]).^2 .+ (vi[phases,bus_id].-vi[n,bus_id]).^2) )
# vm_pv = JuMP.@expression(model, sqrt.((vr[phases,bus_id]).^2 .+ (vi[phases,bus_id]).^2) )
# # JuMP.@constraint(model, qg[phases,gen_id] .== -2*gen["qmax"]/0.2 .* (vm_pv[phases,bus_id] .- 1))

vv_curve(v, qmax) = v <= 0.95 ? qmax : (v >= 1.05 ? -qmax : -2*qmax/0.1*(v-1))
JuMP.@operator(model, vv, 2, vv_curve)

# JuMP.@constraint(model, [i in phases], qg[i,gen_id] == vv(vm_pv[i,bus_id], gen["qmax"][i]) )
JuMP.@constraint(model, qg[1, gen_id] == vv(sum(vm[i,bus_id] for i in ref[:bus][1]["terminals"][1:3])/3, gen["qmax"][1]) )
# JuMP.@constraint(model, qg[1, gen_id] == vv(sum(vmn[i,bus_id] for i in ref[:bus][1]["terminals"][1:3])/3, gen["qmax"][1]) )
JuMP.@constraint(model, qg[2, gen_id] == qg[1, gen_id])
JuMP.@constraint(model, qg[3, gen_id] == qg[1, gen_id])

###
JuMP.optimize!(model)
cost = JuMP.objective_value(model)
v = abs.(JuMP.value.(vm)) * 230.94
pg_values = JuMP.value.(pg)
qg_values = JuMP.value.(qg)
vbase = 230.94
sum(pg_values[:,gen_id])
sum(qg_values[:,gen_id])
round.(abs.(JuMP.value.(vm))[:,1], digits=4)

##
@show [round.(vm2_Control_p2g, digits=4) ; qg_Control_p2g]

@show [round.(abs.(JuMP.value.(vm))[:,1], digits=4); sum(qg_values[:,gen_id])]
