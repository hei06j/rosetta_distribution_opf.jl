using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import OpenDSSDirect
import InfrastructureModels
using Ipopt
using JuMP  # bl/array_nl
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels
const ODSS = OpenDSSDirect


##
ODSS.dss("Redirect ./data/test_gen_3ph_wye.dss")
vm_noControl = ODSS.Circuit.AllBusVolts()
cm_noControl = ODSS.PDElements.AllCurrentsAllCurrents()[5:8]
vm1_noControl = abs.(vm_noControl)[1:3]
vm2_noControl = abs.(vm_noControl)[4:7]
pg_noControl = ODSS.PVsystems.kW()
qg_noControl = ODSS.PVsystems.kvar()

ODSS.dss("Redirect ./data/test_gen_3w_wye_InvControl_v2.dss")
vm_Control = ODSS.Circuit.AllBusVolts()
cm_Control = ODSS.PDElements.AllCurrentsAllCurrents()[4:6]
vm1_Control = abs.(vm_Control)[1:3]
vm2_Control = abs.(vm_Control)[4:6]
pg_Control = ODSS.PVsystems.kW()
qg_Control = ODSS.PVsystems.kvar()

# plot([vm1_noControl[1:3] vm2_noControl[1:3]]', marker=:circle, linestyle=:dash, color=[:blue :red :green], label=["v_a" "v_b" "v_c"])
# plot!([vm1_Control[1:3] vm2_Control[1:3]]', marker=:circle, linestyle=:solid, color=[:blue :red :green], label=["v_a^c" "v_b^c" "v_c^c"])

[pg_noControl pg_Control]
[qg_noControl qg_Control]
qg_Control
round.(vm2_Control./230.94, digits=4)
# qmax = sqrt.((25)^2 - (23)^2)/3 * ones(3)
# # sum(-2*qmax/0.2 .* ((vm2_noControl[1:3] .-  vm2_noControl[4])/230.94.- 1))
# sum(-2*qmax/0.2 .* ((vm2_Control[1:3])/230.94.- 1))

##
data_path = "./data/inverter_3w_wye_unbalanced_loads_v2.dss"

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!, PMD.transform_loops!])
RPMD.pv1_correction!(data_eng)
data_eng["settings"]["sbase_default"] = 1
vbase = data_eng["settings"]["vbases_default"]["b1"]
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)

data_math["bus"]["1"]["vmin"] = 0.9 * ones(3)
data_math["bus"]["1"]["vmax"] = 1.1 * ones(3)

gen_id = 1
gen = data_math["gen"]["$gen_id"]
gen["pmax"] = 23/3 * ones(3)
gen["pmin"] = 23/3 * ones(3)
gen["pg"] = 23/3 * ones(3)
gen["qmax"] = sqrt.(25^2 - 23^2)/3 * ones(3)
gen["qmin"] = -gen["qmax"]

gen["connections"] = gen["connections"][1:3]
data_math["bus"]["1"]["terminals"] = data_math["bus"]["1"]["terminals"][1:3]
data_math["bus"]["1"]["grounded"] = data_math["bus"]["1"]["grounded"][1:3]


gen["cost"] = [0 0]
data_math["gen"]["2"]["cost"] = [1000 0]
##
ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]
objective = "cost"

model = JuMP.Model(Ipopt.Optimizer)
include("./core/variables.jl")
include("./core/constraints.jl")
include("./core/objectives.jl")

## Inverter control Volt-var
gen_id = 1
gen = data_math["gen"]["$gen_id"]
bus_id = gen["gen_bus"]
vmin = 0.9; vmax = 1.1;
phases = gen["connections"]
# n = gen["connections"][end]

terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref[:bus])
vm = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm", lower_bound = 0 ) for i in keys(ref[:bus]))
vm = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm[i].axes[1] ? vm[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
for (i, bus) in ref[:bus]
    JuMP.@constraint(model, [p in ref[:bus][i]["terminals"]], vm[p,i]^2 == vr[p,i]^2 + vi[p,i]^2)
end

# vm_square = JuMP.@expression(model, vr[phases,bus_id].^2 .+ vi[phases,bus_id].^2)
# # vm_pv = JuMP.@expression(model, sqrt.((vr[phases,bus_id].-vr[n,bus_id]).^2 .+ (vi[phases,bus_id].-vi[n,bus_id]).^2) )
# vm_pv = JuMP.@expression(model, sqrt.((vr[phases,bus_id]).^2 .+ (vi[phases,bus_id]).^2) )
# # JuMP.@constraint(model, qg[phases,gen_id] .== -2*gen["qmax"]/0.2 .* (vm_pv[phases,bus_id] .- 1))

vv_curve(v, qmax) = v <= 0.95 ? qmax : (v >= 1.05 ? -qmax : -2*qmax/0.1*(v-1))
JuMP.@operator(model, vv, 2, vv_curve)

# JuMP.@constraint(model, [i in phases], qg[i,gen_id] == vv(vm_pv[i,bus_id], gen["qmax"][i]) )
JuMP.@constraint(model, qg[1, gen_id] == vv(sum(vm[i,bus_id] for i in ref[:bus][1]["terminals"][1:3])/3, gen["qmax"][1]) )
JuMP.@constraint(model, qg[2, gen_id] == qg[1, gen_id])
JuMP.@constraint(model, qg[3, gen_id] == qg[1, gen_id])

##
JuMP.optimize!(model)
cost = JuMP.objective_value(model)
v = abs.(JuMP.value.(vm)) * 230.94
pg_values = JuMP.value.(pg)
qg_values = JuMP.value.(qg)
vbase
sum(pg_values[:,gen_id])
sum(qg_values[:,gen_id])
round.(abs.(JuMP.value.(vm))[:,1], digits=4)


# vm_range = 0.8:0.01:1.2
# plot(vm_range, vv_curve.(vm_range,gen["qmax"][1]), label="VV curve")
# # plot!(vm_range, vv_curve.(vm_range,gen["qmax"][2]))
# # plot!(vm_range, vv_curve.(vm_range,gen["qmax"][3]))
# plot!([sum(v[phases,bus_id])/vbase/1000/3], [qg_values[1,gen_id]], seriestype=:scatter, label="VV point")

# vm_range = 0.8:0.01:1.2
# vw_curve(v, pmax) = v <= 0.95 ? pmax : (v >= 1.05 ? 0 : -pmax/0.1*(v-1.05))
# plot(vm_range, vw_curve.(vm_range,3), label="VW curve")

