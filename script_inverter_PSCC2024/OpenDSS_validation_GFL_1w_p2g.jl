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
ODSS.dss("Redirect ./data/test_gen_1w_wye_InvControl_p2g.dss")
vm_Control_p2g = ODSS.Circuit.AllBusVolts()
cm_Control_p2g = ODSS.PDElements.AllCurrentsAllCurrents()
vm1_Control_p2g = abs.(vm_Control_p2g)[1:3]
vm2_Control_p2g = abs.(vm_Control_p2g)[4:7]./230.94
pg_Control_p2g = ODSS.PVsystems.kW()
qg_Control_p2g = ODSS.PVsystems.kvar()

round.(vm2_Control_p2g, digits=4)

ODSS.PVsystems.First()
pg1 = ODSS.PVsystems.kW()
qg1 = ODSS.PVsystems.kvar()

ODSS.PVsystems.Next()
pg2 = ODSS.PVsystems.kW()
qg2 = ODSS.PVsystems.kvar()

ODSS.PVsystems.Next()
pg3 = ODSS.PVsystems.kW()
qg3 = ODSS.PVsystems.kvar()

qg_dss = [qg1;qg2;qg3]

@show [round.(vm2_Control_p2g, digits=4) ; round.(qg_dss, digits=4) ; round.(sum(qg_dss), digits=4) ]

##
data_path = "./data/test_gen_1w_wye_p2g.dss"

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!, PMD.transform_loops!])
data_eng["settings"]["sbase_default"] = 1
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)

data_math["bus"]["1"]["vmin"] = [0.7 * ones(3) ; 0 ]
data_math["bus"]["1"]["vmax"] = [1.3 * ones(3) ; Inf]

for gen_id in [1, 2, 3]
    gen = data_math["gen"]["$gen_id"]
    smax = 13.33
    pmax = 11.66
    gen["pmax"] = pmax
    gen["pmin"] = pmax
    gen["pg"] = pmax
    gen["qmax"] = sqrt.(smax^2 - pmax^2)
    gen["qmin"] = -gen["qmax"]
    gen["cost"] = [0 0]
end
data_math["gen"]["2"]["cost"] = [1000 0]


ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]
ref[:bus][1]["grounded"][4] = 1
###
model = JuMP.Model(Ipopt.Optimizer)


objective = "cost"
include("./core/variables.jl")
include("./core/constraints.jl")
include("./core/objectives.jl")


### Inverter control Volt-var
terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref[:bus])
vm = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm", lower_bound = 0 ) for i in keys(ref[:bus]))
vm = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm[i].axes[1] ? vm[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
for (i, bus) in ref[:bus]
    JuMP.@constraint(model, [p in ref[:bus][i]["terminals"]], vm[p,i]^2 == vr[p,i]^2 + vi[p,i]^2)
end

vv_curve(v, qmax) = v <= 0.95 ? qmax : (v >= 1.05 ? -qmax : -2*qmax/0.1*(v-1))
JuMP.@operator(model, vv, 2, vv_curve)

for gen_id in [1,2,3]
    gen = data_math["gen"]["$gen_id"]
    bus_id = gen["gen_bus"]
    # vmin = 0.9; vmax = 1.1;
    phases = gen["connections"][1:end-1]
    n = gen["connections"][end]
    @show gen_id, bus_id, phases, n
    @show qg[phases[1], gen_id], vm[phases[1],bus_id], gen["qmax"]

    @show JuMP.@constraint(model, qg[phases[1], gen_id] == vv(vm[phases[1],bus_id], gen["qmax"]) )
    # JuMP.@constraint(model, qg[phases, gen_id] == 0 )
end

###
JuMP.optimize!(model)
cost = JuMP.objective_value(model)
v = abs.(JuMP.value.(vm)) * 230.94
pg_values = JuMP.value.(pg)
qg_values = JuMP.value.(qg)
vbase = 230.94
gen_id = 1
sum(pg_values[:,gen_id])
sum(qg_values[:,gen_id])

vm_pmd = abs.(JuMP.value.(vm))[:,1]
qg_pmd = [qg_values[1,1];qg_values[2,2];qg_values[3,3]]
@show [round.(vm_pmd, digits=4); round.(qg_pmd, digits=4) ]


##
@show [round.(vm2_Control_p2g, digits=4) ; round.([qg1;qg2;qg3], digits=4) ; round.(qg1+qg2+qg3, digits=4) ]

@show [round.(vm_pmd, digits=4); round.(qg_pmd, digits=4) ; round.(sum(qg_pmd), digits=4)]

##
using Plots

plt = plot()
for gen_id in [1, 2, 3]
    gen = data_math["gen"]["$gen_id"]
    smax = 13.33
    pmax = 11.66
    gen["pmax"] = pmax
    gen["pmin"] = pmax
    gen["pg"] = pmax
    gen["qmax"] = sqrt.(smax^2 - pmax^2)
    gen["qmin"] = -gen["qmax"]
    
    plot!(0.9:0.001:1.1, vv_curve.(0.9:0.001:1.1, gen["qmax"]), label=false)
    scatter!([vm2_Control_p2g[gen_id]], [qg_dss[gen_id]], label="dss, phase $gen_id")
    scatter!([vm_pmd[gen_id]], [qg_pmd[gen_id]], label="pmd, phase $gen_id")
end
plt
