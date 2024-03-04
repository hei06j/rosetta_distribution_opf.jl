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

ODSS.dss("Redirect ./data/test_gen_3w_wye_InvControl.dss")
vm_Control = ODSS.Circuit.AllBusVolts()
cm_Control = ODSS.PDElements.AllCurrentsAllCurrents()[5:8]
vm1_Control = abs.(vm_Control)[1:3]
vm2_Control = abs.(vm_Control)[4:7]
pg_Control = ODSS.PVsystems.kW()
qg_Control = ODSS.PVsystems.kvar()

# plot([vm1_noControl[1:3] vm2_noControl[1:3]]', marker=:circle, linestyle=:dash, color=[:blue :red :green], label=["v_a" "v_b" "v_c"])
# plot!([vm1_Control[1:3] vm2_Control[1:3]]', marker=:circle, linestyle=:solid, color=[:blue :red :green], label=["v_a^c" "v_b^c" "v_c^c"])

[pg_noControl pg_Control]
[qg_noControl qg_Control]
qg_Control
round.(vm2_Control./230.94, digits=4)

##
data_path = "./data/inverter_3w_wye_unbalanced_loads.dss"

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!, PMD.transform_loops!])
RPMD.pv1_correction!(data_eng)
data_eng["settings"]["sbase_default"] = 1
vbase = data_eng["settings"]["vbases_default"]["b1"]
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)

data_math["bus"]["1"]["vmin"] = [0.9 * ones(3) ; 0 ]
data_math["bus"]["1"]["vmax"] = [1.1 * ones(3) ; Inf]

gen_id = 1
gen = data_math["gen"]["$gen_id"]
smax = 40
pmax = 35
gen["pmax"] = pmax/3 * ones(3)
gen["pmin"] = pmax/3 * ones(3)
gen["pg"] = pmax/3 * ones(3)
gen["qmax"] = sqrt.(smax^2 - pmax^2)/3 * ones(3)
gen["qmin"] = -gen["qmax"]

gen["connections"] = gen["connections"][1:3]
data_math["bus"]["1"]["terminals"] = data_math["bus"]["1"]["terminals"][1:4]
data_math["bus"]["1"]["grounded"] = data_math["bus"]["1"]["grounded"][1:4]

gen["cost"] = [0 0]
data_math["gen"]["2"]["cost"] = [1000 0]

ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]

model = JuMP.Model(Ipopt.Optimizer)
include("./core/variables.jl")

for (id, generator) in ref[:gen]
    ### if neutral is in generator connections, then false, otherwise it is true
    explicit_neutral = 4 in generator["connections"]
    nphases = RPMD._infer_int_dim_unit(generator, !explicit_neutral)
    # nphases = RPMD._infer_int_dim_unit(generator, false)
    bus_id = generator["gen_bus"]
    bus = ref[:bus][bus_id]
    configuration = generator["configuration"]
    connections = generator["connections"]

    N = length(connections)
    pmin = get(generator, "pmin", fill(-Inf, N))
    pmax = get(generator, "pmax", fill( Inf, N))
    qmin = get(generator, "qmin", fill(-Inf, N))
    qmax = get(generator, "qmax", fill( Inf, N))

    # constraint_mc_generator_current(pm, id)
    if configuration==PMD.WYE || length(pmin)==1 || nphases==1
        if explicit_neutral
            phases = connections[1:end-1]
            n = connections[end]
        else
            phases = connections
            n = 4
        end
            crg_bus[id] = JuMP.Containers.DenseAxisArray([crg[phases,id]..., -sum(crg[phases,id])], connections)
            cig_bus[id] = JuMP.Containers.DenseAxisArray([cig[phases,id]..., -sum(cig[phases,id])], connections)
            # JuMP.@constraint(model,  pg[phases,id] .==  (vr[phases,bus_id] .- vr[n,bus_id]) .* crg[phases,id] .+ (vi[phases,bus_id] .- vi[n,bus_id]) .* cig[phases,id])
            JuMP.@constraint(model,  pg[phases,id] .==  vr[phases,bus_id] .* crg[phases,id] .+ vi[phases,bus_id] .* cig[phases,id])
            
            nonInf_pmin_pmax = [idx for (idx,t) in enumerate(pmin) if pmin[idx].>-Inf || pmax[idx].<Inf]
            JuMP.@constraint(model, pmin[nonInf_pmin_pmax] .<= pg[nonInf_pmin_pmax,id] .<= pmax[nonInf_pmin_pmax] )
    
            nonInf_qmin_qmax = [idx for (idx,t) in enumerate(qmin) if qmin[idx].>-Inf || qmax[idx].<Inf]
            # JuMP.@constraint(model, qg[nonInf_qmin_qmax,id] .== -(vr[nonInf_qmin_qmax,bus_id] .- vr[n,bus_id]) .* cig[nonInf_qmin_qmax,id] .+ (vi[nonInf_qmin_qmax,bus_id] .- vi[n,bus_id]) .* crg[nonInf_qmin_qmax,id])
            JuMP.@constraint(model, qg[nonInf_qmin_qmax,id] .== -vr[nonInf_qmin_qmax,bus_id] .* cig[nonInf_qmin_qmax,id] .+ vi[nonInf_qmin_qmax,bus_id] .* crg[nonInf_qmin_qmax,id])
            JuMP.@constraint(model, qmin[nonInf_qmin_qmax] .<= qg[nonInf_qmin_qmax,id] .<= qmax[nonInf_qmin_qmax] )
    
        # else
        #     phases = connections
        #     crg_bus[id] = crg[:,id]
        #     cig_bus[id] = cig[:,id]
        #     JuMP.@constraint(model,  pg[phases,id] .==  vr[phases,bus_id] .* crg[phases,id] .+ vi[phases,bus_id] .* cig[phases,id])

        #     nonInf_pmin_pmax = [idx for (idx,t) in enumerate(pmin) if pmin[idx].>-Inf || pmax[idx].<Inf]
        #     JuMP.@constraint(model, pmin[nonInf_pmin_pmax] .<= pg[nonInf_pmin_pmax,id] .<= pmax[nonInf_pmin_pmax] )
        #     # JuMP.@constraint(model, pmin .<= pg[phases,id] .<= pmax )

        #     nonInf_qmin_qmax = [idx for (idx,t) in enumerate(qmin) if qmin[idx].>-Inf || qmax[idx].<Inf]
        #     JuMP.@constraint(model, qg[nonInf_qmin_qmax,id] .== -vr[nonInf_qmin_qmax,bus_id] .* cig[nonInf_qmin_qmax,id] .+ vi[nonInf_qmin_qmax,bus_id] .* crg[nonInf_qmin_qmax,id])
        #     JuMP.@constraint(model, qmin[nonInf_qmin_qmax] .<= qg[nonInf_qmin_qmax,id] .<= qmax[nonInf_qmin_qmax] )
        #     # JuMP.@constraint(model, qmin .<= qg[phases,id] .<= qmax )
        # end
        
    else ## configuration==PMD.DELTA

        Md = PMD._get_delta_transformation_matrix(length(connections))
        crg_bus[id] = JuMP.Containers.DenseAxisArray(Md'*Vector{JuMP.AffExpr}(crg[connections,id]), connections)
        cig_bus[id] = JuMP.Containers.DenseAxisArray(Md'*Vector{JuMP.AffExpr}(cig[connections,id]), connections)

        connections2 = [connections[2:end]..., connections[1]]
        vrg = Vector{JuMP.AffExpr}(vr[connections,bus_id]) .- Vector{JuMP.AffExpr}(vr[connections2,bus_id])
        vig = Vector{JuMP.AffExpr}(vi[connections,bus_id]) .- Vector{JuMP.AffExpr}(vi[connections2,bus_id])

        JuMP.@constraint(model, pg[connections,id] .==  vrg .* crg[connections,id] .+ vig .* cig[connections,id])
        JuMP.@constraint(model, qg[connections,id] .== -vrg .* cig[connections,id] .+ vig .* crg[connections,id])

        JuMP.@constraint(model, pmin[connections] .<= pg[connections,id] .<= pmax[connections])
        JuMP.@constraint(model, qmin[connections] .<= qg[connections,id] .<= qmax[connections])
    end
end
crg_bus = JuMP.Containers.DenseAxisArray([t in crg_bus[i].axes[1] ? crg_bus[i][t] : 0 for t in 1:n_ph, i in keys(ref[:gen])], 1:n_ph, keys(ref[:gen]))
cig_bus = JuMP.Containers.DenseAxisArray([t in cig_bus[i].axes[1] ? cig_bus[i][t] : 0 for t in 1:n_ph, i in keys(ref[:gen])], 1:n_ph, keys(ref[:gen]))

objective = "cost"

include("./core/constraints.jl")
include("./core/objectives.jl")


### Inverter control Volt-var
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

vv_curve(v, qmax) = v <= 0.95 ? qmax : (v >= 1.05 ? -qmax : -2*qmax/0.1*(v-1))
JuMP.@operator(model, vv, 2, vv_curve)

# JuMP.@constraint(model, [i in phases], qg[i,gen_id] == vv(vm_pv[i,bus_id], gen["qmax"][i]) )
JuMP.@constraint(model, qg[1, gen_id] == vv(sum(vm[i,bus_id] for i in ref[:bus][1]["terminals"][1:3])/3, gen["qmax"][1]) )
JuMP.@constraint(model, qg[2, gen_id] == qg[1, gen_id])
JuMP.@constraint(model, qg[3, gen_id] == qg[1, gen_id])

###
JuMP.optimize!(model)
cost = JuMP.objective_value(model)
v = abs.(JuMP.value.(vm)) * 230.94
pg_values = JuMP.value.(pg)
qg_values = JuMP.value.(qg)
vbase
sum(pg_values[:,gen_id])
sum(qg_values[:,gen_id])
round.(abs.(JuMP.value.(vm))[:,1], digits=4)

##
@show [round.(abs.(JuMP.value.(vm))[:,1], digits=4) ; sum(qg_values[:,gen_id])]

@show [round.(vm2_Control./230.94, digits=4) ; qg_Control]

