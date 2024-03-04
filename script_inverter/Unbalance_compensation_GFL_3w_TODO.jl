using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP  # bl/array_nl
import LinearAlgebra: diag, diagm
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels


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

for (i, bus) in data_math["bus"]
    bus["vmin"] = [0.9 * ones(3) ; 0 ]
    bus["vmax"] = [1.1 * ones(3) ; Inf]
end


gen_id = 1
gen = data_math["gen"]["$gen_id"]
gen["pmax"] = 23/3 * ones(3)
gen["pmin"] = zeros(3)
# gen["pmax"] = 23/3 * ones(3)
# gen["pmin"] = 23/3 * ones(3)
# gen["pg"] = 23/3 * ones(3)
gen["qmax"] = sqrt.(25^2 - 23^2)/3 * ones(3)
gen["qmin"] = -gen["qmax"]

gen["connections"] = gen["connections"][1:3]
data_math["bus"]["1"]["terminals"] = data_math["bus"]["1"]["terminals"][1:4]
data_math["bus"]["1"]["grounded"] = data_math["bus"]["1"]["grounded"][1:4]

old_gen_bus = gen["gen_bus"]
new_gen_bus = length(data_math["bus"]) + 1
# data_math["bus"]["$gen_bus"]["bus_type"] = 1
data_math["bus"]["$new_gen_bus"] = deepcopy(data_math["bus"]["$old_gen_bus"])
data_math["bus"]["$new_gen_bus"]["bus_i"] = new_gen_bus
data_math["bus"]["$new_gen_bus"]["index"] = new_gen_bus
data_math["bus"]["$new_gen_bus"]["name"] = "inverter_bus"
data_math["bus"]["$new_gen_bus"]["type"] = "GFL"
data_math["bus"]["$new_gen_bus"]["terminals"] = data_math["bus"]["$new_gen_bus"]["terminals"][1:3]  # 3-wire
data_math["bus"]["$new_gen_bus"]["grounded"] = data_math["bus"]["$new_gen_bus"]["grounded"][1:3]    # 3-wire
data_math["bus"]["$new_gen_bus"]["vmin"] = data_math["bus"]["$new_gen_bus"]["vmin"][1:3]            # 3-wire
data_math["bus"]["$new_gen_bus"]["vmax"] = data_math["bus"]["$new_gen_bus"]["vmax"][1:3]            # 3-wire
data_math["bus"]["$new_gen_bus"]["bus_type"] = 1                    # GFL
# data_math["bus"]["$new_gen_bus"]["vm"] = [1. ; 1 ; 1 ; 0]         # GFL
# data_math["bus"]["$new_gen_bus"]["va] = [0 ; -2. ; 1 ; 0]         # GFL
# data_math["bus"]["$new_gen_bus"]["grounded"] = Bool[0, 0, 0, 1]   # GFL
# data_math["bus"]["$old_gen_bus"]["bus_type"] = 1                  # GFL
gen["gen_bus"] = new_gen_bus

Rf = 0.015
Lf = 0.42E-3
Cf = 0.33E-9
new_branch_id = length(data_math["branch"]) + 1
data_math["branch"]["$new_branch_id"] = deepcopy(data_math["branch"]["$(new_branch_id-1)"])
data_math["branch"]["$new_branch_id"]["index"] = new_branch_id
data_math["branch"]["$new_branch_id"]["name"] = "GFL_internal_z"
zbase = 230.94^2 / 1000
data_math["branch"]["$new_branch_id"]["br_r"] = diagm(Rf/zbase * ones(3))
data_math["branch"]["$new_branch_id"]["br_x"] = diagm(Lf*2*pi*50/zbase * ones(3))
data_math["branch"]["$new_branch_id"]["b_to"] = diagm(Cf*2*pi*50*zbase * ones(3))
data_math["branch"]["$new_branch_id"]["b_fr"] = zeros(3,3)
data_math["branch"]["$new_branch_id"]["g_fr"] = zeros(3,3)
data_math["branch"]["$new_branch_id"]["g_to"] = zeros(3,3)
data_math["branch"]["$new_branch_id"]["f_bus"] = new_gen_bus
data_math["branch"]["$new_branch_id"]["t_bus"] = old_gen_bus

data_math["branch"]["$new_branch_id"]["f_connections"] = data_math["branch"]["$new_branch_id"]["f_connections"][1:3]
data_math["branch"]["$new_branch_id"]["t_connections"] = data_math["branch"]["$new_branch_id"]["t_connections"][1:3]

data_math["gen"]["1"]["cost"] = [0 0]
data_math["gen"]["2"]["cost"] = [1000 0]

ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]


##
results_GFL_3w = Dict()

for objective in ["cost", "VUF", "VUF2", "PVUR", "LVUR", "IUF", "IUF2", "PIUR"]

    model = JuMP.Model(Ipopt.Optimizer)
    include("./core/variables.jl")
    include("./core/constraints.jl")
    include("./core/objectives.jl")

    # objective = "VUF"

    # ### objectives: cost, VUF, VUF2, PVUR, LVUR, IUF, IUF2, PIUR, PPUR (x), PQUR (x)
    # if objective == "cost"
    #     JuMP.@objective(model, Min, sum(gen["cost"][1]*sum(pg[:,i]) + gen["cost"][2] for (i,gen) in ref[:gen]))

    # elseif objective in ["VUF" "VUF2" "PVUR" "LVUR"]
    #     # terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref[:bus])
    #     # terminals = Dict(i => collect(1:3) for (i, bus) in ref[:bus])
    #     # vm = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm_$i", lower_bound = 0 ) for i in keys(ref[:bus]))
    #     # vm = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm[i].axes[1] ? vm[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
    #     # vr_012 = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vr_012") for i in keys(ref[:bus]))
    #     # vr_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vr_012[i].axes[1] ? vr_012[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
    #     # vi_012 = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vi_012") for i in keys(ref[:bus]))
    #     # vi_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vi_012[i].axes[1] ? vi_012[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
    #     # vm_012 = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm_012") for i in keys(ref[:bus]))
    #     # vm_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm_012[i].axes[1] ? vm_012[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))

    #     vm = JuMP.@variable(model, [t in 1:3], base_name="vm", lower_bound = 0 )
    #     vr_012 = JuMP.@variable(model, [t in 1:3], base_name="vr_012")
    #     vi_012 = JuMP.@variable(model, [t in 1:3], base_name="vi_012")
    #     vm_012 = JuMP.@variable(model, [t in 1:3], base_name="vm_012", lower_bound = 0)

    #     # for (i, bus) in ref[:bus]
    #         i = 1
    #         bus = ref[:bus][i]
    #         # JuMP.@constraint(model, vm[:,i].^2 .== vr[:,i].^2 .+ vi[:,i].^2)
    #         terminals = ref[:bus][i]["terminals"][1:end-1]
    #         JuMP.@constraint(model, vm[terminals,i].^2 .== vr[terminals,i].^2 .+ vi[terminals,i].^2)
    #         JuMP.@constraint(model, vr_012[terminals,i] .== Tre * Array(vr[terminals,i]) .- Tim * Array(vi[terminals,i]))
    #         JuMP.@constraint(model, vi_012[terminals,i] .== Tre * Array(vi[terminals,i]) .+ Tim * Array(vr[terminals,i]))
    #         JuMP.@constraint(model, vm_012[terminals,i].^2 .== vr_012[terminals,i].^2 .+ vi_012[terminals,i].^2)
    #     # end

    #     if objective == "VUF"
    #         JuMP.@objective(model, Min, vm_012[3,1] / vm_012[2,1])

    #     elseif objective == "VUF2"
    #         JuMP.@objective(model, Min, vm_012[3,1])

    #     elseif objective == "PVUR"
    #         phase_voltage = JuMP.@variable(model, [i=1], base_name="phase_voltage_$i")
    #         # phase_voltage = JuMP.@variable(model, [i in keys(ref[:bus])], base_name="phase_voltage_$i")
    #         # for (i, bus) in ref[:bus]
    #             i = 1
    #             bus = ref[:bus][i]
    #             terminals = ref[:bus][i]["terminals"][1:end-1]
    #             JuMP.@constraint(model, [t in terminals], phase_voltage[i] >= vm[t,i] - sum(vm[terminals,i])/3)
    #         # end
    #         JuMP.@objective(model, Min, phase_voltage[1] / (sum(vm[terminals,1])/3) )

    #     elseif objective == "LVUR"
    #         # line_voltage = JuMP.@variable(model, [i in keys(ref[:bus])], base_name="line_voltage_$i")
    #         # terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref[:bus])
    #         # vm_ll = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm_ll_$i") for i in keys(ref[:bus]))
    #         # vm_ll = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm_ll[i].axes[1] ? vm_ll[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
    #         line_voltage = JuMP.@variable(model, base_name="line_voltage")
    #         vm_ll = JuMP.@variable(model, [t in 1:3], base_name="vm_ll", lower_bound=0)
    #         vr_ll = JuMP.@variable(model, [t in 1:3], base_name="vr_ll")
    #         vi_ll = JuMP.@variable(model, [t in 1:3], base_name="vi_ll")
    #         # for (i, bus) in ref[:bus]
    #             # terminals = ref[:bus][i]["terminals"][1:end-1]
    #             # JuMP.@constraint(model, vm_ll[terminals,i] .== Array(vm[terminals,i]) .- Array(vm[terminals2,i]))
    #             # JuMP.@constraint(model, [t in terminals], line_voltage[i] >= vm_ll[t,i] - sum(vm_ll[terminals,i])/3)
    #             i = 1
    #             bus = ref[:bus][i]
    #             terminals = collect(1:3)
    #             terminals2 = [terminals[2:end]..., terminals[1]]
    #             # JuMP.@constraint(model, vm_ll[terminals] .== vm[terminals] .- vm[terminals2])
    #             JuMP.@constraint(model, vr_ll[terminals] .== Array(vr[terminals,i]) .- Array(vr[terminals2,i]))
    #             JuMP.@constraint(model, vi_ll[terminals] .== Array(vi[terminals,i]) .- Array(vi[terminals2,i]))
    #             JuMP.@constraint(model, vm_ll[terminals].^2 .== vr_ll[terminals].^2 .+ vi_ll[terminals].^2)

    #             JuMP.@constraint(model, line_voltage .>= vm_ll[terminals] .- sum(vm_ll[terminals])/3)
    #         # end
    #         JuMP.@objective(model, Min, line_voltage / (sum(vm_ll[terminals,1])/3) )
    #     end


    # elseif objective in ["IUF" "IUF2" "PIUR"]
    #     int_dim = Dict(i => RPMD._infer_int_dim_unit(gen, !(4 in gen["connections"])) for (i,gen) in ref[:gen])
    #     cmg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cmg_$i", lower_bound=0) for i in keys(ref[:gen]))
    #     cmg = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? cmg[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
    #     crg_012 = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="crg_012_$i") for i in keys(ref[:gen]))
    #     crg_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? crg_012[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
    #     cig_012 = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cig_012_$i") for i in keys(ref[:gen]))
    #     cig_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? cig_012[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
    #     cmg_012 = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cmg_012_$i", lower_bound=0) for i in keys(ref[:gen]))
    #     cmg_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? cmg_012[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))

    #     for (i, gen) in ref[:gen]
    #         if 4 in gen["connections"]
    #             phases = ref[:gen][i]["connections"][1:end-1]
    #         else
    #             phases = ref[:gen][i]["connections"]
    #         end
    #         JuMP.@constraint(model, cmg[phases,i].^2 .== crg_bus[phases,i].^2 .+ crg_bus[phases,i].^2)
    #         JuMP.@constraint(model, crg_012[phases,i] .== Tre * Array(crg_bus[phases,i]) .- Tim * Array(cig_bus[phases,i]))
    #         JuMP.@constraint(model, cig_012[phases,i] .== Tre * Array(cig_bus[phases,i]) .+ Tim * Array(crg_bus[phases,i]))
    #         JuMP.@constraint(model, cmg_012[phases,i].^2 .== crg_012[phases,i].^2 .+ cig_012[phases,i].^2)
    #     end

    #     if objective == "IUF"
    #         JuMP.@objective(model, Min, cmg_012[3,1] / cmg_012[2,1])

    #     elseif objective == "IUF2"
    #         JuMP.@objective(model, Min, cmg_012[3,1])
            
    #     elseif objective == "PIUR"
    #         phase_current = JuMP.@variable(model, [i in keys(ref[:gen])], base_name="phase_current_$i")
    #         # for (i, gen) in ref[:gen]
    #             i = 1
    #             gen = ref[:gen][i]
    #             connections = ref[:gen][i]["connections"][1:end-1]
    #             JuMP.@constraint(model, [t in connections], phase_current[i] >= cmg[t,i] - sum(cmg[connections,i])/3)
    #         # end
    #         JuMP.@objective(model, Min, phase_current[1] / (sum(cmg[connections,1])/3) )
    #     end


    # elseif objective == "PPUR"
    #     phase_p = JuMP.@variable(model, base_name="phase_p", lower_bound=0)
    #     # phase_p = JuMP.@variable(model, [i in keys(ref[:gen])], base_name="phase_p_$i")
    #     # for (i, gen) in ref[:gen]
    #         i = 1
    #         gen = ref[:gen][i]
    #         connections = gen["connections"][1:end-1]
    #         JuMP.@constraint(model, phase_p .>= pg[connections,i] .- sum(pg[connections,i])/3)
    #         # JuMP.@constraint(model, [t in connections], phase_p[i] >= pg[t,i] - sum(pg[connections,i])/3)
    #     # end
    #     # JuMP.@objective(model, Min, phase_p[1] / (sum(pg[connections,1])/3) )
    #     JuMP.@objective(model, Min, phase_p / (sum(pg[connections,1])/3) )


    # elseif objective == "PQUR"
    #     phase_q = JuMP.@variable(model, [i=1], base_name="phase_q_$i")
    #     # phase_q = JuMP.@variable(model, [i in keys(ref[:gen])], base_name="phase_q_$i")
    #     # for (i, gen) in ref[:gen]
    #         i = 1
    #         gen = ref[:gen][i]
    #         connections = ref[:gen][i]["connections"][1:end-1]
    #         JuMP.@constraint(model, [c in connections], phase_q[i] >= qg[connections,i] - sum(qg[connections,i])/3)
    #     # end
    #     JuMP.@objective(model, Min, phase_q[1] / (sum(qg[connections,1])/3) )

    # end

    # include("./core/VV_VW_controls.jl")

    JuMP.optimize!(model)
    cost = JuMP.objective_value(model)

    ##
    gen_id = [parse(Int,i) for (i,gen) in data_math["gen"] if occursin("pv", gen["name"])][1]
    gen_bus_id = data_math["gen"]["$gen_id"]["gen_bus"]

    vr_vals = value.(vr)
    vi_vals = value.(vi)
    v1 = vr_vals[:,gen_bus_id] + im * vi_vals[:,gen_bus_id]
    vm1 = abs.(v1)
    # va1 = angle.(v1).*180/pi
    # va1_pn = va1[1:3] .- va1[4]
    # va1_pn_pp = Array(va1_pn[[1,2,3]]) .- Array(va1_pn[[2,3,1]])
    # va1_pp = Array(va1[[1,2,3]]) .- Array(va1[[2,3,1]])
    v1_012 = T * v1[1:3]
    vm1_012 = abs.(v1_012)

    crg_values = JuMP.value.(crg_bus)
    cig_values = JuMP.value.(cig_bus)
    cg1 = crg_values[:,gen_id] + im * cig_values[:,gen_id]
    cgm1 = abs.(cg1)
    cg1_012 = T * cg1[1:3]
    cgm1_012 = abs.(cg1_012)
    # cga1 = angle.(cg1).*180/pi

    crd_values = JuMP.value.(crd_bus)
    cid_values = JuMP.value.(cid_bus)
    cd1 = sum(Array(crd_values .+ im * cid_values), dims=2)
    # cdm1 = abs.(cd1)
    cd1_012 = T * cd1[1:3]
    cdm1_012 = abs.(cd1_012)

    branch_lij = (1, 2, 1)
    cr1 = value.(cr_bus)[:,branch_lij]
    ci1 = value.(ci_bus)[:,branch_lij]
    c1 = cr1 .+ im * ci1
    # cm1 = abs.(c1)
    c1_012 = T * c1[1:3]
    cm1_012 = abs.(c1_012)

    results_GFL_3w[objective] = Dict()
    results_GFL_3w[objective]["vm_gen_012"] = vm1_012
    results_GFL_3w[objective]["cm_gen_012"] = cgm1_012
    results_GFL_3w[objective]["cm_load_012"] = cdm1_012
    results_GFL_3w[objective]["cm_branch_012"] = cm1_012
end



##
using Plots
objectives = ["cost", "VUF", "VUF2", "PVUR", "LVUR", "IUF", "IUF2", "PIUR"]
vm_zero_seq = [results_GFL_3w[obj]["vm_gen_012"][1] for obj in objectives]
vm_positive_seq = [results_GFL_3w[obj]["vm_gen_012"][2] for obj in objectives]
vm_negative_seq = [results_GFL_3w[obj]["vm_gen_012"][3] for obj in objectives]
cgm_zero_seq = [results_GFL_3w[obj]["cm_gen_012"][1] for obj in objectives]
cgm_positive_seq = [results_GFL_3w[obj]["cm_gen_012"][2] for obj in objectives]
cgm_negative_seq = [results_GFL_3w[obj]["cm_gen_012"][3] for obj in objectives]


scatter(objectives, vm_zero_seq, label="Vm0", xticks=(0.5:7.5, objectives))
scatter!(objectives, vm_negative_seq, label="Vm2")

scatter(objectives, cgm_zero_seq, label="Ig0", xticks=(0.5:7.5, objectives))
scatter!(objectives, cgm_negative_seq, label="Ig2")

##
Array(c1) .- Array(cd1) .+ Array(cg1)
c1_012 .- cd1_012 .+ cg1_012
[c1_012  cd1_012  cg1_012]
cm1_012 .- cdm1_012 .+ cgm1_012

a4w = [1.30852-2.32647im   1.90377-2.1503im    0.595248+0.176166im
0.989932-2.41685im   5.04945-2.59949im    4.05952-0.182641im
0.134892+0.824484im  3.09972-0.260136im   2.96483-1.08462im]

a3w = [ 1.63193-2.74083im   1.91129-2.20162im   0.279359+0.539211im
0.682041-2.16986im   5.01165-2.637im      4.32961-0.467146im
0.021938+0.823998im   3.0776-0.358277im   3.05566-1.18227im]