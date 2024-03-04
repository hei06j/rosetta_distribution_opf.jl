using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP  # bl/array_nl
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels
import LinearAlgebra: diag, diagm


###
data_path = "./data/inverter_GFM_4w_wye_unbalanced_loads.dss"

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!])
# RPMD.pv1_correction!(data_eng)
data_eng["settings"]["sbase_default"] = 1
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)
data_math["bus"]["1"]["vmin"] = [0.9 * ones(3) ; 0 ]
data_math["bus"]["1"]["vmax"] = [1.1 * ones(3) ; Inf]

### - reomve the voltage source gen
### - remove the voltage source bus 
### - remove the voltage source branch that connects to the load bus

### - add the inverter internal bus: connect the inverter to this bus
###   - terminals = [1,2,3,4]
###   - grounded = [0,0,0,1]
### - add the inverter internal branch: f_bus and t_bus are new and old inverter buses

voltage_source_id = "2"
voltage_source_bus = data_math["gen"][voltage_source_id]["gen_bus"]
voltage_source_branch = [i for (i,branch) in data_math["branch"] if branch["f_bus"]==voltage_source_bus || branch["t_bus"]==voltage_source_bus][1]
branch_copy = deepcopy(data_math["branch"][voltage_source_branch])
# delete!(data_math["gen"], "$voltage_source_id")
# delete!(data_math["bus"], "$voltage_source_bus")
# delete!(data_math["branch"], voltage_source_branch)


gen_id = 1
gen = data_math["gen"]["$gen_id"]
smax = 40
pmax = 10
gen["pmax"] = pmax/3 * ones(3)
gen["pmin"] = zeros(3)
gen["qmax"] = sqrt.(smax^2 - pmax^2)/3 * ones(3)
gen["qmin"] = -gen["qmax"]
gen["pg"] = gen["pmax"]./2
gen["qg"] = zeros(3)

data_math["gen"]["1"]["cost"] = [10 0]
data_math["gen"]["2"]["cost"] = [10 0]

gen["Dp"] = 2*diag(data_math["branch"]["1"]["br_r"])[1:3]
gen["Dq"] = 2*diag(data_math["branch"]["1"]["br_x"])[1:3]

# gen["Dp"] = 0.01 * ones(3)
# gen["Dq"] = 0.02 * ones(3)

include("./core/inverter_loss_branch.jl")
add_inverter_losses(data_math, gen_id, GFM=true)

ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]

model =[]
GFM_4w_controls = Dict()
control_forming = []
objective = "cost"

##

for control_forming in ["setpoint", "no_setpoint_droop", "droop"]

    model = JuMP.Model(ipopt_solver)
    include("./core/variables.jl")
    include("./core/constraints.jl")
    include("./core/objectives.jl")


    JuMP.optimize!(model)
    @assert(JuMP.termination_status(model) == LOCALLY_SOLVED)
    obj_value = JuMP.objective_value(model)

    ###

    pg_vals = JuMP.value.(pg)
    qg_vals = JuMP.value.(qg)
    vr_vals = JuMP.value.(vr)
    vi_vals = JuMP.value.(vi)
    v = vr_vals .+ im*vi_vals
    vmag = abs.(v)
    va = angle.(v) .* 180/pi 
    # v1_012 = T * v[1:3,1]
    # v1_012_m = abs.(v1_012)
    crg_values = JuMP.value.(crg_bus)
    cig_values = JuMP.value.(cig_bus)
    cg = Array(crg_values[:,1] + im * cig_values[:,1])

    cost = sum(ref[:gen][1]["cost"][1]*sum(pg_vals[:,1]) +  ref[:gen][1]["cost"][2])

    branch_idx = (1, 2, 1)
    cr_vals = value.(cr)
    ci_vals = value.(ci)
    c = cr_vals[:,branch_idx] .+ im *  ci_vals[:,branch_idx]
    c012 = T * Array(c[1:3])
    c2 = abs.(c012)[3]

    vmag = abs.(v)
    vangle = angle.(v).*180/pi

    GFM_4w_controls[control_forming] = Dict()
    GFM_4w_controls[control_forming]["pg"] = Array(pg_vals[:,1])
    GFM_4w_controls[control_forming]["qg"] = Array(qg_vals[:,1])
    GFM_4w_controls[control_forming]["v"] = v
    GFM_4w_controls[control_forming]["cr"] = c
    GFM_4w_controls[control_forming]["cr_2"] = c2
    GFM_4w_controls[control_forming]["cg"] = cg
    GFM_4w_controls[control_forming]["cost"] = cost
    GFM_4w_controls[control_forming]["objective"] = obj_value

end

##

vm_inverter = round.([abs.(GFM_4w_controls["setpoint"]["v"][:,3]) abs.(GFM_4w_controls["no_setpoint_droop"]["v"][:,3]) abs.(GFM_4w_controls["droop"]["v"][:,3])]', digits=3)
va_inverter = round.([angle.(GFM_4w_controls["setpoint"]["v"][:,3]) angle.(GFM_4w_controls["no_setpoint_droop"]["v"][:,3]) angle.(GFM_4w_controls["droop"]["v"][:,3])]'.*180/pi, digits=1)

pg_inverter = round.([GFM_4w_controls["setpoint"]["pg"] GFM_4w_controls["no_setpoint_droop"]["pg"] GFM_4w_controls["droop"]["pg"]], digits=2)
qg_inverter = round.([GFM_4w_controls["setpoint"]["qg"] GFM_4w_controls["no_setpoint_droop"]["qg"] GFM_4w_controls["droop"]["qg"]], digits=2)
cg_inverter = round.([GFM_4w_controls["setpoint"]["cg"] GFM_4w_controls["no_setpoint_droop"]["cg"] GFM_4w_controls["droop"]["cg"]]', digits=2)

costs = [GFM_4w_controls["setpoint"]["cost"] GFM_4w_controls["no_setpoint_droop"]["cost"] GFM_4w_controls["droop"]["cost"]]
costs = [GFM_4w_controls["setpoint"]["objective"] GFM_4w_controls["no_setpoint_droop"]["objective"] GFM_4w_controls["droop"]["objective"]]
costs = round.([GFM_4w_controls["setpoint"]["cr_2"] GFM_4w_controls["no_setpoint_droop"]["cr_2"] GFM_4w_controls["droop"]["cr_2"]], digits=2)

cg_inverter = abs.(round.([GFM_4w_controls["setpoint"]["cg"] GFM_4w_controls["no_setpoint_droop"]["cg"] GFM_4w_controls["droop"]["cg"]]', digits=2))

##
pg_gen = Array(value.(pg))[:,1]
qg_gen = Array(value.(qg))[:,1]
v_gen = Array(value.(vr) .+ im*value.(vi))[:,3]
abs.(v_gen)
angle.(v_gen).*180/pi
pmax = data_math["gen"]["1"]["pmax"]
qmax = data_math["gen"]["1"]["pmax"]

Dp = (abs.(v_gen[1:3]) .- ones(3)) ./ (pmax .- pg_gen[1:3])
Dq = (abs.(v_gen[1:3]) .- ones(3)) ./ (qmax .- qg_gen[1:3])

##
function plot_phasors(phasor, Imax, color; plus=false, equal=false, labeled=false)
    plt = Plots.plot([0,real.(phasor[1])],[0,imag.(phasor[1])], arrow=true, color=color, linewidth=2, linestyle=:solid, label="a", border=:none)
    Plots.plot!([0,real.(phasor[2])],[0,imag.(phasor[2])], arrow=true, color=color, linewidth=2, linestyle=:solid, label="b", border=:none)
    Plots.plot!([0,real.(phasor[3])],[0,imag.(phasor[3])], arrow=true, color=color, linewidth=2, linestyle=:solid, label="c", border=:none)
    if phasor[4] !==  0 + 0im
        Plots.plot!([0,real.(phasor[4])],[0,imag.(phasor[4])], arrow=true, color=color, linewidth=2, linestyle=:solid, label="n", border=:none)
    end
    Plots.plot!([0,0],[0,1.1*Imax], arrow=true, color=:grey, linestyle=:dot, label=false)
    Plots.plot!([0,1.1*Imax*real(exp(im*210/180*pi))],[0,1.1*Imax*imag(exp(im*210/180*pi))], arrow=true, color=:grey, linestyle=:dot, label=false)
    Plots.plot!([0,1.1*Imax*real(exp(im*330/180*pi))],[0,1.1*Imax*imag(exp(im*330/180*pi))], arrow=true, color=:grey, linestyle=:dot, label=false)
    if labeled
        Plots.plot!(Imax*exp.(im*(0:0.01:2pi)), color=:black, border=:none, label=false, markersize=10, size=(300,300), legend=:topleft)
    else
        Plots.plot!(Imax*exp.(im*(0:0.01:2pi)), color=:black, border=:none, label=false, markersize=10, size=(300,300), legend=false)
    end
    return plt
end


Vmax = maximum([abs.(c_set_nodroop)  abs.(c_notset_nodroop) abs.(c_notset_droop)])
plot_phasors(c_set_nodroop, Vmax, :blue; plus=false, equal=false, labeled=false)
plot_phasors(c_notset_nodroop, Vmax, :red; plus=false, equal=false, labeled=false)
plot_phasors(c_notset_droop, Vmax, :red; plus=false, equal=false, labeled=false)



##
crg_values = JuMP.value.(crg_bus)
cig_values = JuMP.value.(cig_bus)

val_cg = crg_values[1:3,1] .+ im .* cig_values[1:3,1]
val_cg_012 = T * Array(val_cg)
val_cmg_012 = abs.(val_cg_012)

abs.(val_cg)

abs.(crg_values[:,1] .+ im .* cig_values[:,1])
