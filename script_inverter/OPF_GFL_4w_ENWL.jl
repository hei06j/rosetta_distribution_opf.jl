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
data_path = "./data/ENWL_4w_Network1_Feeder1/Master.dss"

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>1, "sb"=>"yes","warm_start_init_point"=>"yes")
# data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!, PMD.transform_loops!, PMD.reduce_lines!])
data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!, PMD.transform_loops!])
# RPMD.pv1_correction!(data_eng)
data_eng["settings"]["sbase_default"] = 1
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)

for (i, bus) in data_math["bus"]
    bus["vmin"] = [0.9 * ones(3) ; 0 ]
    bus["vmax"] = [1.1 * ones(3) ; Inf]
end

for (i, load) in data_math["load"]
    load["pd"] *= 4
    load["qd"] *= 4
end


data_math["gen"]["1"]["cost"] = [1000 0]

include("./inverter_loss_branch.jl")

for i = 1:10:length(data_math["load"])
    load = data_math["load"]["$i"]
    pd = load["pd"][1]
    gen_id = length(data_math["gen"]) + 1
    data_math["gen"]["$gen_id"] = deepcopy(data_math["gen"]["1"])
    Smax = 2 * ceil(pd) * ones(3)
    gen = data_math["gen"]["$gen_id"]
    gen["gen_bus"] = copy(load["load_bus"])
    gen["smax"] = Smax
    gen["pmax"] = 0.8 * Smax
    gen["pmin"] = 0.0 * Smax
    gen["qmax"] = sqrt.(Smax.^2 - gen["pmax"].^2)
    gen["qmin"] = -gen["qmax"]
    gen["cost"] = [10 0]
    gen["type"] = "GFL-4w"

    add_inverter_losses(data_math, gen_id)
end


ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]

##
results_GFL_4w = Dict()

objective = "cost"

model = JuMP.Model(Ipopt.Optimizer)
include("./variables.jl")
include("./constraints.jl")
include("./objectives.jl")
# include("./VV_VW_controls.jl")

###
result = JuMP.optimize!(model)
@assert(JuMP.termination_status(model) == LOCALLY_SOLVED)
obj_val_GFL = JuMP.objective_value(model)
solve_time_GFL = JuMP.solve_time(model)
iter_GFL = 20 # 42

v = value.(vr) .+ im * value.(vi)
v_axes2 = v.axes[2]
v012 = T * Array(v[1:3,:])
v2 = v012[3,:]
v2m_GFL = abs.(v2)
v0m_GFL = abs.(v012[1,:])
(v2m_max, idx) = findmax(v2m_GFL)
idx_bus = v_axes2[idx]
round.(abs.(v[:,idx_bus]), digits=4)
round.(angle.(v[:,idx_bus])*180/pi, digits=2)

c = value.(cr) .+ im * value.(ci)
c_axes2 = c.axes[2]
c012 = T * Array(c[1:3,:])
c2 = c012[3,:]
c2m_GFL = abs.(c2)

[(i, gen["pmax"][1]) for (i, gen) in data_math["gen"]]
pg_vals = Array(value.(pg)[:,1:7])
qg_vals = Array(value.(qg)[:,1:7])

gen_buses = [(parse(Int,i), gen["gen_bus"]) for (i, gen) in data_math["gen"]]
vm_gens = abs.(v[:,last.(gen_buses)])
round.(vm_gens[1:3,:] .- 1, digits=3) ./ 2


# cgens = zeros(4)
# for i = 2:7
#     gen_bus = data_math["gen"]["$i"]["gen_bus"]
#     t_bus = [branch["t_bus"] for (i,branch) in data_math["branch"] if branch["f_bus"]==gen_bus]
#     arc = [(l,i,j) for (l,i,j) in ref[:arcs_branch] if i == gen_bus]
#     cgens = [cgens abs.(c[:,arc])]
# end

# cgens_ratings= [3*data_math["gen"]["$i"]["smax"][1] * 1000 / (253*3) / 4.33 for i in 2:7]


##
# # using Plots
# histogram(c2m.*100, label=false, xlabel="Current Negative Sequence (%)")

# histogram(v2m.*100, label=false, xlabel="Voltage Negative Sequence (%)")
# vline!([2], label="Grid Code Limit", linewidth=2, legend=:topright)


# plot(abs.(v))
