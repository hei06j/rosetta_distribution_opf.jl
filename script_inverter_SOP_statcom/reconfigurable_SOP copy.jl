using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
# import InfrastructureModels
# using AppleAccelerate
# using HSL_jll
using Ipopt
# using Juniper
# using HiGHS
using JuMP
# using LinearAlgebra
# import LinearAlgebra: diag, diagm
using LaTeXStrings

const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
# const IM = InfrastructureModels

# ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes", "hsllib"=>HSL_jll.libhsl_path, "linear_solver"=>"ma86")
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")

setting = Dict("conventional"=>true, "reconfigurable" => false, "ideal" => false, "dc_link" => true, "induction_motor" => true)

function modify_sourcebus_voltage!(data_math; ground_all_buses=false, feeder="1")
    for (i,bus) in data_math["bus"]
        if bus["vbase"] == 0.2
            bus["vbase"] = 0.2309
        end
        if ground_all_buses
            if bus["grounded"] == [0, 0, 0, 0]
                bus["grounded"][end] = true
            end
        end
    end

    sourcebuses = [parse(Int,i) for (i,bus) in data_math["bus"] if bus["bus_type"]==3]

    data_math["bus"]["$(sourcebuses[1])"]["vm"] = [1.0918, 1.0445, 1.0445, 0]  # [1.1193025, 1.1193025, 1.069328, 0, 0]
    data_math["bus"]["$(sourcebuses[1])"]["va"] = [0, -121.511, 121.511, 0] .* pi/180  # [-28.533731, -151.466269, 90, 0, 0] .* pi/180
    data_math["bus"]["$(sourcebuses[1])"]["vmin"] = copy(data_math["bus"]["$(sourcebuses[1])"]["vm"])
    data_math["bus"]["$(sourcebuses[1])"]["vmax"] = copy(data_math["bus"]["$(sourcebuses[1])"]["vm"])

    data_math["bus"]["$(sourcebuses[2])"]["terminals"] = [1,2,3,4]
    data_math["gen"]["2"]["connections"] = [1,2,3,4]
    data_math["branch"]["192"]["t_connections"] = [1,2,3,4]
    data_math["branch"]["192"]["f_connections"] = [1,2,3,4]
end

function add_induction_motors!(data_math)
    sbase = data_math["settings"]["sbase"]                          # p.u.
    sbace_factor = data_math["settings"]["power_scale_factor"]      # 
    vbase = [v for v in values(data_math["settings"]["vbases_default"])][1]
    vbase_factor = data_math["settings"]["voltage_scale_factor"]
    Ibase = (sbase * sbace_factor) / (vbase * vbase_factor)  #[kA]
    zbase = (vbase * vbase_factor)^2 / (sbase * sbace_factor)
    
    IM1_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("1_276", bus["name"])][1] #"F1_882.1.2.3.4"
    IM2_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("1_556", bus["name"])][1] #"F1_882.1.2.3.4"
    IM3_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("1_899", bus["name"])][1] #"F1_882.1.2.3.4"
    
    for bus_id in [IM1_bus, IM2_bus, IM3_bus]
        load_id = [i for (i,load) in data_math["load"] if load["load_bus"]==bus_id][1]
        # load_id = length(data_math["load"]) + 1
        # data_math["load"]["$load_id"] = deepcopy(data_math["load"]["1"])
        pd = bus_id == IM3_bus ?  8*[1e3, 1e3, 1e3] : 4*[1e3, 1e3, 1e3]
        data_math["load"]["$load_id"]["connections"] = [1, 2, 3, 4]
        data_math["load"]["$load_id"]["vbase"] = 0.4
        data_math["load"]["$load_id"]["index"] = load_id
        data_math["load"]["$load_id"]["load_bus"] = bus_id
        data_math["load"]["$load_id"]["name"] = "IM_$bus_id"
        data_math["load"]["$load_id"]["source_id"] = "load.IM_$bus_id"
        data_math["load"]["$load_id"]["pd"] = pd / (sbase * sbace_factor)
        data_math["load"]["$load_id"]["pf"] = 0.85
        data_math["load"]["$load_id"]["qd"] = data_math["load"]["$load_id"]["pd"] * tan(acos(data_math["load"]["$load_id"]["pf"]))
    end

    return data_math, IM1_bus, IM2_bus, IM3_bus
end

# function add_voltage_source_feeder!(data_math)
#     # sourcebus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("sourcebus", bus["name"])]
#     sourcebus = [parse(Int,i) for (i,bus) in data_math["bus"] if bus["bus_type"]==3][1]
#     sourcebus_branches = [i for (i, branch) in data_math["branch"] if branch["f_bus"]==sourcebus]
#     sourceZ_branch_id = sourcebus_branches[1]
#     feeder2_branch_id = sourcebus_branches[2]

#     tx_branch =  [i for (i, branch) in data_math["branch"] if occursin("tr", branch["name"])]

#     sourcebus_new = 10000
#     data_math["bus"]["$sourcebus_new"] = deepcopy(data_math["bus"]["$sourcebus"])
#     data_math["bus"]["$sourcebus_new"]["bus_i"] = sourcebus_new
#     data_math["bus"]["$sourcebus_new"]["index"] = sourcebus_new
#     data_math["bus"]["$sourcebus_new"]["source_id"] = "bus.sourcebus2"
#     data_math["bus"]["$sourcebus_new"]["name"] = "sourcebus_2"
#     data_math["bus"]["$sourcebus_new"]["vm"] = [1.02, 1.02, 1.02, 0, 0]
#     data_math["bus"]["$sourcebus_new"]["va"] = [0, -120, 120, 0, 0] .* pi/180
#     data_math["bus"]["$sourcebus_new"]["vmin"] = copy(data_math["bus"]["$sourcebus_new"]["vm"])
#     data_math["bus"]["$sourcebus_new"]["vmax"] = copy(data_math["bus"]["$sourcebus_new"]["vm"])


#     data_math["branch"]["$feeder2_branch_id"]["f_bus"] = sourcebus_new # = deepcopy(data_math["branch"]["$sourceZ_branch_id"])
#     # data_math["branch"]["$feeder2_branch_id"]["vbase"] = 0.2309
#     data_math["branch"]["$feeder2_branch_id"]["name"] = "source2_Z"
#     data_math["branch"]["$feeder2_branch_id"]["br_r"] = data_math["branch"]["$sourceZ_branch_id"]["br_r"]
#     data_math["branch"]["$feeder2_branch_id"]["br_x"] = data_math["branch"]["$sourceZ_branch_id"]["br_x"]
#     data_math["branch"]["$feeder2_branch_id"]["b_fr"] = data_math["branch"]["$sourceZ_branch_id"]["b_fr"]
#     data_math["branch"]["$feeder2_branch_id"]["b_to"] = data_math["branch"]["$sourceZ_branch_id"]["b_to"]
#     data_math["branch"]["$feeder2_branch_id"]["g_fr"] = data_math["branch"]["$sourceZ_branch_id"]["g_fr"]
#     data_math["branch"]["$feeder2_branch_id"]["g_to"] = data_math["branch"]["$sourceZ_branch_id"]["g_to"]

#     data_math["gen"]["2"] = deepcopy(data_math["gen"]["1"])
#     data_math["gen"]["2"]["gen_bus"] = sourcebus_new
#     data_math["gen"]["2"]["index"] = 2
# end

##

# # data_path = "./data/ENWL_4w_Network1_Feeders1and2/Master.dss"
# data_path = "./data/feeder_12/Master.dss"
# data_path = "./data/feeder_12/Master_updated_v1.dss"

# # data_path = "./data/ENWL_4w_Network1_Feeder1/Master.dss"
# data_path = "./data/ENWL_4w_Network1_Feeders1and2 copy/Master.dss"
data_path = "./data/ENWL_4w_Network1_Feeders1and2 copy/Master1.dss"
# data_path = "./data/ENWL_4w_Network1_Feeders1and2 copy/Master2.dss"
data_path = "./data/ENWL_4w_Network1_Feeders1and2 copy/Master_updated.dss"

###
### parse data
data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!, PMD.reduce_lines!])
data_eng["settings"]["sbase_default"] = 1
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)

modify_sourcebus_voltage!(data_math; ground_all_buses=false)
data_math, IM1_bus, IM2_bus, IM3_bus = add_induction_motors!(data_math)
# add_voltage_source_feeder!(data_math)

# f1_1_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if bus["name"]=="f1_1"][1]
# data_math["bus"]["$f1_1_bus"]["grounded"][4] = 1

## ##################### No SOP #####################
PMD.add_start_vrvi!(data_math)
# model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, PMD.build_mc_opf);
model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, RPMD.build_mc_opf_mx_sop; setting=setting);
result = PMD.optimize_model!(model, optimizer=ipopt_solver)

IM1_bus_seq = abs.(RPMD.sequence(result["solution"]["bus"]["$IM1_bus"]["vr"][1:3] .+ im*result["solution"]["bus"]["$IM1_bus"]["vi"][1:3])) .* 100
IM2_bus_seq = abs.(RPMD.sequence(result["solution"]["bus"]["$IM2_bus"]["vr"][1:3] .+ im*result["solution"]["bus"]["$IM2_bus"]["vi"][1:3])) .* 100
IM3_bus_seq = abs.(RPMD.sequence(result["solution"]["bus"]["$IM3_bus"]["vr"][1:3] .+ im*result["solution"]["bus"]["$IM3_bus"]["vi"][1:3])) .* 100
[IM1_bus_seq[3] ; IM2_bus_seq[3] ; IM3_bus_seq[3]]


## ##################### Conventional SOP #####################
data_math_sop = deepcopy(data_math)
fbus = [parse(Int,i) for (i,bus) in data_math_sop["bus"] if occursin("1_882", bus["name"])][1] #"F1_882.1.2.3.4"  25?
tbus = [parse(Int,i) for (i,bus) in data_math_sop["bus"] if occursin("2_396", bus["name"])][1] #"F2_396.1.2.3.4"  123?
# data_math_sop["bus"]["$fbus"]["grounded"][4] = 1
# data_math_sop["bus"]["$tbus"]["grounded"][4] = 1
RPMD.add_sop_inverter_losses_v2!(data_math_sop, fbus, tbus; c_rating_a=25*ones(3), dc_link=setting["dc_link"])

load_id1 = [i for (i,load) in data_math_sop["load"] if load["load_bus"]==fbus][1]
load_id2 = [i for (i,load) in data_math_sop["load"] if load["load_bus"]==tbus][1]
delete!(data_math_sop["load"], load_id1)
delete!(data_math_sop["load"], load_id2)

PMD.add_start_vrvi!(data_math_sop)
model_sop = PMD.instantiate_mc_model(data_math_sop, PMD.IVRENPowerModel, RPMD.build_mc_opf_mx_sop; setting=setting);
result_conv = PMD.optimize_model!(model_sop, optimizer=ipopt_solver)

IMbus_vmneg = [result_conv["solution"]["bus"]["$bus_id"]["vmneg"] for bus_id in [IM1_bus, IM2_bus, IM3_bus]] .*100
# IM1_bus_seq_sop = abs.(RPMD.sequence(result_conv["solution"]["bus"]["$IM1_bus"]["vr"][1:3] .+ im*result_conv["solution"]["bus"]["$IM1_bus"]["vi"][1:3])) .* 100
# IM2_bus_seq_sop = abs.(RPMD.sequence(result_conv["solution"]["bus"]["$IM2_bus"]["vr"][1:3] .+ im*result_conv["solution"]["bus"]["$IM2_bus"]["vi"][1:3])) .* 100
# IM3_bus_seq_sop = abs.(RPMD.sequence(result_conv["solution"]["bus"]["$IM3_bus"]["vr"][1:3] .+ im*result_conv["solution"]["bus"]["$IM3_bus"]["vi"][1:3])) .* 100

IMbus_vmneg_nosop = [IM1_bus_seq[3] ; IM2_bus_seq[3] ; IM3_bus_seq[3]]


##
using Plots

a0 = 0.033125*100
a1 = 2.75*100
a2 = 56.25*100
Db_curve(vmneg, vmnegsqr) = vmneg <= 0.01 ? 100.0 : 
        ((vmneg >= 0.01 && vmneg <= 0.05) ? 100.0 - (a2 * vmnegsqr + a1 * vmneg - a0) : 
        0.0)

IM_load_ids = [(load["load_bus"], sum(sqrt.(load["pd"].^2 .+ load["qd"].^2))) for (i, load) in data_math_sop["load"] if startswith(load["name"], "IM")]

induction_obj = [sd * (100 - Db_curve(result_conv["solution"]["bus"]["$bus_id"]["vmneg"], result_conv["solution"]["bus"]["$bus_id"]["vmnegsqr"])) for (bus_id, sd) in IM_load_ids]

IMbus_vmneg ./= 100
IMbus_vmneg_nosop ./= 100


vmneg = collect(0.0:0.001:0.05)
plot(vmneg, Db_curve.(vmneg, vmneg.^2), title="IM derating factors", label="")
plot!([IMbus_vmneg[1]], [Db_curve.(IMbus_vmneg[1], IMbus_vmneg[1].^2)], seriestype=:scatter, lw=2, label="IM1 w SOP", color=1)
plot!([IMbus_vmneg[2]], [Db_curve.(IMbus_vmneg[2], IMbus_vmneg[2].^2)], seriestype=:scatter, lw=2, label="IM2 w SOP", color=2)
plot!([IMbus_vmneg[3]], [Db_curve.(IMbus_vmneg[3], IMbus_vmneg[3].^2)], seriestype=:scatter, lw=2, label="IM3 w SOP", color=3)

plot!([IMbus_vmneg_nosop[1]], [Db_curve.(IMbus_vmneg_nosop[1], IMbus_vmneg_nosop[1].^2)], seriestype=:scatter, lw=2, label="IM1 w/o SOP", color=1)
plot!([IMbus_vmneg_nosop[2]], [Db_curve.(IMbus_vmneg_nosop[2], IMbus_vmneg_nosop[2].^2)], seriestype=:scatter, lw=2, label="IM2 w/o SOP", color=2)
plot!([IMbus_vmneg_nosop[3]], [Db_curve.(IMbus_vmneg_nosop[3], IMbus_vmneg_nosop[3].^2)], seriestype=:scatter, lw=2, label="IM3 w/o SOP", color=3)



IM_load_ids = [(load["load_bus"], sum(sqrt.(load["pd"].^2 .+ load["qd"].^2))) 
    for (i, load) in PMD.ref(pm, 0, :load) if startswith(load["name"], "IM")]
# @show IM_load_ids

induction_obj = JuMP.@expression(pm.model,   
    sum(sd * (100 - Db(PMD.var(pm, 0, :vmneg)[bus_id], PMD.var(pm, 0, :vmnegsqr)[bus_id])) 
        for (bus_id, sd) in IM_load_ids)
    )


##
using Plots
using DataFrames
using PlotlyJS

vmneg_kron_nosop= [2.809919103201267 ; 2.7841868755472556 ; 2.7217788833144763]
vmneg_4wire_nosop = [2.799275615927309 ; 2.7728415848857093 ; 2.7089820744665585]
vmneg_kron_25Asop = [1.4085789302562775 ; 1.2122047070786937 ; 0.9739545517667308]
vmneg_4wire_25Asop_grounded = [1.4007101992333153 ; 1.2034357267112943 ; 0.9697191178395723]
vmneg_4wire_25Asop_noneutral = [1.3961531961902303 ; 1.1984924060969946 ; 0.9708437850973898]
vmneg_4wire_25Asop_ungrounded = [1.3188554627946196 ; 1.118154085967026 ; 0.9757769579360701]

vmneg_4wire_20Asop_ungrounded = [1.596146442694053 ; 1.4101093907976945 ; 0.9999910163340194]

vmneg_4wire_20Asop_grounded = []
vmneg_4wire_20Asop_nonuetral = []


group_labels = ["IM1 276 (12 kW)", "IM2 556 (12 kW)", "IM3 899 (24 kW)"]
trace1 = PlotlyJS.bar(;x=group_labels, y=vmneg_kron_nosop, name="Before: Kron")
trace2 = PlotlyJS.bar(;x=group_labels, y=vmneg_4wire_nosop, name="Before: 4wire")
trace3 = PlotlyJS.bar(;x=group_labels, y=vmneg_kron_25Asop, name="After: 25A sop Kron")
trace4 = PlotlyJS.bar(;x=group_labels, y=vmneg_4wire_25Asop_grounded, name="After: 4wire 25A sop grounded")
trace5 = PlotlyJS.bar(;x=group_labels, y=vmneg_4wire_25Asop_ungrounded, name="After: 4wire 25A sop ungrounded")
data = [trace1, trace2, trace3, trace4, trace5]
layout = Layout(yaxis_title="Voltage unbalance %";barmode="group")
vmneg_plt = PlotlyJS.plot(data, layout)
PlotlyJS.savefig(vmneg_plt, "vmneg_plot.html")


Sd = [14.12 ; 14.12 ; 28.24] # kVA
Pd = [12, 12, 24] # kW
derating_after = Pd .*  (1 .- Db_curve.(vmneg_4wire_25Asop_ungrounded./100, (vmneg_4wire_25Asop_ungrounded./100).^2)./100)
derating_before = Pd .* (1 .- Db_curve.(vmneg_kron_nosop./100, (vmneg_kron_nosop./100).^2)./100)
group_labels = ["IM1 276 (12 kW)", "IM2 556 (12 kW)", "IM3 899 (24 kW)"]
trace1 = PlotlyJS.bar(;x=group_labels, y=derating_after, name="Before: 4wire")
trace2 = PlotlyJS.bar(;x=group_labels, y=derating_before, name="Before: 4wire")
data = [trace1, trace2]
layout = Layout(yaxis_title="Required derating (kW)"; barmode="group")
derating_plt = PlotlyJS.plot(data, layout)
PlotlyJS.savefig(derating_plt, "derating_plot.html")
## ############################################################
###############################################################
data_path = "./data/ENWL_4w_Network1_Feeder1/Master.dss"
data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
data_eng["settings"]["sbase_default"] = 1
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)

sourcebus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("sourcebus", bus["name"])][1]
data_math["bus"]["$sourcebus"]["vm"] = [1.0918, 1.0445, 1.0445, 0]
data_math["bus"]["$sourcebus"]["vmin"] = copy(data_math["bus"]["$sourcebus"]["vm"])
data_math["bus"]["$sourcebus"]["vmax"] = copy(data_math["bus"]["$sourcebus"]["vm"])
data_math["bus"]["$sourcebus"]["va"] = [0, -121.511, 121.511, 0] .* pi/180

IM1_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("276", bus["name"])][1] #"F1_882.1.2.3.4"
IM2_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("556", bus["name"])][1] #"F1_882.1.2.3.4"
IM3_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("899", bus["name"])][1] #"F1_882.1.2.3.4"

sbase = data_math["settings"]["sbase"]                          # p.u.
sbace_factor = data_math["settings"]["power_scale_factor"]      # 
vbase = [v for v in values(data_math["settings"]["vbases_default"])][1]
vbase_factor = data_math["settings"]["voltage_scale_factor"]
Ibase = (sbase * sbace_factor) / (vbase * vbase_factor)  #[kA]
zbase = (vbase * vbase_factor)^2 / (sbase * sbace_factor)

for bus_id in [IM1_bus, IM2_bus, IM3_bus]
    load_id = length(data_math["load"]) + 1
    data_math["load"]["$load_id"] = deepcopy(data_math["load"]["1"])
    # load_id = [i for (i,load) in data_math["load"] if load["load_bus"] == bus_id][1]
    pd = bus_id == IM3_bus ?  8*[1e3, 1e3, 1e3] : 4*[1e3, 1e3, 1e3]
    data_math["load"]["$load_id"]["connections"] = [1, 2, 3, 4]
    data_math["load"]["$load_id"]["vbase"] = 0.4
    data_math["load"]["$load_id"]["index"] = load_id
    data_math["load"]["$load_id"]["load_bus"] = bus_id
    data_math["load"]["$load_id"]["name"] = "IM_$bus_id"
    data_math["load"]["$load_id"]["source_id"] = "load.IM_$bus_id"
    data_math["load"]["$load_id"]["pd"] = pd / (sbase * sbace_factor)
    data_math["load"]["$load_id"]["pf"] = 0.85
    data_math["load"]["$load_id"]["qd"] = data_math["load"]["$load_id"]["pd"] * tan(acos(data_math["load"]["$load_id"]["pf"]))
end

for (i, bus) in data_math["bus"]
    if length(bus["vmax"]) < 4
        @show bus
    end
end

setting = Dict("conventional"=>false, "reconfigurable" => false, "ideal" => false, "dc_link" => false, "induction_motor" => true)
PMD.add_start_vrvi!(data_math)
model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, RPMD.build_mc_opf_mx_sop; setting=setting)
result_noconv = PMD.optimize_model!(model, optimizer=ipopt_solver)

for (i, bus) in data_math["bus"]
    if length(bus["vmax"]) < 4
        @show bus
    end
end


##
model = JuMP.Model(ipopt_solver)
JuMP.@variable(model, w, lower_bound=5.0, upper_bound=10.0)
JuMP.@variable(model, x, lower_bound=5.0, upper_bound=10.0)
JuMP.@variable(model, y, lower_bound=-10, upper_bound=0.0)
JuMP.@variable(model, z, lower_bound=-10.0, upper_bound=10.0)
JuMP.@constraint(model, w + y + z == 4)  # 7 + -1 + -2 = 4
JuMP.@constraint(model, x + y + z == 2)  # 5 + -1 + -2 = 2
# JuMP.@constraint(model, w + z == 5)     # 7 + -2 = 5
# JuMP.@constraint(model, x + y == 4)     # 5 + -1 = 4

JuMP.@objective(model, Max, w^2 + x^2 + y^2 + z^2)
# JuMP.@objective(model, Max, w + x + y + z)
# JuMP.@objective(model, Min, z)
# JuMP.@objective(model, Max, z)

JuMP.optimize!(model)
@show JuMP.termination_status(model)

JuMP.value.(JuMP.all_variables(model))

