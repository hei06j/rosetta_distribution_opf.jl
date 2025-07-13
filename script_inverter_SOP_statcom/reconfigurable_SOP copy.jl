using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using AppleAccelerate
# using HSL_jll
using Ipopt
using Juniper
using HiGHS
using Gurobi
using JuMP
using LinearAlgebra
import LinearAlgebra: diag, diagm
using LaTeXStrings

const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

# ENV["GUROBI_HOME"]="/Library/gurobi1103/macos_universal2"
# ENV["GRB_LICENSE_FILE"]="/Users/hei06j/gurobi/gurobi_11.lic"

# ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes", "hsllib"=>HSL_jll.libhsl_path, "linear_solver"=>"ma86")
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
highs_solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
gurobi_solver = optimizer_with_attributes(Gurobi.Optimizer, "output_flag" => false)
juniper_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt_solver, "mip_solver" => highs_solver)

data_path = "./data/ENWL_4w_Network1_Feeders1and2/Master.dss"

data_path = "./data/feeder_12/Master.dss"

# data_path = "./data/ENWL_4w_Network1_Feeder1/Master.dss"

# data_path = "./data/ENWL_4w_Network1_Feeders1and2_grounded/Master.dss"

##
### parse data
data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
data_eng["settings"]["sbase_default"] = 1
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)

sourcebus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("sourcebus", bus["name"])][1]
# sourcebus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("tx1", bus["name"])][1]
sourcebus_branches = [i for (i, branch) in data_math["branch"] if branch["f_bus"]==sourcebus]

source2branch = "657"
new_bus_id = 10000
data_math["bus"]["$new_bus_id"] = deepcopy(data_math["bus"]["$sourcebus"])
data_math["bus"]["$new_bus_id"]["bus_i"] = new_bus_id
data_math["bus"]["$new_bus_id"]["index"] = new_bus_id
data_math["bus"]["$new_bus_id"]["source_id"] = "bus.sourcebus2"
data_math["bus"]["$new_bus_id"]["name"] = "sourcebus_2"

data_math["branch"]["657"]["f_bus"] = new_bus_id

data_math["gen"]["2"] = deepcopy(data_math["gen"]["1"])
data_math["gen"]["2"]["gen_bus"] = new_bus_id
data_math["gen"]["2"]["index"] = 2


##
function make_case(data_math; combined=true)
    # sourcebus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("sourcebus", bus["name"])][1]
    # data_math["bus"]["$sourcebus"]["vm"] = [1.0918, 1.0445, 1.0445, 0, 0]
    # data_math["bus"]["$sourcebus"]["vmin"] = copy(data_math["bus"]["$sourcebus"]["vm"])
    # data_math["bus"]["$sourcebus"]["vmax"] = copy(data_math["bus"]["$sourcebus"]["vm"])
    # data_math["bus"]["$sourcebus"]["va"] = [0, -121.511, 121.511, 0, 0] .* pi/180

    # for (i, load) in data_math["load"]
    #     load["pd"] = [1]
    # end

    
    sbase = data_math["settings"]["sbase"]                          # p.u.
    sbace_factor = data_math["settings"]["power_scale_factor"]      # 
    vbase = [v for v in values(data_math["settings"]["vbases_default"])][1]
    vbase_factor = data_math["settings"]["voltage_scale_factor"]
    Ibase = (sbase * sbace_factor) / (vbase * vbase_factor)  #[kA]
    zbase = (vbase * vbase_factor)^2 / (sbase * sbace_factor)
    
    if combined
        IM1_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("1_276", bus["name"])][1] #"F1_882.1.2.3.4"
        IM2_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("1_556", bus["name"])][1] #"F1_882.1.2.3.4"
        IM3_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("1_899", bus["name"])][1] #"F1_882.1.2.3.4"
    else
        IM1_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("276", bus["name"])][1] #"F1_882.1.2.3.4"
        IM2_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("556", bus["name"])][1] #"F1_882.1.2.3.4"
        IM3_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("899", bus["name"])][1] #"F1_882.1.2.3.4"
    end
    
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

        @show data_math["bus"]["$bus_id"]["grounded"]
    end

    return data_math, IM1_bus, IM2_bus, IM3_bus
end

### ##################### No SOP #####################
setting = Dict("conventional"=>true, "reconfigurable" => false, "ideal" => false, "dc_link" => true, "induction_motor" => false)
data_math_no_sop, IM1_bus, IM2_bus, IM3_bus = make_case(deepcopy(data_math); combined=true)
PMD.add_start_vrvi!(data_math_no_sop)
# model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, PMD.build_mc_opf);
model = PMD.instantiate_mc_model(data_math_no_sop, PMD.IVRENPowerModel, RPMD.build_mc_opf_mx_sop; setting=setting);
result = PMD.optimize_model!(model, optimizer=ipopt_solver)


IM1_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("1_276", bus["name"])][1] #"F1_882.1.2.3.4"
IM2_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("1_556", bus["name"])][1] #"F1_882.1.2.3.4"
IM3_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("1_899", bus["name"])][1] #"F1_882.1.2.3.4"

IM1_bus_seq = abs.(RPMD.sequence(result["solution"]["bus"]["$IM1_bus"]["vr"][1:3] .+ im*result["solution"]["bus"]["$IM1_bus"]["vi"][1:3])) .* 100
IM2_bus_seq = abs.(RPMD.sequence(result["solution"]["bus"]["$IM2_bus"]["vr"][1:3] .+ im*result["solution"]["bus"]["$IM2_bus"]["vi"][1:3])) .* 100
IM3_bus_seq = abs.(RPMD.sequence(result["solution"]["bus"]["$IM3_bus"]["vr"][1:3] .+ im*result["solution"]["bus"]["$IM3_bus"]["vi"][1:3])) .* 100

[IM1_bus_seq[3] ; IM2_bus_seq[3] ; IM3_bus_seq[3]]

# [result["solution"]["bus"]["$bus_id"]["vmneg"] for bus_id in [IM1_bus, IM2_bus, IM3_bus]] .* 100

## ##################### Conventional SOP #####################
setting = Dict("conventional"=>true, "reconfigurable" => false, "ideal" => false, "dc_link" => true, "induction_motor" => true)
data_math_sop, IM1_bus, IM2_bus, IM3_bus = make_case(deepcopy(data_math); combined=true)
fbus = [parse(Int,i) for (i,bus) in data_math_sop["bus"] if occursin("f1_882", bus["name"])][1] #"F1_882.1.2.3.4"
tbus = [parse(Int,i) for (i,bus) in data_math_sop["bus"] if occursin("f2_396", bus["name"])][1] #"F2_396.1.2.3.4"
RPMD.add_sop_inverter_losses_v2!(data_math_sop, fbus, tbus; c_rating_a=20*ones(3), dc_link=setting["dc_link"])

###
PMD.add_start_vrvi!(data_math_sop)
model_sop = PMD.instantiate_mc_model(data_math_sop, PMD.IVRENPowerModel, RPMD.build_mc_opf_mx_sop; setting=setting);
result_conv = PMD.optimize_model!(model_sop, optimizer=ipopt_solver)

## inspect results 
IMbus_vmneg = [result_conv["solution"]["bus"]["$bus_id"]["vmneg"] for bus_id in [IM1_bus, IM2_bus, IM3_bus]] .*100

IMbus_vmneg_nosop = [IM1_bus_seq[3] ; IM2_bus_seq[3] ; IM3_bus_seq[3]]./100

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
