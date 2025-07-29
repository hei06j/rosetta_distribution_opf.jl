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
using JuMP
using LinearAlgebra
import LinearAlgebra: diag, diagm
using LaTeXStrings

const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

# ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes", "hsllib"=>HSL_jll.libhsl_path, "linear_solver"=>"ma86")
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")

function modify_sourcebus_voltage!(data_math)
    # for (i,bus) in data_math["bus"]
    #     bus["vmin"] .= 0.5
    #     bus["vmax"] .= 1.5
    # end
    # txbus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("tx1", bus["name"])][1]
    sourcebus = [parse(Int,i) for (i,bus) in data_math["bus"] if bus["bus_type"]==3][1]
    # sourcebus = 956
    # data_math["bus"]["$sourcebus"]["vm"] = [1.0918, 1.0445, 1.0445, 0, 0]
    data_math["bus"]["$sourcebus"]["vm"] = [1.0918, 1.0445, 1.0445, 0]
    data_math["bus"]["$sourcebus"]["va"] = [0, -121.511, 121.511, 0] .* pi/180
    # data_math["bus"]["$sourcebus"]["vm"] = [1.1193025, 1.1193025, 1.069328, 0, 0]
    # data_math["bus"]["$sourcebus"]["va"] = [-28.533731, -151.466269, 90, 0, 0] .* pi/180
    data_math["bus"]["$sourcebus"]["vmin"] = copy(data_math["bus"]["$sourcebus"]["vm"])
    data_math["bus"]["$sourcebus"]["vmax"] = copy(data_math["bus"]["$sourcebus"]["vm"])
end

function add_voltage_source_feeder!(data_math)
    # sourcebus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("sourcebus", bus["name"])]
    sourcebus = [parse(Int,i) for (i,bus) in data_math["bus"] if bus["bus_type"]==3][1]
    sourcebus_branches = [i for (i, branch) in data_math["branch"] if branch["f_bus"]==sourcebus]
    sourceZ_branch_id = sourcebus_branches[1]
    feeder2_branch_id = sourcebus_branches[2]

    tx_branch =  [i for (i, branch) in data_math["branch"] if occursin("tr", branch["name"])]

    sourcebus_new = 10000
    data_math["bus"]["$sourcebus_new"] = deepcopy(data_math["bus"]["$sourcebus"])
    data_math["bus"]["$sourcebus_new"]["bus_i"] = sourcebus_new
    data_math["bus"]["$sourcebus_new"]["index"] = sourcebus_new
    data_math["bus"]["$sourcebus_new"]["source_id"] = "bus.sourcebus2"
    data_math["bus"]["$sourcebus_new"]["name"] = "sourcebus_2"
    data_math["bus"]["$sourcebus_new"]["vm"] = [1.02, 1.02, 1.02, 0, 0]
    data_math["bus"]["$sourcebus_new"]["va"] = [0, -120, 120, 0, 0] .* pi/180
    data_math["bus"]["$sourcebus_new"]["vmin"] = copy(data_math["bus"]["$sourcebus_new"]["vm"])
    data_math["bus"]["$sourcebus_new"]["vmax"] = copy(data_math["bus"]["$sourcebus_new"]["vm"])


    data_math["branch"]["$feeder2_branch_id"]["f_bus"] = sourcebus_new # = deepcopy(data_math["branch"]["$sourceZ_branch_id"])
    # data_math["branch"]["$feeder2_branch_id"]["vbase"] = 0.2309
    data_math["branch"]["$feeder2_branch_id"]["name"] = "source2_Z"
    data_math["branch"]["$feeder2_branch_id"]["br_r"] = data_math["branch"]["$sourceZ_branch_id"]["br_r"]
    data_math["branch"]["$feeder2_branch_id"]["br_x"] = data_math["branch"]["$sourceZ_branch_id"]["br_x"]
    data_math["branch"]["$feeder2_branch_id"]["b_fr"] = data_math["branch"]["$sourceZ_branch_id"]["b_fr"]
    data_math["branch"]["$feeder2_branch_id"]["b_to"] = data_math["branch"]["$sourceZ_branch_id"]["b_to"]
    data_math["branch"]["$feeder2_branch_id"]["g_fr"] = data_math["branch"]["$sourceZ_branch_id"]["g_fr"]
    data_math["branch"]["$feeder2_branch_id"]["g_to"] = data_math["branch"]["$sourceZ_branch_id"]["g_to"]

    data_math["gen"]["2"] = deepcopy(data_math["gen"]["1"])
    data_math["gen"]["2"]["gen_bus"] = sourcebus_new
    data_math["gen"]["2"]["index"] = 2
end

function add_induction_motors!(data_math; combined=true)
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

## ######################## Feeder 1 : No SOP #####################
data_path = "./data/ENWL_4w_Network1_Feeders1and2 copy/Master1.dss"

data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!, PMD.reduce_lines!])
data_eng["settings"]["sbase_default"] = 1
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)
modify_sourcebus_voltage!(data_math)
data_math, IM1_bus, IM2_bus, IM3_bus = add_induction_motors!(data_math; combined=true)
# f1_1_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if bus["name"]=="f1_1"][1]
# data_math["bus"]["$f1_1_bus"]["grounded"][4] = 1

data_eng1 = deepcopy(data_eng)

### ##################### No SOP #####################
# setting = Dict("conventional"=>true, "reconfigurable" => false, "ideal" => false, "dc_link" => true, "induction_motor" => true)
# PMD.add_start_vrvi!(data_math)
# # model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, PMD.build_mc_opf);
# model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, RPMD.build_mc_opf_mx_sop; setting=setting);
# result = PMD.optimize_model!(model, optimizer=ipopt_solver)

# IM1_bus_seq = abs.(RPMD.sequence(result["solution"]["bus"]["$IM1_bus"]["vr"][1:3] .+ im*result["solution"]["bus"]["$IM1_bus"]["vi"][1:3])) .* 100
# IM2_bus_seq = abs.(RPMD.sequence(result["solution"]["bus"]["$IM2_bus"]["vr"][1:3] .+ im*result["solution"]["bus"]["$IM2_bus"]["vi"][1:3])) .* 100
# IM3_bus_seq = abs.(RPMD.sequence(result["solution"]["bus"]["$IM3_bus"]["vr"][1:3] .+ im*result["solution"]["bus"]["$IM3_bus"]["vi"][1:3])) .* 100
# [IM1_bus_seq[3] ; IM2_bus_seq[3] ; IM3_bus_seq[3]]


## ######################## Feeder 2 : No SOP #####################
data_path = "./data/ENWL_4w_Network1_Feeders1and2 copy/Master2.dss"

data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!, PMD.reduce_lines!])
data_eng["settings"]["sbase_default"] = 1
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)
# f1_1_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if bus["name"]=="f1_1"][1]
# data_math["bus"]["$f1_1_bus"]["grounded"][4] = 1

data_eng2 = deepcopy(data_eng)


# PMD.add_start_vrvi!(data_math)
# # model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, PMD.build_mc_opf);
# model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, RPMD.build_mc_opf_mx_sop; setting=setting);
# result = PMD.optimize_model!(model, optimizer=ipopt_solver)

##
data_eng["bus"]["sourcebus2"] = deepcopy(data_eng["bus"]["sourcebus"])
source_branch2_id = [i for (i,line) in data_eng["line"] if line["f_bus"] == "sourcebus"][1]
data_eng["line"][source_branch2_id]["f_bus"] = "sourcebus2"
data_eng["voltage_source"]["source2"] = deepcopy(data_eng["voltage_source"]["source"])



merge!(data_eng["bus"], data_eng1["bus"])
merge!(data_eng["line"], data_eng1["line"])
merge!(data_eng["voltage_source"], data_eng1["voltage_source"])
data_eng["switch"] = deepcopy(data_eng1["switch"])
data_eng["switch"]["switch_1"]["t_bus"] = "f2_396"


data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)
