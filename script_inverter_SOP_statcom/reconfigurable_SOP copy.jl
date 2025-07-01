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


# function add_solar_gen!(data_math, dss_includes_gens; counter=20)
#     if !dss_includes_gens
#         sop_gen_ids = []
#         for i = 1:counter:length(data_math["load"])
#             load = data_math["load"]["$i"]
#             pd = load["pd"][1]
#             gen_id = length(data_math["gen"]) + 1
#             data_math["gen"]["$gen_id"] = deepcopy(data_math["gen"]["1"])
#             Smax = 2 * ceil(pd) * ones(3)
#             gen = data_math["gen"]["$gen_id"]
#             gen["gen_bus"] = copy(load["load_bus"])
#             gen["index"] = gen_id
#             gen["smax"] = Smax
#             gen["pmax"] = Smax
#             gen["pmin"] = 0.0 * Smax
#             gen["qmax"] = Smax #sqrt.(Smax.^2 - gen["pmax"].^2)
#             gen["qmin"] = -Smax #-gen["qmax"]
#             gen["cost"] = [10 0]
#             gen["type"] = "GFL-4w"
#             gen["name"] = "GFL-4w-bus-$gen_id"
#             gen["cost"] = [1, 0]
#             push!(sop_gen_ids, gen_id)
#         end
#         data_math["gen"]["1"]["cost"][1] = 1000
#     else
#         sop_gen_ids = [i for (i,gen) in data_math["gen"] if !occursin("voltage_source", gen["name"])][1:2]
#         source_gen_id = [i for (i,gen) in data_math["gen"] if occursin("voltage_source", gen["name"])][1]
#         data_math["gen"]["$source_gen_id"]["cost"][1] = 1000
#     end

#     return sop_gen_ids
# end


# function build_sop_case_v2(data_eng, fbus, tbus, setting)

#     @assert setting["conventional"] + setting["reconfigurable"] + setting["ideal"] <= 1   "Choose only one type of SOP: conventional, reconfigurable, or ideal"

#     ### transform data_eng to data_math
#     data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)

#     ### if the case does not have two pv gens to replace with a sop, add two pv gens to the network
#     if setting["conventional"] || setting["reconfigurable"] || setting["ideal"]
#         ### replace the pv gens by a sop, and add the branch
#         # RPMD.add_sop_inverter_losses!(data_math, fbus, tbus reconfigurable=setting["reconfigurable"])
#         RPMD.add_sop_inverter_losses_v2!(data_math, fbus, tbus; reconfigurable=setting["conventional"], dc_link=setting["dc_link"])
#     end
    
#     ### Plot network graph
#     # bus_coordinates_file = "./data/ENWL_4w_Network1_Feeder1/BusCoords.txt"
#     # RPMD.plot_network_graph(graph, data_eng, data_math)

#     return data_math
# end


function run_sop_case(data_math, setting)
    ### build optimisation model and solve opf
    PMD.add_start_vrvi!(data_math)
    model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, RPMD.build_mc_opf_mx_sop; setting=setting)
    # result = PMD.optimize_model!(model, optimizer=juniper_solver)
    # results_mx_dict = RPMD.get_solutions(model, result)
    if setting["reconfigurable"]
        result = PMD.optimize_model!(model, optimizer=juniper_solver)
    else
        result = PMD.optimize_model!(model, optimizer=ipopt_solver)
    end
    return result
end


###
### parse data
data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
data_eng["settings"]["sbase_default"] = 1
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)

### ##################### No inverter #####################
# # data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)
# PMD.add_start_vrvi!(data_math)
# model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, PMD.build_mc_opf)
# result = PMD.optimize_model!(model, optimizer=ipopt_solver)


### ##################### Conventional inverter #####################
### 4-leg inverters: set conventional, reconfigurable and ideal true or false
data_math_conv = deepcopy(data_math)
setting = Dict("conventional"=>true, "reconfigurable" => false, "ideal" => false, "dc_link" => true, "induction_motor" => true)
fbus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("f1_882", bus["name"])][1] #"F1_882.1.2.3.4"
tbus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("f2_396", bus["name"])][1] #"F2_396.1.2.3.4"
RPMD.add_sop_inverter_losses_v2!(data_math_conv, fbus, tbus; c_rating_a=20*ones(3), dc_link=setting["dc_link"])

# edit vsource.source pu=1.0918 angle=0 phases=1
# new vsource.ENWL_network_1_Feeder_1_4wire_phaseB bus1=sourcebus.2 BasekV=6.35085 pu=1.0445 angle=-121.511 phases=1 ISC3=100000 ISC1=100000
# new vsource.ENWL_network_1_Feeder_1_4wire_phaseC bus1=sourcebus.3 BasekV=6.35085 pu=1.0445 angle=121.511 phases=1 ISC3=100000 ISC1=100000

sourcebus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("sourcebus", bus["name"])][1]
data_math["bus"]["$sourcebus"]["vm"] = [1, 1.0445, 1.0445]
data_math["bus"]["$sourcebus"]["vmin"] = copy(data_math["bus"]["$sourcebus"]["vm"])
data_math["bus"]["$sourcebus"]["vmax"] = copy(data_math["bus"]["$sourcebus"]["vm"])
data_math["bus"]["$sourcebus"]["va"] = [0, -121.511, 121.511] .* pi/180

IM1_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("f1_276", bus["name"])][1] #"F1_882.1.2.3.4"
IM2_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("f1_556", bus["name"])][1] #"F1_882.1.2.3.4"
IM3_bus = [parse(Int,i) for (i,bus) in data_math["bus"] if occursin("f1_899", bus["name"])][1] #"F1_882.1.2.3.4"

sbase = data_math["settings"]["sbase"]                          # p.u.
sbace_factor = data_math["settings"]["power_scale_factor"]      # 
vbase = [v for v in values(data_math["settings"]["vbases_default"])][1]
vbase_factor = data_math["settings"]["voltage_scale_factor"]
# vbase = 0.2309      # [kV]  data_math["settings"]["vbases_default"]["5"]
Ibase = (sbase * sbace_factor) / (vbase * vbase_factor)  #[kA]
zbase = (vbase * vbase_factor)^2 / (sbase * sbace_factor)

for bus_id in [IM1_bus, IM2_bus, IM3_bus]
    load_id = length(data_math_conv["load"]) + 1
    data_math_conv["load"]["$load_id"] = deepcopy(data_math_conv["load"]["1"])
    # load_id = [i for (i,load) in data_math_conv["load"] if load["load_bus"] == bus_id][1]
    pd = bus_id == IM3_bus ?  8*[1e3, 1e3, 1e3] : 4*[1e3, 1e3, 1e3]
    data_math_conv["load"]["$load_id"]["connections"] = [1, 2, 3, 4]
    data_math_conv["load"]["$load_id"]["vbase"] = 0.4
    data_math_conv["load"]["$load_id"]["index"] = load_id
    data_math_conv["load"]["$load_id"]["load_bus"] = bus_id
    data_math_conv["load"]["$load_id"]["name"] = "IM_$bus_id"
    data_math_conv["load"]["$load_id"]["source_id"] = "load.IM_$bus_id"
    data_math_conv["load"]["$load_id"]["pd"] = pd / (sbase * sbace_factor)
    data_math_conv["load"]["$load_id"]["pf"] = 0.85
    data_math_conv["load"]["$load_id"]["qd"] = data_math_conv["load"]["$load_id"]["pd"] * tan(acos(data_math_conv["load"]["$load_id"]["pf"]))
end

###
# PMD.add_start_vrvi!(data_math_conv)
# model = PMD.instantiate_mc_model(data_math_conv, PMD.IVRENPowerModel, RPMD.build_mc_opf_mx_sop; setting=setting)
# result = PMD.optimize_model!(model, optimizer=ipopt_solver)
result_conv = run_sop_case(data_math_conv, setting)

## inspect results 
[result_conv["solution"]["bus"]["$bus_id"]["vmneg"] for bus_id in [IM1_bus, IM2_bus, IM3_bus]] .* 100


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
