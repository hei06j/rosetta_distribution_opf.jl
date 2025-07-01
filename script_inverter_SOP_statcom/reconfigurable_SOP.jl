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

"""
TODOs:
- add DC voltage source or load
"""

# ENV["GUROBI_HOME"]="/Library/gurobi1103/macos_universal2"
# ENV["GRB_LICENSE_FILE"]="/Users/hei06j/gurobi/gurobi_11.lic"

# ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes", "hsllib"=>HSL_jll.libhsl_path, "linear_solver"=>"ma86")
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
highs_solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
gurobi_solver = optimizer_with_attributes(Gurobi.Optimizer, "output_flag" => false)
juniper_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt_solver, "mip_solver" => highs_solver)
# juniper_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt_solver, "mip_solver" => gurobi_solver)

# set_attribute(model, "hsllib", HSL_jll.libhsl_path)
# set_attribute(model, "linear_solver", "ma86")

# data_path = "./data/case5_gen_3ph_wye_v1.dss"
# dss_includes_gens = true
# sop_branch_id = "6"

data_path = "./data/ENWL_4w_Network1_Feeders1and2/Master.dss"
# data_path = "./data/European_LV_network/Master.dss"
# data_path = "./data/ENWL_4w_Network1_Feeder1/Master.dss"
dss_includes_gens = false
sop_branch_id = "907"

function add_solar_gen!(data_math, dss_includes_gens; counter=20)
    if !dss_includes_gens
        sop_gen_ids = []
        for i = 1:counter:length(data_math["load"])
            load = data_math["load"]["$i"]
            pd = load["pd"][1]
            gen_id = length(data_math["gen"]) + 1
            data_math["gen"]["$gen_id"] = deepcopy(data_math["gen"]["1"])
            Smax = 2 * ceil(pd) * ones(3)
            gen = data_math["gen"]["$gen_id"]
            gen["gen_bus"] = copy(load["load_bus"])
            gen["index"] = gen_id
            gen["smax"] = Smax
            gen["pmax"] = Smax
            gen["pmin"] = 0.0 * Smax
            gen["qmax"] = Smax #sqrt.(Smax.^2 - gen["pmax"].^2)
            gen["qmin"] = -Smax #-gen["qmax"]
            gen["cost"] = [10 0]
            gen["type"] = "GFL-4w"
            gen["name"] = "GFL-4w-bus-$gen_id"
            gen["cost"] = [1, 0]
            push!(sop_gen_ids, gen_id)
        end
        data_math["gen"]["1"]["cost"][1] = 1000
    else
        sop_gen_ids = [i for (i,gen) in data_math["gen"] if !occursin("voltage_source", gen["name"])][1:2]
        source_gen_id = [i for (i,gen) in data_math["gen"] if occursin("voltage_source", gen["name"])][1]
        data_math["gen"]["$source_gen_id"]["cost"][1] = 1000
    end

    return sop_gen_ids
end


function build_sop_case(data_eng, setting)

    @assert setting["conventional"] + setting["reconfigurable"] + setting["ideal"] == 1   "Choose only one type of SOP: conventional, reconfigurable, or ideal"

    ### transform data_eng to data_math
    data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)

    ### if the case does not have two pv gens to replace with a sop, add two pv gens to the network
    sop_gen_ids = add_solar_gen!(data_math, dss_includes_gens)

    ### Plot network graph
    # bus_coordinates_file = "./data/ENWL_4w_Network1_Feeder1/BusCoords.txt"
    # RPMD.plot_network_graph(graph, data_eng, data_math)

    ### replace the pv gens by a sop, and add the branch
    RPMD.add_sop_inverter_losses!(data_math, sop_gen_ids[1], sop_gen_ids[2]; reconfigurable=setting["reconfigurable"])

    return data_math
end


function run_sop_case(data_math, setting)
    ### build optimisation model and solve opf
    PMD.add_start_vrvi!(data_math)
    model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, RPMD.build_mc_opf_mx_sop; setting=setting)
    # result = PMD.optimize_model!(model, optimizer=juniper_solver)
    # results_mx_dict = RPMD.get_solutions(model, result)
    if setting["conventional"] || setting["ideal"]
        result = PMD.optimize_model!(model, optimizer=ipopt_solver)
    elseif setting["reconfigurable"]
        result = PMD.optimize_model!(model, optimizer=juniper_solver)
    end
    return result
end


##
### parse data
data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
data_eng["settings"]["sbase_default"] = 1
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0


## ##################### Conventional inverter #####################
### 4-leg inverters: set conventional, reconfigurable and ideal true or false
setting = Dict("conventional"=>true, "reconfigurable" => false, "ideal" => false, "dc_link" => true)
data_math_conv = build_sop_case(data_eng, setting)
result_conv = run_sop_case(data_math_conv, setting)


## ##################### Reconfigurable inverter #####################
### 4-leg inverters: set conventional, reconfigurable and ideal true or false
setting = Dict("conventional"=>false, "reconfigurable" => true, "ideal" => false, "dc_link" => true)
data_math_mx = build_sop_case(data_eng, setting)
result_mx = run_sop_case(data_math_mx, setting)


## ##################### Ideal inverter #####################
### 4-leg inverters: set conventional, reconfigurable and ideal true or false
setting = Dict("conventional"=>false, "reconfigurable" => false, "ideal" => true, "dc_link" => true)
data_math_ideal = build_sop_case(data_eng, setting)
result_ideal = run_sop_case(data_math_ideal, setting)


## inspect results
result_conv["solution"]["branch"]["$sop_branch_id"]["pdc_link"]
result_ideal["solution"]["branch"]["$sop_branch_id"]["pdc_link"]
result_mx["solution"]["branch"]["$sop_branch_id"]["pdc_link"]

round.(result_mx["solution"]["branch"]["$sop_branch_id"]["bg"])
result_conv["solution"]["branch"]["$sop_branch_id"]["cr_fr"]
result_conv["solution"]["branch"]["$sop_branch_id"]["cr_to"]
result_mx["solution"]["branch"]["$sop_branch_id"]["cr_fr"] .+ result_mx["solution"]["branch"]["$sop_branch_id"]["cr_to"]

result_mx["solution"]["branch"]["$sop_branch_id"]["pf_idx"]
result_mx["solution"]["branch"]["$sop_branch_id"]["pf"]
sum(result_mx["solution"]["branch"]["$sop_branch_id"]["pf_idx"])
sum(result_mx["solution"]["branch"]["$sop_branch_id"]["pf"])

result_mx["solution"]["branch"]["$sop_branch_id"]["pt_idx"]
result_mx["solution"]["branch"]["$sop_branch_id"]["pt"]
sum(result_mx["solution"]["branch"]["$sop_branch_id"]["pt_idx"])
sum(result_mx["solution"]["branch"]["$sop_branch_id"]["pt"])

sum(result_mx["solution"]["branch"]["$sop_branch_id"]["pf_idx"]) + sum(result_mx["solution"]["branch"]["$sop_branch_id"]["pt_idx"])
sum(result_mx["solution"]["branch"]["$sop_branch_id"]["pf"]) + sum(result_mx["solution"]["branch"]["$sop_branch_id"]["pt"])


##
# # result_mx_sols = RPMD.get_solutions(model, result_mx)
# result_mx["solution"]["gen"]["1"]
# round.(abs.(result_mx_sols["c_source_branch_012"]), digits=4)
# round.(result_mx["solution"]["branch"]["6"]["c_rating"], digits=4)
# round.(abs.(result_mx["solution"]["branch"]["6"]["cr_fr"] .+ im*result_mx["solution"]["branch"]["6"]["ci_fr"]), digits=4)

# result_mx["solution"]["branch"]["6"]
# Int.(round.(result_mx["solution"]["branch"]["7"]["bg"].+1e-6))
# # result_mx["solution"]["gen"]["2"]["pg"]


# result_conv_sols = RPMD.get_solutions(model, result_conv)
# result_conv["solution"]["gen"]["1"]
# round.(abs.(result_conv_sols["c_source_branch_012"]), digits=4)
# round.(data_math_conv["branch"]["6"]["c_rating_a"], digits=4)
# round.(abs.(result_conv["solution"]["branch"]["6"]["cr_fr"] .+ im*result_conv["solution"]["branch"]["6"]["ci_fr"]), digits=4)
# result_conv["solution"]["branch"]["6"]
# Int.(round.(result_conv["solution"]["branch"]["7"]["bg"].+1e-6))
# # result_conv["solution"]["gen"]["2"]["pg"]


# result_ideal_sols = RPMD.get_solutions(model, result_ideal)
# result_ideal["solution"]["gen"]["1"]
# round.(abs.(result_ideal_sols["c_source_branch_012"]), digits=4)
# round.(data_math_ideal["branch"]["6"]["c_rating_a"], digits=4)
# round.(abs.(result_ideal["solution"]["branch"]["6"]["cr_fr"] .+ im*result_ideal["solution"]["branch"]["6"]["ci_fr"]), digits=4)
# result_ideal["solution"]["branch"]["6"]
# Int.(round.(result_ideal["solution"]["branch"]["7"]["bg"].+1e-6))
# # result_ideal["solution"]["gen"]["2"]["pg"]

##

# ### load, inverters, and source currents
# results_GFL_4w_mx["c_load"]
# results_GFL_4w_mx["c_inverter_branch"]
# results_GFL_4w_mx["c_source_branch"]

# ### load, inverters, and source currents sequence values in complex
# results_GFL_4w_mx["c_load_012"]
# results_GFL_4w_mx["c_inverter_branch_012"]
# results_GFL_4w_mx["c_source_branch_012"]

# ### load, inverters, and source currents sequence values in magnitude
# round.(abs.(results_GFL_4w_mx["c_load_012"]), digits=6)
# round.(abs.(results_GFL_4w_mx["c_inverter_branch_012"]), digits=6)
# round.(abs.(results_GFL_4w_mx["c_source_branch_012"]), digits=6)