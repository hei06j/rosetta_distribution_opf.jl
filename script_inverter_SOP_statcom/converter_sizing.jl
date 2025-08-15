using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
# using Juniper
# using HiGHS
using JuMP  # bl/array_nl
import LinearAlgebra: diag, diagm
using LaTeXStrings
using DataFrames
using CSV
using Plots

const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

PMD.silence!()

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes", "warm_start_init_point"=>"yes", "max_iter"=>100000)

data_path = "./data/inverter_4w_wye_unbalanced_loads_2bus.dss"


function build_inverter_case(data_eng, setting)
    ### transform data_eng to data_math
    data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)

    ### add inverter lossy branches
    pv_gen_ids = [i for (i, gen) in data_math["gen"] if !occursin("source", gen["name"])]
    for gen_id in pv_gen_ids
        RPMD.add_inverter_losses!(data_math, gen_id; c_rating_a=30*ones(3))
        data_math["gen"][gen_id]["srating_min"] = 0.0
        data_math["gen"][gen_id]["srating_max"] = 10
    end

    return data_math
end


function build_data_math(data_path; sbase=0.8, setting=nothing)
    data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
    data_eng["settings"]["sbase_default"] = sbase
    data_eng["voltage_source"]["source"]["rs"] *= 0
    data_eng["voltage_source"]["source"]["xs"] *= 0

    data_math = build_inverter_case(data_eng, setting);
    data_math["gen"]["1"]["qmax"] = copy(data_math["gen"]["1"]["pmax"])
    data_math["gen"]["1"]["qmin"] = - copy(data_math["gen"]["1"]["qmax"])
    data_math["gen"]["1"]["pmax"] .= 0             # changing the inverter source to 0, not to provide any power

    for (i,bus) in data_math["bus"]
        bus["vmin"] = [0.9, 0.9, 0.9, 0] * 0.1
        bus["vmax"] = [1.1, 1.1, 1.1, 1.1] * 1.1
    end

    ### make voltage source capacity larger
    for (i, gen) in data_math["gen"]
        if gen["gen_bus"] ∈ [bus["index"] for (i, bus) in data_math["bus"] if bus["bus_type"] == 3]
            gen["pmax"] = [100, 100, 100]
            gen["pmin"] = -[100, 100, 100]
            gen["qmax"] = [100, 100, 100]
            gen["qmin"] = -[100, 100, 100]
        end
    end
    
    # Sbase = data_math["settings"]["sbase"]                        # p.u.
    sbace_factor = data_math["settings"]["power_scale_factor"]      # 
    vbase = [v for v in values(data_math["settings"]["vbases_default"])][1]
    # vbase = 0.2309      # [kV]  data_math["settings"]["vbases_default"]["5"]
    vbase_factor = data_math["settings"]["voltage_scale_factor"]
    Ibase = (sbase * sbace_factor) / (vbase * vbase_factor)  #[kA]
    
    return data_math, Ibase
end


## ##################### Conventional inverter #####################
setting = Dict("dc_link" => true)

sbase = 1
data_math, Ibase = build_data_math(data_path; setting=setting, sbase=sbase)
data_math["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 100] / Ibase  # inverter branch conductor ratings

PMD.add_start_vrvi!(data_math)
model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, RPMD.build_mc_opf_sizing; setting=setting)
result = PMD.optimize_model!(model, optimizer=ipopt_solver)


result["solution"]["gen"]["1"]

## TODOs
"""
    - make multi-period
    - add different objective functions 
    - how to include batteries?
    - scenarios
    - (how to deal with the converter branch? We should do the study with a lossless converter)
    - 
"""