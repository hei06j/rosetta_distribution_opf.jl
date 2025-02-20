using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using Juniper
using HiGHS
using JuMP  # bl/array_nl
using LinearAlgebra
import LinearAlgebra: diag, diagm
using LaTeXStrings
using DataFrames
using CSV
using Plots

const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

PMD.silence!()

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer)#, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "sb"=>"yes","warm_start_init_point"=>"yes")
# ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer)
highs_solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
juniper_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt_solver, "mip_solver" => highs_solver)

data_path = "./data/inverter_4w_wye_unbalanced_loads_2bus.dss"
dss_includes_gens= true

function add_solar_gen!(data_math, dss_includes_gens; counter=100)
    if !dss_includes_gens
        pv_gen_ids = []
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
            push!(pv_gen_ids, gen_id)
        end
        data_math["gen"]["1"]["cost"][1] = 1000
    else
        pv_gen_ids = [i for (i,gen) in data_math["gen"] if !occursin("voltage_source", gen["name"])][1:2]
        source_gen_id = [i for (i,gen) in data_math["gen"] if occursin("voltage_source", gen["name"])][1]
        data_math["gen"]["$source_gen_id"]["cost"][1] = 1000
    end

    return pv_gen_ids
end


function build_inverter_case(data_eng, setting)
    @assert setting["conventional"] + setting["reconfigurable"] + setting["ideal"] == 1   "Choose only one type of inverter: conventional, reconfigurable, or ideal"

    ### transform data_eng to data_math
    data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)

    ### add inverter lossy branches
    pv_gen_ids = [i for (i, gen) in data_math["gen"] if !occursin("source", gen["name"])]
    for gen_id in pv_gen_ids
        RPMD.add_inverter_losses!(data_math, gen_id; c_rating_a=30*ones(3), reconfigurable=setting["reconfigurable"])
    end

    return data_math
end


function build_case_without_inverter(data_eng)
    ### transform data_eng to data_math
    data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)

    ### find pv generators
    pv_gen_ids = [i for (i, gen) in data_math["gen"] if !occursin("source", gen["name"])]

    ### remove pv generators
    for id in pv_gen_ids
        delete!(data_math["gen"], id)
    end

    return data_math
end


function run_inverter_case(data_math, setting)
    ### build optimisation model and solve opf
    PMD.add_start_vrvi!(data_math)
    model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, RPMD.build_mc_opf_mx; setting=setting)

    if setting["conventional"] || setting["ideal"]
        result = PMD.optimize_model!(model, optimizer=ipopt_solver)
    elseif setting["reconfigurable"]
        result = PMD.optimize_model!(model, optimizer=juniper_solver)
    end

    # result_summary = RPMD.get_solutions(model, result)

    return result
end


function update_load_data!(data_math, load_data, timestep)
    data_math["load"]["3"]["pd"][1] = load_data[timestep, 1] * data_math["load"]["3"]["vnom_kv"]*data_math["load"]["3"]["vbase"]
    data_math["load"]["2"]["pd"][1] = load_data[timestep, 2] * data_math["load"]["2"]["vnom_kv"]*data_math["load"]["2"]["vbase"]
    data_math["load"]["1"]["pd"][1] = load_data[timestep, 3] * data_math["load"]["1"]["vnom_kv"]*data_math["load"]["1"]["vbase"]

end


function plot_load_duration_curve!(data; plt=plot(), title="Load Duration Curve", xlabel="Percentage of Time (%)", ylabel="Load (MW)", label=false)
    # Sort load data in descending order
    # sorted_data = sort(data, rev=true)
    sorted_data = sort(data[:,1])

    # Generate time percentages (normalized x-axis for duration curve)
    n = length(sorted_data)
    time_percent = (1:n) ./ n * 100  # Percentage of time

    # Plot the Load Duration Curve
    plot!(plt, time_percent, sorted_data, 
        label=label, 
        xlabel=xlabel, 
        ylabel=ylabel, 
        title=title, 
        lw=2, 
        grid=true)
end

function build_data_math(data_path, load_data, i; sbase=0.8, setting=nothing)
    data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
    data_eng["settings"]["sbase_default"] = sbase
    data_eng["voltage_source"]["source"]["rs"] *= 0
    data_eng["voltage_source"]["source"]["xs"] *= 0

    if isnothing(setting)
        @show setting
        data_math = build_case_without_inverter(data_eng);
    else
        data_math = build_inverter_case(data_eng, setting);
        data_math["gen"]["1"]["qmax"] = copy(data_math["gen"]["1"]["pmax"])
        data_math["gen"]["1"]["qmin"] = - copy(data_math["gen"]["1"]["qmax"])
        data_math["gen"]["1"]["pmax"] .= 0             # changing the inverter source to 0, not to provide any power
    end

    for (i,bus) in data_math["bus"]
        # if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
            bus["vmin"] = [0.9, 0.9, 0.9, 0]*0.9
            bus["vmax"] = [1.1, 1.1, 1.1, 1.1]*1.1
            # bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
            # bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
            # bus["grounded"] .=  0
        # end
    end

    for (i, gen) in data_math["gen"]
        if gen["gen_bus"] ∈ [bus["index"] for (i, bus) in data_math["bus"] if bus["bus_type"] == 3]
            @show gen["gen_bus"], [bus["index"] for (i, bus) in data_math["bus"] if bus["bus_type"] == 3]
            gen["pmax"] = [100, 100, 100]
            gen["pmin"] = -[100, 100, 100]
            gen["qmax"] = [100, 100, 100]
            gen["qmin"] = -[100, 100, 100]
        end
    end
    
    update_load_data!(data_math, load_data, i) ## Update load data

    # Sbase = data_math["settings"]["sbase"]                        # p.u.
    sbace_factor = data_math["settings"]["power_scale_factor"]      # 
    vbase = [v for v in values(data_math["settings"]["vbases_default"])][1]
    # vbase = 0.2309      # [kV]  data_math["settings"]["vbases_default"]["5"]
    vbase_factor = data_math["settings"]["voltage_scale_factor"]
    Ibase = (sbase * sbace_factor) / (vbase * vbase_factor)  #[kA]
    
    return data_math, Ibase
end


## Read load data
laod_data_path = "./data/Load profiles/network_example.csv"
load_data = CSV.read(laod_data_path, DataFrame, header=false)

setting = Dict("conventional"=>true, "reconfigurable" => false, "ideal" => false, "dc_link" => true)
i = 1


## sc 0
# sbase = 1
# data_math, Ibase = build_data_math(data_path, load_data, i; sbase=sbase)

### sc 1
sbase = 1
data_math, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
# data_math["gen"]["1"]["pdcmin"] = -Inf
# data_math["gen"]["1"]["pdcmax"] = Inf / sbase
data_math["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 100] / Ibase

# ### sc 2
# sbase = 1
# data_math, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
# # data_math["gen"]["1"]["pdcmin"] = 0
# # data_math["gen"]["1"]["pdcmax"] = 100 / sbase
# data_math["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 0.0] / Ibase
# # data_math["branch"]["2"]["c_rating_a"] = copy(data_math["gen"]["1"]["c_rating"])

### sc3
# data_math, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
# data_math["gen"]["1"]["pdcmin"] = 0
# data_math["gen"]["1"]["pdcmax"] = 0
# data_math["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 100] / Ibase

### sc4
# data_math, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
# data_math["gen"]["1"]["pdcmin"] = 0
# data_math["gen"]["1"]["pdcmax"] = 3 / sbase
# data_math["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 30] / Ibase

##
PMD.add_start_vrvi!(data_math)
model = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, RPMD.build_mc_opf_mx; setting=setting)
result_conv = PMD.optimize_model!(model, optimizer=ipopt_solver)

# println(model.model)
##
open("model.txt", "w") do f
    println(f, model)
 end
##
crd = [result_conv["solution"]["load"]["3"]["crd"][1] result_conv["solution"]["load"]["2"]["crd"][1] result_conv["solution"]["load"]["1"]["crd"][1]]'/data_math["load"]["3"]["vbase"]
cid = [result_conv["solution"]["load"]["3"]["cid"][1] result_conv["solution"]["load"]["2"]["cid"][1] result_conv["solution"]["load"]["1"]["cid"][1]]'/data_math["load"]["3"]["vbase"]
cd = [crd .+ im*cid ; -sum(crd)-im*sum(cid)]
cdm = abs.([crd .+ im*cid ; sum(crd) + im*sum(cid)])

vsource_crg = result_conv["solution"]["gen"]["2"]["crg"]./data_math["load"]["2"]["vbase"]
vsource_cig = result_conv["solution"]["gen"]["2"]["cig"]./data_math["load"]["2"]["vbase"]
vsource_cg = [vsource_crg.+im*vsource_cig ; -sum(vsource_crg)-im*sum(vsource_cig)]
vsource_cgm = abs.(vsource_cg)

branch_cr = result_conv["solution"]["branch"]["1"]["cr_fr"]./data_math["load"]["1"]["vbase"]
branch_ci = result_conv["solution"]["branch"]["1"]["ci_fr"]./data_math["load"]["1"]["vbase"]
branch_c = branch_cr .+ im*branch_ci
branch_cm = abs.(branch_cr .+ im*branch_ci)

inverter_crg = result_conv["solution"]["gen"]["1"]["crg"]./data_math["load"]["2"]["vbase"]
inverter_cig = result_conv["solution"]["gen"]["1"]["cig"]./data_math["load"]["2"]["vbase"]
inverter_cg = [inverter_crg.+im*inverter_cig ; -sum(inverter_crg)-im*sum(inverter_cig)]
inverter_cgm = abs.(inverter_cg)

inverter_cr = result_conv["solution"]["branch"]["2"]["cr_fr"]./data_math["load"]["2"]["vbase"]
inverter_ci = result_conv["solution"]["branch"]["2"]["ci_fr"]./data_math["load"]["2"]["vbase"]
inverter_c = inverter_cr .+ im*inverter_ci
inverter_cm = abs.(inverter_cr .+ im*inverter_ci)

# [cdm[1:3] branch_cm[1:3] inverter_cm[1:3]]
# [cdm branch_cm inverter_cm]
# [cd branch_c inverter_c]
# sum([-cd branch_c inverter_c], dims=2)

[cd vsource_cg inverter_cg]
sum([-cd vsource_cg inverter_cg], dims=2)

abs.([cd vsource_cg inverter_cg])

##
# [[crd+im*cid;-sum(crd)-im*sum(cid)] branch_cr+im*branch_ci [inverter_cr .+ im*inverter_ci ; -sum(inverter_cr) - im*sum(inverter_ci)]]
# sum([-[crd+im*cid;-sum(crd)-im*sum(cid)] branch_cr+im*branch_ci [inverter_cr .+ im*inverter_ci ; sum(inverter_cr) + im*sum(inverter_ci)]], dims=2)

I_seq_re, I_seq_im, I_seq_m_load = RPMD.get_sequence_components(crd, cid)
I_seq_re, I_seq_im, I_seq_m_source = RPMD.get_sequence_components(branch_cr[1:3], branch_ci[1:3])
I_seq_re, I_seq_im, I_seq_m_statcom = RPMD.get_sequence_components(inverter_cr[1:3], inverter_ci[1:3])
[I_seq_m_load I_seq_m_source I_seq_m_statcom]

I_seq_re, I_seq_im, I_seq_m_load = RPMD.get_sequence_components(crd, cid)
I_seq_re, I_seq_im, I_seq_m_source = RPMD.get_sequence_components(vsource_crg[1:3], vsource_cig[1:3])
I_seq_re, I_seq_im, I_seq_m_statcom = RPMD.get_sequence_components(inverter_crg[1:3], inverter_cig[1:3])
[I_seq_m_load I_seq_m_source I_seq_m_statcom]

##
result_conv["solution"]["gen"]["1"]["pdc_link"]

##
vr_1 = result_conv["solution"]["bus"]["1"]["vr"]
vi_1 = result_conv["solution"]["bus"]["1"]["vi"]
vm_1 = abs.(vr_1 + im*vi_1)

vr_2 = result_conv["solution"]["bus"]["2"]["vr"]
vi_2 = result_conv["solution"]["bus"]["2"]["vi"]
vm_2 = abs.(vr_2 + im*vi_2)

vr_3 = result_conv["solution"]["bus"]["3"]["vr"]
vi_3 = result_conv["solution"]["bus"]["3"]["vi"]
vm_3 = abs.(vr_3 + im*vi_3)

[vm_2 vm_1 vm_3]