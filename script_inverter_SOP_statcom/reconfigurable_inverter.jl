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

const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
highs_solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
juniper_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt_solver, "mip_solver" => highs_solver)

""" TODO
- assymetrical alpha's
- change the testcase
- timeseries analysis
"""

data_path = "./data/inverter_4w_wye_unbalanced_loads.dss"
# data_path = "./data/case5_gen_3ph_wye.dss"

# data_path = "./data/ENWL_4w_Network1_Feeders1and2/Master.dss"
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
        RPMD.add_inverter_losses!(data_math, gen_id; reconfigurable=setting["reconfigurable"])
        # data_math["gen"]["1"]["pdcmax"] = 1
    end
    # data_math["branch"]["2"]["c_rating_a"] = ones(4)

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

    result_summary = RPMD.get_solutions(model, result)

    return result, result_summary
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
data_math_conv = build_inverter_case(data_eng, setting)
result_conv, result_conv_summary = run_inverter_case(data_math_conv, setting)

result_conv["solution"]["gen"]["1"]
round.(abs.(result_conv_summary["c_source_branch_012"]), digits=4)
round.(data_math_conv["branch"]["2"]["c_rating_a"], digits=4)
round.(abs.(result_conv["solution"]["branch"]["2"]["cr_fr"] .+ im*result_conv["solution"]["branch"]["2"]["ci_fr"]), digits=4)

## testing whether to put the pv gen in the opendss file or add in the script
# result_conv_dss = deepcopy(result_conv)
# result_conv_fcn = deepcopy(result_conv)

# result_conv_dss["solution"]["gen"]["2"]
# result_conv_fcn["solution"]["gen"]["2"]

## ##################### Reconfigurable inverter #####################
### 4-leg inverters: set conventional, reconfigurable and ideal true or false
setting = Dict("conventional"=>false, "reconfigurable" => true, "ideal" => false, "dc_link" => true)
data_math_mx = build_inverter_case(data_eng, setting)
result_mx, result_mx_summary = run_inverter_case(data_math_mx, setting)

result_mx["solution"]["gen"]["1"]
round.(abs.(result_mx_summary["c_source_branch_012"]), digits=4)
round.(result_mx["solution"]["branch"]["2"]["c_rating"], digits=4)
round.(abs.(result_mx["solution"]["branch"]["2"]["cr_fr"] .+ im*result_mx["solution"]["branch"]["2"]["ci_fr"]), digits=4)


## ##################### Ideal inverter #####################
### 4-leg inverters: set conventional, reconfigurable and ideal true or false
setting = Dict("conventional"=>false, "reconfigurable" => false, "ideal" => true, "dc_link" => true)
data_math_ideal = build_inverter_case(data_eng, setting)
result_ideal, result_ideal_summary = run_inverter_case(data_math_ideal, setting)

result_ideal["solution"]["gen"]["1"]
round.(abs.(result_ideal_summary["c_source_branch_012"]), digits=4)
round.(sum(data_math_ideal["branch"]["2"]["c_rating_a"]), digits=4) * ones(4)
round.(abs.(result_ideal["solution"]["branch"]["2"]["cr_fr"] .+ im*result_ideal["solution"]["branch"]["2"]["ci_fr"]), digits=4)


## inspect more results from the model - Not working here as the model is not available in the script
ref_gen, ref_bus, ref_arc, ref_branch = RPMD.get_ref_bus_branch(PMD.ref(model, 0))

ref_branch_current_mx = result_mx["solution"]["branch"]["$ref_branch"]["cr_fr"] .+ im*result_mx["solution"]["branch"]["$ref_branch"]["ci_fr"]
ref_branch_current_conv = result_conv["solution"]["branch"]["$ref_branch"]["cr_fr"] .+ im*result_conv["solution"]["branch"]["$ref_branch"]["ci_fr"]
ref_branch_current_ideal = result_ideal["solution"]["branch"]["$ref_branch"]["cr_fr"] .+ im*result_ideal["solution"]["branch"]["$ref_branch"]["ci_fr"]

ref_branch_current_abs_mx = sum(abs.(ref_branch_current_mx))
ref_branch_current_abs_conv = sum(abs.(ref_branch_current_conv))
ref_branch_current_abs_ideal = sum(abs.(ref_branch_current_ideal))

ref_branch_current_012_mx = abs.(RPMD.sequence(ref_branch_current_mx[1:3]))
ref_branch_current_012_conv = abs.(RPMD.sequence(ref_branch_current_conv[1:3]))
ref_branch_current_012_ideal = abs.(RPMD.sequence(ref_branch_current_ideal[1:3]))

Imax = maximum([ref_branch_current_abs_mx ; ref_branch_current_abs_conv ; ref_branch_current_abs_ideal] )
ref_branch_current_plot = plot_phasors!(ref_branch_current_mx, Imax)
plot_phasors!(ref_branch_current_conv, Imax)
plot_phasors!(ref_branch_current_ideal, Imax)

using Plots
function plot_phasors!(phasor, Imax; labeled=false, I2=[], I0=[])
    plt = Plots.plot([0,imag.(phasor[1])], [0,real.(phasor[1])], arrow=true, color=:blue, linewidth=3, linestyle=:solid, label="a", border=:none)
    Plots.plot!([0,imag.(phasor[2])], [0,real.(phasor[2])], arrow=true, color=:red, linewidth=3, linestyle=:solid, label="b", border=:none)
    Plots.plot!([0,imag.(phasor[3])], [0,real.(phasor[3])], arrow=true, color=:green, linewidth=3, linestyle=:solid, label="c", border=:none)
    if phasor[4] !==  0 + 0im
        Plots.plot!([0,imag.(phasor[4])], [0,real.(phasor[4])], arrow=true, color=:black, linewidth=3, linestyle=:solid, label="n", border=:none)
    end
    Plots.plot!([0,0], [0,1.1*Imax], arrow=true, color=:grey, linestyle=:dot, label=false)
    Plots.plot!([0,1.1*Imax*real(exp(im*210/180*pi))], [0,1.1*Imax*imag(exp(im*210/180*pi))], arrow=true, color=:grey, linestyle=:dot, label=false)
    Plots.plot!([0,1.1*Imax*real(exp(im*330/180*pi))], [0,1.1*Imax*imag(exp(im*330/180*pi))], arrow=true, color=:grey, linestyle=:dot, label=false)
    if labeled
        Plots.plot!(Imax*exp.(im*(0:0.01:2pi)), color=:black, border=:none, label=false, markersize=10, legend=:bottom, legendcolumns=4, legendfontsize=30)
    else
        Plots.plot!(Imax*exp.(im*(0:0.01:2pi)), color=:black, border=:none, label=false, markersize=10, legend=false)
    end
    if !isempty(I2)
        annotate!([-7], [-Imax], text(latexstring("I_2= $(I2)"), :black, 40))
    end
    if !isempty(I0)
        annotate!([-7], [-Imax+4], text(latexstring("I_0= $(I0)"), :black, 40))
    end
    return plt
end

vuf_ref_branch_current = abs(ref_branch_current_012[3]) / abs(ref_branch_current_012[3]) * 100

vm = [sqrt.(bus["vr"].^2 .+ bus["vi"].^2) for (i,bus) in result_mx["solution"]["bus"]]

##

round.(value.(result["solution"]["gen"]["1"]["bg"]), digits=2)
round.(value.(result["solution"]["gen"]["2"]["bg"]), digits=2)
round.(value.(result["solution"]["gen"]["3"]["bg"]), digits=2)

results_GFL_4w_mx = RPMD.get_solutions(model, result)

### load, inverters, and source currents
results_GFL_4w_mx["c_load"]
results_GFL_4w_mx["c_inverter_branch"]
results_GFL_4w_mx["c_source_branch"]

### load, inverters, and source currents sequence values in complex
results_GFL_4w_mx["c_load_012"]
results_GFL_4w_mx["c_inverter_branch_012"]
results_GFL_4w_mx["c_source_branch_012"]

### load, inverters, and source currents sequence values in magnitude
round.(abs.(results_GFL_4w_mx["c_load_012"]), digits=6)
round.(abs.(results_GFL_4w_mx["c_inverter_branch_012"]), digits=6)
round.(abs.(results_GFL_4w_mx["c_source_branch_012"]), digits=6)