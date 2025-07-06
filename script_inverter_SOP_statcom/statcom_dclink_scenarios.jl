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


# ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes", "warm_start_init_point"=>"yes", "max_iter"=>100000)
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes", "warm_start_init_point"=>"yes", "max_iter"=>100000)
# ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer)
# highs_solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
# juniper_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt_solver, "mip_solver" => highs_solver)

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

function plot_load_duration_curve!(data; plt=plot(), title="Load Duration Curve", xlabel="Percentage of Time (%)", ylabel="Load (MW)", label=false, markershape=:none)
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
        # lw=2, 
        grid=true,
        markershape=markershape)
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
            bus["vmin"] = [0.9, 0.9, 0.9, 0] * 0.1
            bus["vmax"] = [1.1, 1.1, 1.1, 1.1] * 1.1
            # bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
            # bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
            # bus["grounded"] .=  0
        # end
    end
    update_load_data!(data_math, load_data, i) ## Update load data

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


## Read load data
laod_data_path = "./data/Load profiles/network_example.csv"
load_data = CSV.read(laod_data_path, DataFrame, header=false)

timesteps = 1:2

## ##################### Conventional inverter #####################
setting = Dict("conventional"=>true, "reconfigurable" => false, "ideal" => false, "dc_link" => true)

i = 1
data_math_NO_conv = build_data_math(data_path, load_data, i; sbase=0.8)
data_math_conv = build_data_math(data_path, load_data, i; setting=setting, sbase=1)
@show data_math_conv["load"]

data_math_conv_sc0_mn = Dict("nw"=>Dict(string(i)=>deepcopy(data_math_NO_conv) for i in 1:length(timesteps)))
data_math_conv_sc1_mn = Dict("nw"=>Dict(string(i)=>deepcopy(data_math_conv) for i in 1:length(timesteps)))
data_math_conv_sc2_mn = Dict("nw"=>Dict(string(i)=>deepcopy(data_math_conv) for i in 1:length(timesteps)))
data_math_conv_sc3_mn = Dict("nw"=>Dict(string(i)=>deepcopy(data_math_conv) for i in 1:length(timesteps)))
data_math_conv_sc4_mn = Dict("nw"=>Dict(string(i)=>deepcopy(data_math_conv) for i in 1:length(timesteps)))

result_conv_sc0_mn = Dict("solution"=>Dict("nw"=>Dict(string(i)=>Dict() for i in 1:length(timesteps))))
result_conv_sc1_mn = Dict("solution"=>Dict("nw"=>Dict(string(i)=>Dict() for i in 1:length(timesteps))))
result_conv_sc2_mn = Dict("solution"=>Dict("nw"=>Dict(string(i)=>Dict() for i in 1:length(timesteps))))
result_conv_sc3_mn = Dict("solution"=>Dict("nw"=>Dict(string(i)=>Dict() for i in 1:length(timesteps))))
result_conv_sc4_mn = Dict("solution"=>Dict("nw"=>Dict(string(i)=>Dict() for i in 1:length(timesteps))))

scenarios_dict = Dict(0 => "No inverter",
                      1 => "Free ripple limit, Free Istc neutral rating",
                      2 => "Free ripple limit, Zero Istc neutral rating", 
                      3 => "Zero ripple limit, Free Istc neutral rating", 
                      4 => "Small ripple limit, Small Istc neutral rating")

for i in collect(timesteps)
    @show i

    # Scenario 0 - No inverter/statcom to compensate unbalance currents
    # sbase = 0.8
    sbase = 1
    data_math_NO_conv, Ibase = build_data_math(data_path, load_data, i; sbase=sbase)
    result_conv_sc0 = run_inverter_case(data_math_NO_conv, setting)
    result_conv_sc0_mn["solution"]["nw"]["$i"] = result_conv_sc0["solution"]
    result_conv_sc0_mn["solution"]["nw"]["$i"]["objective"] = result_conv_sc0["objective"]
    result_conv_sc0_mn["solution"]["nw"]["$i"]["termination_status"] = result_conv_sc0["termination_status"]
    @show "SC0", result_conv_sc0["termination_status"]

    ## Scenario 1 - No ripple or neutral constraints – i.e., both are unconstrained and only phase currents are constrained.
    sbase = 1
    data_math_conv, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
    # data_math_conv["gen"]["1"]["pdcmin"] = 0
    # data_math_conv["gen"]["1"]["pdcmax"] = 100 / sbase
    data_math_conv["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 100] / Ibase
    result_conv_sc1 = run_inverter_case(data_math_conv, setting)
    result_conv_sc1_mn["solution"]["nw"]["$i"] = result_conv_sc1["solution"]
    result_conv_sc1_mn["solution"]["nw"]["$i"]["objective"] = result_conv_sc1["objective"]
    result_conv_sc1_mn["solution"]["nw"]["$i"]["termination_status"] = result_conv_sc1["termination_status"]
    @show "SC1", result_conv_sc1["termination_status"]
    
    ## Scenario 2 - 2w ripple is unconstrained, but the neutral current is set to be fully constrained (ie no neutral current).
    sbase = 1
    data_math_conv, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
    # data_math_conv["gen"]["1"]["pdcmin"] = 0
    # data_math_conv["gen"]["1"]["pdcmax"] = 100 / sbase
    data_math_conv["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 0.0] / Ibase
    result_conv_sc2 = run_inverter_case(data_math_conv, setting)
    result_conv_sc2_mn["solution"]["nw"]["$i"] = result_conv_sc2["solution"]
    result_conv_sc2_mn["solution"]["nw"]["$i"]["objective"] = result_conv_sc2["objective"]
    result_conv_sc2_mn["solution"]["nw"]["$i"]["termination_status"] = result_conv_sc2["termination_status"]
    @show "SC2", result_conv_sc2["termination_status"]

    ## Scenario 3 - The neutral current is unconstrained, but, there is constraint saying there can be no 2w ripple.
    # sbase = 0.8
    sbase = 1
    data_math_conv, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
    data_math_conv["gen"]["1"]["pdcmin"] = 0
    data_math_conv["gen"]["1"]["pdcmax"] = 0
    data_math_conv["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 100] / Ibase
    result_conv_sc3 = run_inverter_case(data_math_conv, setting)
    result_conv_sc3_mn["solution"]["nw"]["$i"] = result_conv_sc3["solution"]
    result_conv_sc3_mn["solution"]["nw"]["$i"]["objective"] = result_conv_sc3["objective"]
    result_conv_sc3_mn["solution"]["nw"]["$i"]["termination_status"] = result_conv_sc3["termination_status"]
    @show "SC3", result_conv_sc3["termination_status"]
    
    ## Scenario 4 - A ripple *and* a neutral constraint. Exactly what value to choose for the constraint values could be based on a few options.
    # sbase = 0.8
    sbase = 1
    data_math_conv, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
    data_math_conv["gen"]["1"]["pdcmin"] = 0
    data_math_conv["gen"]["1"]["pdcmax"] = 5 / sbase
    data_math_conv["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 30] / Ibase
    result_conv_sc4 = run_inverter_case(data_math_conv, setting)
    result_conv_sc4_mn["solution"]["nw"]["$i"] = result_conv_sc4["solution"]
    result_conv_sc4_mn["solution"]["nw"]["$i"]["objective"] = result_conv_sc4["objective"]
    result_conv_sc4_mn["solution"]["nw"]["$i"]["termination_status"] = result_conv_sc4["termination_status"]
    @show "SC4", result_conv_sc4["termination_status"]
    if result_conv_sc4["termination_status"] ∉ [LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
        sbase = 1.2
        data_math_conv, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
        data_math_conv["gen"]["1"]["pdcmin"] = 0
        data_math_conv["gen"]["1"]["pdcmax"] = 5 / sbase
        data_math_conv["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 30] / Ibase
        result_conv_sc4 = run_inverter_case(data_math_conv, setting)
        result_conv_sc4_mn["solution"]["nw"]["$i"] = result_conv_sc4["solution"]
        result_conv_sc4_mn["solution"]["nw"]["$i"]["objective"] = result_conv_sc4["objective"]
        result_conv_sc4_mn["solution"]["nw"]["$i"]["termination_status"] = result_conv_sc4["termination_status"]
        @show "SC4", result_conv_sc4["termination_status"]
    end
end


##
function get_currents(result, comp, id, idx, timesteps; vnom=0.23094)
    # vm_nom = data_math_conv["load"]["1"]["vbase"]
    if comp == "load"
        key_re = "crd"
        key_im = "cid"
    elseif comp == "branch"
        key_re = "cr_fr"
        key_im = "ci_fr"
    elseif comp == "gen"
        key_re = "crg"
        key_im = "cig"
    end

    c_re = [result["solution"]["nw"]["$i"][comp][id][key_re][idx]/vm_nom for i in timesteps]
    c_im = [result["solution"]["nw"]["$i"][comp][id][key_im][idx]/vm_nom for i in timesteps]
    c = c_re .+ im*c_im

    return c_re, c_im, c
end


function get_results(result, timesteps; no_inverter=false)
    load_1_c_re, load_1_c_im, load_1_c = get_currents(result, "load", "1", 1, timesteps)
    load_2_c_re, load_2_c_im, load_2_c = get_currents(result, "load", "2", 1, timesteps)
    load_3_c_re, load_3_c_im, load_3_c = get_currents(result, "load", "3", 1, timesteps)
    load_c_re = [load_3_c_re  load_2_c_re  load_1_c_re]
    load_c_im = [load_3_c_im  load_2_c_im  load_1_c_im]
    load_c = [load_3_c  load_2_c  load_1_c]

    branch_cr_1, branch_ci_1, branch_c_1 = get_currents(result, "branch", "1", 1, timesteps)
    branch_cr_2, branch_ci_2, branch_c_2 = get_currents(result, "branch", "1", 2, timesteps)
    branch_cr_3, branch_ci_3, branch_c_3 = get_currents(result, "branch", "1", 3, timesteps)
    branch_cr_4, branch_ci_4, branch_c_4 = get_currents(result, "branch", "1", 4, timesteps)
    branch_cr = [branch_cr_1  branch_cr_2  branch_cr_3  branch_cr_4]
    branch_ci = [branch_ci_1  branch_ci_2  branch_ci_3  branch_ci_4]
    branch_c = [branch_c_1  branch_c_2  branch_c_3  branch_c_4]

    vsource_crg_1, vsource_cig_1, vsource_cg_1 = get_currents(result, "gen", "2", 1, timesteps)
    vsource_crg_2, vsource_cig_2, vsource_cg_2 = get_currents(result, "gen", "2", 2, timesteps)
    vsource_crg_3, vsource_cig_3, vsource_cg_3 = get_currents(result, "gen", "2", 3, timesteps)
    vsource_crg_4 = -(vsource_crg_1 + vsource_crg_2 + vsource_crg_3)
    vsource_cig_4 = -(vsource_cig_1 + vsource_cig_2 + vsource_cig_3)
    vsource_cg_4 = -(vsource_cg_1 + vsource_cg_2 + vsource_cg_3)
    vsource_crg = [vsource_crg_1 vsource_crg_2  vsource_crg_3  vsource_crg_4]
    vsource_cig = [vsource_cig_1  vsource_cig_2  vsource_cig_3  vsource_cig_4]
    vsource_cg = [vsource_cg_1  vsource_cg_2  vsource_cg_3  vsource_cg_4]

    inverter_cg = []
    p_dclink = []
    if !no_inverter
        inverter_cr_1, inverter_ci_1, inverter_c_1 = get_currents(result, "branch", "2", 1, timesteps)
        inverter_cr_2, inverter_ci_2, inverter_c_2 = get_currents(result, "branch", "2", 2, timesteps)
        inverter_cr_3, inverter_ci_3, inverter_c_3 = get_currents(result, "branch", "2", 3, timesteps)
        inverter_cr_4, inverter_ci_4, inverter_c_4 = get_currents(result, "branch", "2", 4, timesteps)
        inverter_cr = [inverter_cr_1  inverter_cr_2  inverter_cr_3  inverter_cr_4]
        inverter_ci = [inverter_ci_1  inverter_ci_2  inverter_ci_3  inverter_ci_4]
        inverter_c = [inverter_c_1  inverter_c_2  inverter_c_3  inverter_c_4]

        inverter_crg_1, inverter_cig_1, inverter_cg_1 = get_currents(result, "gen", "1", 1, timesteps)
        inverter_crg_2, inverter_cig_2, inverter_cg_2 = get_currents(result, "gen", "1", 2, timesteps)
        inverter_crg_3, inverter_cig_3, inverter_cg_3 = get_currents(result, "gen", "1", 3, timesteps)
        inverter_crg_4 = -(inverter_crg_1 + inverter_crg_2 + inverter_crg_3)
        inverter_cig_4 = -(inverter_cig_1 + inverter_cig_2 + inverter_cig_3)
        inverter_cg_4 = -(inverter_cg_1 + inverter_cg_2 + inverter_cg_3)
        inverter_crg = [inverter_crg_1  inverter_crg_2  inverter_crg_3  inverter_crg_4]
        inverter_cig = [inverter_cig_1  inverter_cig_2  inverter_cig_3  inverter_cig_4]
        inverter_cg = [inverter_cg_1  inverter_cg_2  inverter_cg_3  inverter_cg_4]

        p_dclink = [result["solution"]["nw"]["$i"]["gen"]["1"]["pdc_link"][1] for i in timesteps]
    end

    return load_c, vsource_cg, inverter_cg, p_dclink
    
end

result_sc = deepcopy(result_conv_sc2_mn)
vm_nom = data_math_conv["load"]["1"]["vbase"]
load_c, vsource_cg, inverter_cg, p_dclink = get_results(result_sc, timesteps)

currents = plot(minimum(abs.(load_c), dims=2), fillrange = maximum(abs.(load_c), dims=2), fillalpha = 0.2, c = :grey, label = false, legend = :best)
plot!(abs.(load_c[:,1]), color=1, linewidth=1, label="phase a")
plot!(abs.(load_c[:,2]), color=2, linewidth=1, label="phase b")
plot!(abs.(load_c[:,3]), color=3, linewidth=1, label="phase c")
# plot!((load_1+load_2+load_3)/3, color=:black, linewidth=1, label="mean abc")
# savefig(currents, "Figures/STATCOM_load.pdf")

##

function plot_results(result, k, timesteps; Iseq_plt=plot(), ldc_plt=plot(), no_inverter=false)
    load_c, vsource_cg, inverter_cg, p_dclink = get_results(result, timesteps, no_inverter=no_inverter)

    ######################## plot magnitudes
    source_plt = plot(minimum(abs.(load_c), dims=2), fillrange = maximum(abs.(load_c), dims=2), fillalpha = 0.2, c = :grey, label = false, legend = :topleft)
    # plot!(branch_cm_1, color=1, linewidth=1, linestyle=:dash, label="phase a")
    # plot!(branch_cm_2, color=2, linewidth=1, linestyle=:dash, label="phase b")
    # plot!(branch_cm_3, color=3, linewidth=1, linestyle=:dash, label="phase c")
    # plot!(branch_cm_4, color=4, linewidth=1, linestyle=:dash, label="neutral")
    # plot!((branch_cm_1+branch_cm_2+branch_cm_3)/3, color=:black, linewidth=1, label="mean abc")
    plot!(abs.(vsource_cg[:,1]), color=1, linewidth=1, linestyle=:dash, label="phase a")
    plot!(abs.(vsource_cg[:,2]), color=2, linewidth=1, linestyle=:dash, label="phase b")
    plot!(abs.(vsource_cg[:,3]), color=3, linewidth=1, linestyle=:dash, label="phase c")
    plot!(abs.(vsource_cg[:,4]), color=4, linewidth=1, linestyle=:dash, label="neutral")
    plot!((abs.(vsource_cg[:,1])+abs.(vsource_cg[:,2])+abs.(vsource_cg[:,3]))/3, color=:black, linewidth=1, label="mean abc")
    ylabel!("Im src")
    title!(scenarios_dict[k])

    if !no_inverter
        # statcom_plt = plot(inverter_cm_1, color=1, linewidth=1, linestyle=:dash, label=false)
        # plot!(inverter_cm_2, color=2, linewidth=1, linestyle=:dash, label=false)
        # plot!(inverter_cm_3, color=3, linewidth=1, linestyle=:dash, label=false)
        # plot!(inverter_cm_4, color=4, linewidth=1, linestyle=:dash, label=false)
        # plot!((inverter_cm_1+inverter_cm_2+inverter_cm_3)/3, color=:black, linewidth=1, label=false, size=(600, 200))
        statcom_plt = plot(abs.(inverter_cg[:,1]), color=1, linewidth=1, linestyle=:dash, label=false)
        plot!(abs.(inverter_cg[:,2]), color=2, linewidth=1, linestyle=:dash, label=false)
        plot!(abs.(inverter_cg[:,3]), color=3, linewidth=1, linestyle=:dash, label=false)
        plot!(abs.(inverter_cg[:,4]), color=4, linewidth=1, linestyle=:dash, label=false)
        plot!((abs.(inverter_cg[:,1])+abs.(inverter_cg[:,2])+abs.(inverter_cg[:,3]))/3, color=:black, linewidth=1, label=false, size=(600, 200))
        ylabel!("Im stc")

        p_dclink_plt = plot(p_dclink, label=false, ylabel="Pcap")
    else
        statcom_plt = plot()
        p_dclink_plt = plot()
    end

    l = @layout [a{0.5h} ; b{0.3h} ; c{0.2h}]
    currents_sc = plot(source_plt, statcom_plt, p_dclink_plt, layout = l, legend=false)
    savefig(currents_sc, "Figures/STATCOM_load_sc$k.pdf")


    # #################### plot real and imaginary
    # source_re_plt = plot(minimum([load_1_c_re load_2_c_re load_3_c_re], dims=2), fillrange = maximum([load_1_c_re load_2_c_re load_3_c_re], dims=2), fillalpha = 0.2, c = :grey, label =false, legend = :topleft)
    # plot!(branch_cr_1, color=1, linewidth=1, linestyle=:solid, label="phase a")
    # plot!(branch_cr_2, color=2, linewidth=1, linestyle=:solid, label="phase b")
    # plot!(branch_cr_3, color=3, linewidth=1, linestyle=:solid, label="phase c")
    # plot!(branch_cr_4, color=4, linewidth=1, linestyle=:solid, label="neutral")
    # # plot!((branch_cr_1+branch_cr_2+branch_cr_3)/3, color=:black, linewidth=1, label="mean abc")
    # ylabel!("Ire src")
    # title!(scenarios_dict[k])

    # source_im_plt = plot(minimum([load_1_c_im load_2_c_im load_3_c_im], dims=2), fillrange = maximum([load_1_c_im load_2_c_im load_3_c_im], dims=2), fillalpha = 0.2, c = :grey, label =false, legend = :topleft)
    # plot!(branch_ci_1, color=1, linewidth=1, linestyle=:solid, label="phase a")
    # plot!(branch_ci_2, color=2, linewidth=1, linestyle=:solid, label="phase b")
    # plot!(branch_ci_3, color=3, linewidth=1, linestyle=:solid, label="phase c")
    # plot!(branch_ci_4, color=4, linewidth=1, linestyle=:solid, label="neutral")
    # # plot!((branch_ci_1+branch_ci_2+branch_ci_3)/3, color=:black, linewidth=1, label="mean abc")
    # ylabel!("Iim src")

    # if !no_inverter
    #     statcom_re_plt = plot(inverter_cr_1, color=1, linewidth=1, linestyle=:solid, label=false)
    #     plot!(inverter_cr_2, color=2, linewidth=1, linestyle=:solid, label=false)
    #     plot!(inverter_cr_3, color=3, linewidth=1, linestyle=:solid, label=false)
    #     plot!(inverter_cr_4, color=4, linewidth=1, linestyle=:solid, label=false)
    #     # plot!((inverter_cr_1+inverter_cr_2+inverter_cr_3)/3, color=:black, linewidth=1, label=false, size=(600, 200))
    #     ylabel!("Ire stc")

    #     statcom_im_plt = plot(inverter_ci_1, color=1, linewidth=1, linestyle=:solid, label=false)
    #     plot!(inverter_ci_2, color=2, linewidth=1, linestyle=:solid, label=false)
    #     plot!(inverter_ci_3, color=3, linewidth=1, linestyle=:solid, label=false)
    #     plot!(inverter_ci_4, color=4, linewidth=1, linestyle=:solid, label=false)
    #     # plot!((inverter_ci_1+inverter_ci_2+inverter_ci_3)/3, color=:black, linewidth=1, label=false, size=(600, 200))
    #     ylabel!("Iim stc")
    # else
    #     statcom_re_plt = plot()
    #     statcom_im_plt = plot()
    # end


    # l = @layout [a b; c d]
    # currents_sc = plot(source_re_plt, statcom_re_plt, source_im_plt, statcom_im_plt, layout = l, legend=false)
    # savefig(currents_sc, "Figures/STATCOM_re_im_load_sc$k.pdf")
    
    #################### plot negative current sequence
    I_seq_m_load_mn = Vector{Float64}(undef, 3)
    I_seq_m_source_mn = Vector{Float64}(undef, 3)
    for i in timesteps
        I_seq_re, I_seq_im, I_seq_m_load = RPMD.get_sequence_components(load_c[i,:])
        I_seq_re, I_seq_im, I_seq_m_source = RPMD.get_sequence_components(vsource_cg[i,1:3])
        I_seq_m_load_mn = [I_seq_m_load_mn I_seq_m_load]
        I_seq_m_source_mn = [I_seq_m_source_mn I_seq_m_source]
    end
    # if k==1   # only plot the load current negative sequence ONCE
    #     plot!(Iseq_plt, I_seq_m_load_mn[3,2:end], color=1, label="Load")
    # end
    plot!(Iseq_plt, I_seq_m_source_mn[3,2:end], label="Sc. $k")
    savefig(Iseq_plt, "Figures/STATCOM_load_Iseq.pdf")


    plot_load_duration_curve!(I_seq_m_source_mn; plt=Iseq_ldc_plt, label="Sc $k", title="I neg seq duration curve", xlabel="Percentage of Time (%)", ylabel="I neg seq (A)")
    savefig(Iseq_ldc_plt, "Figures/STATCOM_Iseq_LDC.pdf")

    # vr_1 = [result["solution"]["nw"]["$i"]["bus"]["1"]["vr"] for i in timesteps] # load bus
    # vi_1 = [result["solution"]["nw"]["$i"]["bus"]["1"]["vi"][1] for i in timesteps] # load bus
    # vm_1 = abs.(vr_1 + im*vi_1)
    # vm_1 = [abs.(result["solution"]["nw"]["$i"]["bus"]["1"]["vr"].+im*result["solution"]["nw"]["$i"]["bus"]["1"]["vi"]) for i in timesteps] # load bus

    # vr_2 = [result["solution"]["nw"]["$i"]["bus"]["2"]["vr"][1] for i in timesteps] # source bus
    # vi_2 = [result["solution"]["nw"]["$i"]["bus"]["2"]["vi"][1] for i in timesteps] # source bus
    # vm_2 = abs.(vr_2 + im*vi_2)
    # vm_2 = [abs.(result["solution"]["nw"]["$i"]["bus"]["2"]["vr"].+im*result["solution"]["nw"]["$i"]["bus"]["2"]["vi"]) for i in timesteps] # source bus

    # vr_3 = [result["solution"]["nw"]["$i"]["bus"]["3"]["vr"][1] for i in timesteps] # inverter bus
    # vi_3 = [result["solution"]["nw"]["$i"]["bus"]["3"]["vi"][1] for i in timesteps] # inverter bus
    # vm_3 = abs.(vr_3 + im*vi_3)

    # vm_3 = [abs.(result["solution"]["nw"]["$i"]["bus"]["3"]["vr"].+im*result["solution"]["nw"]["$i"]["bus"]["3"]["vi"]) for i in timesteps] # inverter bus

    # plot(vm_1, vm_2, vm_3)

    return currents_sc
end

# ldc_plt = plot()
# data_no_comp = [maximum([k for k in load_data[i,:]]) for i in timesteps]
# plot_load_duration_curve!(data_no_comp; plt=ldc_plt, label="No compensation")

Iseq_ldc_plt = plot()
Iseq_plt = plot()
currents_sc_plt = Dict()
for (k, result) in enumerate([result_conv_sc0_mn result_conv_sc1_mn result_conv_sc2_mn result_conv_sc3_mn result_conv_sc4_mn])
# for (k, result) in enumerate([result_conv_sc2_mn])
    if k==1 # no inverter => no compensation
        currents_sc_plt[k] = plot_results(result, k-1, timesteps; Iseq_plt=Iseq_plt, ldc_plt=Iseq_ldc_plt, no_inverter=true)
    else
        currents_sc_plt[k] = plot_results(result, k-1, timesteps; Iseq_plt=Iseq_plt, ldc_plt=Iseq_ldc_plt)
    end
end


##
ldc_plt = plot()
for (k, result) in enumerate([result_conv_sc0_mn result_conv_sc1_mn result_conv_sc2_mn result_conv_sc3_mn result_conv_sc4_mn])
# for (k, result) in enumerate([result_conv_sc2_mn])
    markershape = :none
    if k == 1
        markershape = :circle
    end

    branch_cr_1, branch_ci_1, branch_c_1 = get_currents(result, "branch", "1", 1, timesteps)
    branch_cr_2, branch_ci_2, branch_c_2 = get_currents(result, "branch", "1", 2, timesteps)
    branch_cr_3, branch_ci_3, branch_c_3 = get_currents(result, "branch", "1", 3, timesteps)
    branch_cr_4, branch_ci_4, branch_c_4 = get_currents(result, "branch", "1", 4, timesteps)
    
    #################### plot load duration curves
    data = maximum(abs.([branch_c_1 branch_c_2 branch_c_3 branch_c_4]), dims=2)
    plot_load_duration_curve!(data; plt=ldc_plt, label="Sc $(k-1)", markershape=markershape, title="Max source phase current", xlabel="Percentage of Time (%)", ylabel="max phase I src (A)")
end
display(ldc_plt)
savefig(ldc_plt, "Figures/STATCOM_Isrc_max_LDC.pdf")