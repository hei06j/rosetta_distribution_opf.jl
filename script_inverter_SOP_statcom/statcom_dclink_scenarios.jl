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
using JLD2

const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

PMD.silence!()


# ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes", "warm_start_init_point"=>"yes", "max_iter"=>100000)
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes", "warm_start_init_point"=>"yes", "max_iter"=>100000)
# ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer)
# highs_solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
# juniper_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt_solver, "mip_solver" => highs_solver)

data_path = "./data/small networks/inverter_4w_wye_unbalanced_loads_2bus.dss"
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
    plot!(plt, collect(1:n), sorted_data, 
        label=label, 
        xlabel=xlabel, 
        ylabel=ylabel, 
        # title=title, 
        # lw=0.5,
        markersize=1, 
        grid=true,
        markershape=markershape)
end

function build_data_math(data_path, load_data, i; sbase=0.8, setting=nothing)
    data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
    data_eng["settings"]["sbase_default"] = sbase
    data_eng["voltage_source"]["source"]["rs"] *= 0
    data_eng["voltage_source"]["source"]["xs"] *= 0

    if isnothing(setting)
        # @show setting
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

function get_currents(result, comp, id, idx, timesteps; vm_nom=0.23094)
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

        p_dclink = [sqrt(result["solution"]["nw"]["$i"]["gen"]["1"]["pdc_link_sqr"][1]) for i in timesteps]
    end

    return load_c, vsource_cg, inverter_cg, p_dclink
    
end

function plot_results(result, k, timesteps; Iseq_plt=plot(), ldc_plt=plot(), no_inverter=false)
    load_c, vsource_cg, inverter_cg, p_dclink = get_results(result, timesteps, no_inverter=no_inverter)
    
    ######################## plot magnitudes
    source_plt = plot(minimum(abs.(load_c), dims=2), fillrange = maximum(abs.(load_c), dims=2), fillalpha = 0.4, c = :grey, label = "Before STATCOM", legend = :topright, size=(600, 400))
    plot!((abs.(load_c[:,1])+abs.(load_c[:,2])+abs.(load_c[:,3]))/3, color=:black, linewidth=1, label="mean before")
    plot!(minimum(abs.(vsource_cg[:,1:3]), dims=2), fillrange = maximum(abs.(vsource_cg[:,1:3]), dims=2), fillalpha = 0.4, c = :red, label = "After STATCOM", legend = :topright, size=(600, 400))
    plot!((abs.(vsource_cg[:,1])+abs.(vsource_cg[:,2])+abs.(vsource_cg[:,3]))/3, color=:red, linewidth=1, label="mean after", xticks=(collect(1:12*6:289), ["00:00", "12:00", "00:00", "12:00", "00:00"]))
    # plot!(abs.(vsource_cg[:,1]), color=1, linewidth=1, linestyle=:solid, label="phase a")
    # plot!(abs.(vsource_cg[:,2]), color=2, linewidth=1, linestyle=:solid, label="phase b")
    # plot!(abs.(vsource_cg[:,3]), color=3, linewidth=1, linestyle=:solid, label="phase c")
    # plot!(abs.(vsource_cg[:,4]), color=4, linewidth=1, linestyle=:solid, label="neutral")
    ylabel!("Source Current (A)")
    xlabel!("Time (h)")
    # ylabel!(L"|I_{source}|")
    # title!(scenarios_dict[k])

    if !no_inverter
        statcom_plt = plot(abs.(inverter_cg[:,1]), color=1, linewidth=1, linestyle=:dash, label=false)
        plot!(abs.(inverter_cg[:,2]), color=2, linewidth=1, linestyle=:dash, label=false)
        plot!(abs.(inverter_cg[:,3]), color=3, linewidth=1, linestyle=:dash, label=false)
        plot!(abs.(inverter_cg[:,4]), color=4, linewidth=1, linestyle=:dash, label=false)
        plot!((abs.(inverter_cg[:,1])+abs.(inverter_cg[:,2])+abs.(inverter_cg[:,3]))/3, color=:black, linewidth=1, label=false, size=(600, 200))
        ylabel!(L"|I_{statcom}|")

        p_dclink_plt = plot(p_dclink, label=false, ylabel=L"P_{\mathrm{dc\,link}}^{2w}")
    else
        statcom_plt = plot()
        p_dclink_plt = plot()
    end

    l = @layout [a{0.5h} ; b{0.3h} ; c{0.2h}]
    currents_sc = plot(source_plt, statcom_plt, p_dclink_plt, layout = l, legend=false)
    Plots.savefig(currents_sc, "Figures/STATCOM_load_sc$k.pdf")
    Plots.savefig(source_plt, "Figures/STATCOM_Isrc_sc$k.pdf")
    
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
    Plots.savefig(Iseq_plt, "Figures/STATCOM_load_Iseq.pdf")


    # plot_load_duration_curve!(I_seq_m_source_mn; plt=Iseq_ldc_plt, label="Sc $k", xlabel="Percentage of Time (%)", ylabel=L"I^{-} \mathrm{(A)}")
    # Plots.savefig(Iseq_ldc_plt, "Figures/STATCOM_Iseq_LDC.pdf")

    return currents_sc
end


### Read load data
load_data_path = "./data/Load profiles/network_example.csv"
load_data = CSV.read(load_data_path, DataFrame, header=false)

timesteps = 1:8806

timesteps = 1:1

## ##################### Conventional inverter #####################
setting = Dict("conventional"=>true, "reconfigurable" => false, "ideal" => false, "dc_link" => true)

i = 1
data_math_NO_conv = build_data_math(data_path, load_data, i; sbase=0.8)
data_math_conv = build_data_math(data_path, load_data, i; setting=setting, sbase=1)
# @show data_math_conv[1]["load"]

# data_math_conv_sc0_mn = Dict("nw"=>Dict(string(i)=>deepcopy(data_math_NO_conv) for i in 1:length(timesteps)))
# data_math_conv_sc1_mn = Dict("nw"=>Dict(string(i)=>deepcopy(data_math_conv) for i in 1:length(timesteps)))
# data_math_conv_sc2_mn = Dict("nw"=>Dict(string(i)=>deepcopy(data_math_conv) for i in 1:length(timesteps)))
# data_math_conv_sc3_mn = Dict("nw"=>Dict(string(i)=>deepcopy(data_math_conv) for i in 1:length(timesteps)))
# data_math_conv_sc4_mn = Dict("nw"=>Dict(string(i)=>deepcopy(data_math_conv) for i in 1:length(timesteps)))

result_conv_sc0_mn = Dict("solution"=>Dict("nw"=>Dict(string(i)=>Dict() for i in 1:length(timesteps))))
result_conv_sc1_mn = Dict("solution"=>Dict("nw"=>Dict(string(i)=>Dict() for i in 1:length(timesteps))))
result_conv_sc2_mn = Dict("solution"=>Dict("nw"=>Dict(string(i)=>Dict() for i in 1:length(timesteps))))
result_conv_sc3_mn = Dict("solution"=>Dict("nw"=>Dict(string(i)=>Dict() for i in 1:length(timesteps))))
result_conv_sc4_mn = Dict("solution"=>Dict("nw"=>Dict(string(i)=>Dict() for i in 1:length(timesteps))))

scenarios_dict = Dict(0 => "No STATCOM",
                      1 => "Free 2w ripple limit, Free STATCOM neutral rating",
                      2 => "Free 2w ripple limit, Zero STATCOM neutral rating", 
                      3 => "Zero 2w ripple limit, Free STATCOM neutral rating", 
                      4 => "Small 2w ripple limit, Small STATCOM neutral rating")
data_math_conv = nothing
for i in collect(timesteps)
    @show i

    # # Scenario 0 - No inverter/statcom to compensate unbalance currents
    # # sbase = 0.8
    # sbase = 1
    # data_math_NO_conv, Ibase = build_data_math(data_path, load_data, i; sbase=sbase)
    # result_conv_sc0 = run_inverter_case(data_math_NO_conv, setting)
    # result_conv_sc0_mn["solution"]["nw"]["$i"] = result_conv_sc0["solution"]
    # result_conv_sc0_mn["solution"]["nw"]["$i"]["objective"] = result_conv_sc0["objective"]
    # result_conv_sc0_mn["solution"]["nw"]["$i"]["termination_status"] = result_conv_sc0["termination_status"]
    # # @show "SC0", result_conv_sc0["termination_status"]

    # ## Scenario 1 - No ripple or neutral constraints – i.e., both are unconstrained and only phase currents are constrained.
    # sbase = 1
    # data_math_conv, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
    # # data_math_conv["gen"]["1"]["pdcmin"] = 0
    # # data_math_conv["gen"]["1"]["pdcmax"] = 100 / sbase
    # data_math_conv["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 100] / Ibase
    # result_conv_sc1 = run_inverter_case(data_math_conv, setting)
    # result_conv_sc1_mn["solution"]["nw"]["$i"] = result_conv_sc1["solution"]
    # result_conv_sc1_mn["solution"]["nw"]["$i"]["objective"] = result_conv_sc1["objective"]
    # result_conv_sc1_mn["solution"]["nw"]["$i"]["termination_status"] = result_conv_sc1["termination_status"]
    # # @show "SC1", result_conv_sc1["termination_status"]
    
    # ## Scenario 2 - 2w ripple is unconstrained, but the neutral current is set to be fully constrained (ie no neutral current).
    # sbase = 1
    # data_math_conv, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
    # # data_math_conv["gen"]["1"]["pdcmin"] = 0
    # # data_math_conv["gen"]["1"]["pdcmax"] = 100 / sbase
    # data_math_conv["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 0.0] / Ibase
    # result_conv_sc2 = run_inverter_case(data_math_conv, setting)
    # result_conv_sc2_mn["solution"]["nw"]["$i"] = result_conv_sc2["solution"]
    # result_conv_sc2_mn["solution"]["nw"]["$i"]["objective"] = result_conv_sc2["objective"]
    # result_conv_sc2_mn["solution"]["nw"]["$i"]["termination_status"] = result_conv_sc2["termination_status"]
    # # @show "SC2", result_conv_sc2["termination_status"]
    
    ## Scenario 3 - The neutral current is unconstrained, but, there is constraint saying there can be no 2w ripple.
    # sbase = 0.8
    sbase = 1
    data_math_conv, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
    data_math_conv["gen"]["1"]["pdcmin"] = 0
    data_math_conv["gen"]["1"]["pdcmax"] = 0
    data_math_conv["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 100] / Ibase
    data_math_conv["bus"]["3"]["vmax"] = [1.1, 1,1, 1.1, 0]
    data_math_conv["bus"]["3"]["vmin"] = [0.9, 0.9, 0.9, 0]
    result_conv_sc3 = run_inverter_case(data_math_conv, setting)
    result_conv_sc3_mn["solution"]["nw"]["$i"] = result_conv_sc3["solution"]
    result_conv_sc3_mn["solution"]["nw"]["$i"]["objective"] = result_conv_sc3["objective"]
    result_conv_sc3_mn["solution"]["nw"]["$i"]["termination_status"] = result_conv_sc3["termination_status"]
    # @show "SC3", result_conv_sc3["termination_status"]


    sbase = data_math_conv["settings"]["sbase"]                          # p.u.
    sbace_factor = data_math_conv["settings"]["power_scale_factor"]      # 
    vbase_factor = data_math_conv["settings"]["voltage_scale_factor"]
    vbase = 0.2309      # [kV]  data_math["settings"]["vbases_default"]["5"]
    Ibase = (sbase * sbace_factor) / (vbase * vbase_factor)  #[kA]
    
    cd = [result_conv_sc3["solution"]["load"]["3"]["crd"] + im*result_conv_sc3["solution"]["load"]["3"]["cid"]
          result_conv_sc3["solution"]["load"]["2"]["crd"] + im*result_conv_sc3["solution"]["load"]["2"]["cid"]
          result_conv_sc3["solution"]["load"]["1"]["crd"] + im*result_conv_sc3["solution"]["load"]["1"]["cid"]]
    _, _, seqd = RPMD.get_sequence_components(cd)
    # seqd = RPMD.sequence(cd)
    # cd = [cd; -sum(cd)]
    # @show abs.(cd) * Ibase
    @show seqd

    cg1 = result_conv_sc3["solution"]["gen"]["1"]["crg"] + im*result_conv_sc3["solution"]["gen"]["1"]["cig"]
    _, _, seq_g1 = RPMD.get_sequence_components(cg1)
    # seq_g1 = RPMD.sequence(cg1)

    cg2 = result_conv_sc3["solution"]["gen"]["2"]["crg"] + im*result_conv_sc3["solution"]["gen"]["2"]["cig"]
    _, _, seq_g2 = RPMD.get_sequence_components(cg2)
    # seq_g2 = RPMD.sequence(cg2)

    # cd = [cd; -sum(cd)]
    # @show abs.(cd) * Ibase
    @show seq_g1, seq_g2
    @show sum(cg1), sum(cg2)


    c1 = result_conv_sc3["solution"]["branch"]["1"]["cr_fr"] + im*result_conv_sc3["solution"]["branch"]["1"]["ci_fr"]
    _, _, seq1 = RPMD.get_sequence_components(c1[1:3])
    # seq1 = RPMD.sequence(c1[1:3])

    c2 = result_conv_sc3["solution"]["branch"]["2"]["cr_fr"] + im*result_conv_sc3["solution"]["branch"]["2"]["ci_fr"]
    _, _, seq2 = RPMD.get_sequence_components(c2[1:3])
    # seq2 = RPMD.sequence(c2[1:3])
    # @show abs.(c1), abs.(c2)
    @show seq2, seq1

    v1 = abs.(result_conv_sc3["solution"]["bus"]["1"]["vr"] + im*result_conv_sc3["solution"]["bus"]["1"]["vi"])
    v2 = abs.(result_conv_sc3["solution"]["bus"]["2"]["vr"] + im*result_conv_sc3["solution"]["bus"]["2"]["vi"])
    v3 = abs.(result_conv_sc3["solution"]["bus"]["3"]["vr"] + im*result_conv_sc3["solution"]["bus"]["3"]["vi"])
    _, _, v3 = RPMD.get_sequence_components(result_conv_sc3["solution"]["bus"]["3"]["vr"][1:3] + im*result_conv_sc3["solution"]["bus"]["3"]["vi"][1:3])
    @show v3

    # ## Scenario 4 - A ripple *and* a neutral constraint. Exactly what value to choose for the constraint values could be based on a few options.
    # # sbase = 0.8
    # sbase = 1
    # data_math_conv, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
    # data_math_conv["gen"]["1"]["pdcmin"] = 0
    # data_math_conv["gen"]["1"]["pdcmax"] = 5 / sbase
    # data_math_conv["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 30] / Ibase
    # result_conv_sc4 = run_inverter_case(data_math_conv, setting)
    # result_conv_sc4_mn["solution"]["nw"]["$i"] = result_conv_sc4["solution"]
    # result_conv_sc4_mn["solution"]["nw"]["$i"]["objective"] = result_conv_sc4["objective"]
    # result_conv_sc4_mn["solution"]["nw"]["$i"]["termination_status"] = result_conv_sc4["termination_status"]
    # # @show "SC4", result_conv_sc4["termination_status"]
    # if result_conv_sc4["termination_status"] ∉ [LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
    #     sbase = 1.2
    #     data_math_conv, Ibase = build_data_math(data_path, load_data, i; setting=setting, sbase=sbase)
    #     data_math_conv["gen"]["1"]["pdcmin"] = 0
    #     data_math_conv["gen"]["1"]["pdcmax"] = 5 / sbase
    #     data_math_conv["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 30] / Ibase
    #     result_conv_sc4 = run_inverter_case(data_math_conv, setting)
    #     result_conv_sc4_mn["solution"]["nw"]["$i"] = result_conv_sc4["solution"]
    #     result_conv_sc4_mn["solution"]["nw"]["$i"]["objective"] = result_conv_sc4["objective"]
    #     result_conv_sc4_mn["solution"]["nw"]["$i"]["termination_status"] = result_conv_sc4["termination_status"]
    #     # @show "SC4", result_conv_sc4["termination_status"]
    # end
end

# @save "result_conv_sc0_mn.jld2" result_conv_sc0_mn
# @save "result_conv_sc1_mn.jld2" result_conv_sc1_mn
# @save "result_conv_sc2_mn.jld2" result_conv_sc2_mn
# @save "result_conv_sc3_mn.jld2" result_conv_sc3_mn
# @save "result_conv_sc4_mn.jld2" result_conv_sc4_mn


# @load "result_conv_sc0_mn.jld2" result_conv_sc0_mn
# @load "result_conv_sc1_mn.jld2" result_conv_sc1_mn
# @load "result_conv_sc2_mn.jld2" result_conv_sc2_mn
# @load "result_conv_sc3_mn.jld2" result_conv_sc3_mn
# @load "result_conv_sc4_mn.jld2" result_conv_sc4_mn

##
timesteps_plt = 1:48*6


result_sc = deepcopy(result_conv_sc2_mn)
# vm_nom = data_math_conv["load"]["1"]["vbase"]
load_c, vsource_cg, inverter_cg, p_dclink = get_results(result_sc, timesteps_plt)

currents = plot(minimum(abs.(load_c), dims=2), fillrange = maximum(abs.(load_c), dims=2), fillalpha = 0.2, c = :grey, label = false, legend = :best)
plot!(abs.(load_c[:,1]), color=1, linewidth=1, label="phase a")
plot!(abs.(load_c[:,2]), color=2, linewidth=1, label="phase b")
plot!(abs.(load_c[:,3]), color=3, linewidth=1, label="phase c")
# plot!((load_1+load_2+load_3)/3, color=:black, linewidth=1, label="mean abc")
# savefig(currents, "Figures/STATCOM_load.pdf")

##
plotlyjs()


# ldc_plt = plot()
# data_no_comp = [maximum([k for k in load_data[i,:]]) for i in timesteps]
# plot_load_duration_curve!(data_no_comp; plt=ldc_plt, label="No compensation")

Iseq_ldc_plt = Plots.plot()
Iseq_plt = Plots.plot()
currents_sc_plt = Dict()
for (k, result) in enumerate([result_conv_sc0_mn result_conv_sc1_mn result_conv_sc2_mn result_conv_sc3_mn result_conv_sc4_mn])
# for (k, result) in enumerate([result_conv_sc2_mn])
    if k==1 # no inverter => no compensation
        currents_sc_plt[k] = plot_results(result, k-1, timesteps_plt; Iseq_plt=Iseq_plt, ldc_plt=Iseq_ldc_plt, no_inverter=true)
    else
        currents_sc_plt[k] = plot_results(result, k-1, timesteps_plt; Iseq_plt=Iseq_plt, ldc_plt=Iseq_ldc_plt)
    end
end


##
all_labels = ["Case 1a) uncst. 2w & neutral", "Case 1b) uncst. 2w & no neutral", "Case 1c) no 2w & uncst. neutral", "Case 1d) cst. 2w & neutral", "Non-Mitigated F2 Substation Current"]
tickvals = [1, 72, 144, 216, 288]
ticktext = ["00:00", "12:00", "00:00", "12:00", "00:00"]

function plot_results_plotlyjs(result, k, timesteps; Iseq_plt=nothing, ldc_plt=nothing, no_inverter=false)
    load_c, vsource_cg, inverter_cg, p_dclink = get_results(result, timesteps, no_inverter=no_inverter)

    n = size(load_c, 1)
    time = 1:n

    # Font and style settings
    latex_font = PlotlyJS.attr(family="serif", size=22, color="black")
    axis_font = PlotlyJS.attr(size=22, family="serif")
    tick_font = PlotlyJS.attr(size=18, family="serif")
    legend_style = PlotlyJS.attr(
        font=axis_font,
        orientation="h",
        y=1.3,
        x=0.0,
        xanchor="left",
        yanchor="top",
        bgcolor="rgba(255,255,255,0.7)",
        bordercolor="black"
    )
    # Custom tickvals and ticktext for hours
    tickvals = [1, 72, 144, 216, 288]
    ticktext = ["00:00", "12:00", "00:00", "12:00", "00:00"]
    # Source current traces (before and after STATCOM)
    min_load = minimum(abs.(load_c), dims=2)[:]
    max_load = maximum(abs.(load_c), dims=2)[:]
    mean_load = ((abs.(load_c[:,1]) .+ abs.(load_c[:,2]) .+ abs.(load_c[:,3])) ./ 3)
    min_vsrc = minimum(abs.(vsource_cg[:,1:3]), dims=2)[:]
    max_vsrc = maximum(abs.(vsource_cg[:,1:3]), dims=2)[:]
    mean_vsrc = ((abs.(vsource_cg[:,1]) .+ abs.(vsource_cg[:,2]) .+ abs.(vsource_cg[:,3])) ./ 3)

    # traces for load (shaded area + mean)
    trace_minmax_load = PlotlyJS.scatter(
        x=time, y=min_load,
        mode="lines",
        fill="none",
        line=attr(color="rgba(31,119,180,0.9)", width=1),
        name="Before STATCOM (min)", showlegend=false
    )
    trace_max_load = PlotlyJS.scatter(
        x=time, y=max_load,
        mode="lines",
        fill="tonexty",
        fillcolor="rgba(31,119,180,0.25)",
        line=attr(color="rgba(31,119,180,0.9)", width=1),
        name="Before STATCOM (range)", showlegend=true
    )
    trace_mean_load = PlotlyJS.scatter(
        x=time, y=mean_load,
        line=attr(color="#1f77b4", width=2),
        name="Mean before", showlegend=true
    )

    # traces for vsource (after STATCOM)
    trace_minmax_vsrc = PlotlyJS.scatter(
        x=time, y=min_vsrc,
        fill="none",
        line=attr(color="rgba(255,127,14,0.9)", width=1),
        name="After STATCOM (min)", showlegend=false
    )
    trace_max_vsrc = PlotlyJS.scatter(
        x=time, y=max_vsrc,
        fill="tonexty",
        fillcolor="rgba(255,127,14,0.4)",
        line=attr(color="rgba(255,127,14,0.9)", width=1),
        name="After STATCOM (range)", showlegend=true
    )
    trace_mean_vsrc = PlotlyJS.scatter(
        x=time, y=mean_vsrc,
        line=attr(color="#ff7f0e", width=1),
        name="Mean after", showlegend=true
    )

    # phase traces for load (phases a,b,c)
    phase_traces_load = [
        PlotlyJS.scatter(x=time, y=abs.(load_c[:,1]), mode="lines", name="Phase a", line=attr(color="#1f77b4", width=1.5), showlegend=true),
        PlotlyJS.scatter(x=time, y=abs.(load_c[:,2]), mode="lines", name="Phase b", line=attr(color="#ff7f0e", width=1.5), showlegend=true),
        PlotlyJS.scatter(x=time, y=abs.(load_c[:,3]), mode="lines", name="Phase c", line=attr(color="#2ca02c", width=1.5), showlegend=true)
    ]

    # every index = 10 minutes (6 steps = 1 hour)
    minutes_per_step = 10

    # Use fixed overlay indices (you can change these indices as needed)
    n = length(time)
    start_idx = clamp(31, 1, n)
    end_idx = clamp(43, 1, n)
    if end_idx < start_idx
        tmp = start_idx; start_idx = end_idx; end_idx = tmp
    end

    # prepare inset (zoom) window data
    time_win = time[start_idx:end_idx]
    min_load_win = min_load[start_idx:end_idx]
    max_load_win = max_load[start_idx:end_idx]
    mean_load_win = mean_load[start_idx:end_idx]
    min_vsrc_win = min_vsrc[start_idx:end_idx]
    max_vsrc_win = max_vsrc[start_idx:end_idx]
    mean_vsrc_win = mean_vsrc[start_idx:end_idx]

    inset_traces = [
        PlotlyJS.scatter(x=time_win, y=min_load_win, mode="lines", line=attr(color="rgba(31,119,180,0.9)", width=1), showlegend=false),
        PlotlyJS.scatter(x=time_win, y=max_load_win, mode="lines", fill="tonexty", fillcolor="rgba(31,119,180,0.25)", line=attr(color="rgba(31,119,180,0.9)", width=1), showlegend=false),
        PlotlyJS.scatter(x=time_win, y=mean_load_win, mode="lines", line=attr(color="#1f77b4", width=2), showlegend=false),
        PlotlyJS.scatter(x=time_win, y=min_vsrc_win, mode="lines", line=attr(color="rgba(255,127,14,0.9)", width=1), showlegend=false),
        PlotlyJS.scatter(x=time_win, y=max_vsrc_win, mode="lines", fill="tonexty", fillcolor="rgba(255,127,14,0.4)", line=attr(color="rgba(255,127,14,0.9)", width=1), showlegend=false),
        PlotlyJS.scatter(x=time_win, y=mean_vsrc_win, mode="lines", line=attr(color="#ff7f0e", width=1), showlegend=false)
    ]

        # inset placement in paper coords
        inset_x0, inset_x1 = 0.74, 0.96
        inset_y0, inset_y1 = 0.75, 0.98
        inset_xc = (inset_x0 + inset_x1) / 2
        inset_yc = (inset_y0 + inset_y1) / 2

        # map time indices to paper x
        x0_p = (start_idx - 1) / max(1, (n - 1))
        x1_p = (end_idx - 1) / max(1, (n - 1))

        margin_before_inset = 0.04
        x1_p_clamped = min(x1_p, inset_x0 - margin_before_inset)
        x0_p_clamped = min(x0_p, x1_p_clamped - 0.01)
        xcenter_p = (x0_p_clamped + x1_p_clamped) / 2

        # map y to paper using global displayed y-range
        y_min_all = minimum(vcat(min_load, min_vsrc))
        y_max_all = maximum(vcat(max_load, max_vsrc))
        yspan = max(1e-6, y_max_all - y_min_all)
        y0_p = (minimum([min_load_win; min_vsrc_win]) - y_min_all) / yspan
        y1_p = (maximum([max_load_win; max_vsrc_win]) - y_min_all) / yspan
        y0_p_c = clamp(y0_p, 0.0, 1.0)
        y1_p_c = clamp(y1_p, 0.0, 1.0)
        ycenter_p = (y0_p_c + y1_p_c) / 2

        # Tick positions for inset: show whole-hour marks within the selected index window
        hour_indices = [i for i in start_idx:end_idx if ((i-1) * minutes_per_step) % 60 == 0]
        if isempty(hour_indices)
            hour_indices = [start_idx, end_idx]
        end
        tick_idx = unique(hour_indices)
        function hhmm_label(idx)
            total_minutes = (idx - 1) * minutes_per_step
            h = Int(fld(total_minutes, 60)) % 24
            m = Int(total_minutes % 60)
            return @sprintf("%02d:%02d", h, m)
        end
        tick_labels = [hhmm_label(i) for i in tick_idx]

        # For scenario 0: do not show inline annotations; show legend labels only.
        if k == 0
            # create a non-legend version of the shaded 'range' trace so it doesn't appear in the legend for sc0
            trace_max_load_sc0 = PlotlyJS.scatter(
                x=time, y=max_load,
                mode="lines",
                fill="tonexty",
                fillcolor="rgba(31,119,180,0.25)",
                line=attr(color="rgba(31,119,180,0.9)", width=1),
                name="", showlegend=false
            )

            # main traces: keep shaded area (no legend) and show only phases + mean (no neutral at all)
            # reverse the legend order by placing the mean trace before the phase traces and reversing the phase traces order
            source_traces = vcat([trace_minmax_load, trace_max_load_sc0, trace_mean_load], reverse(phase_traces_load))

            layout_source = Layout(
                yaxis_title="Feeder F2 Substation Current (A)<br>Before Unbalance Mitigation",
                xaxis_title="Time (h)",
                plot_bgcolor="white",
                paper_bgcolor="white",
                showlegend=true,
                legend=legend_style,
                font=axis_font,
                xaxis=attr(
                    titlefont=axis_font,
                    tickfont=tick_font,
                    tickvals=tickvals,
                    ticktext=ticktext,
                    showline=true, linecolor="black", linewidth=2, mirror=true,
                    showgrid=true, gridcolor="rgba(120,120,120,0.6)", gridwidth=1.2,
                    domain=[0.0, 1.0]
                ),
                yaxis=attr(
                    titlefont=axis_font,
                    tickfont=tick_font,
                    showline=true, linecolor="black", linewidth=2, mirror=true,
                    showgrid=true, gridcolor="rgba(120,120,120,0.6)", gridwidth=1.2,
                    domain=[0.0, 1.0]
                ),
                width=820, height=400
            )

            source_plt = PlotlyJS.Plot(source_traces, layout_source)
            PlotlyJS.savefig(source_plt, "Figures/STATCOM_Isrc_sc$k.pdf")
        else
            # show both before and after, plus inset/overlay
            # ensure before mean and range labels are shown in legend (trace_max_load, trace_mean_load showlegend=true above)
            source_traces = [trace_minmax_load, trace_max_load, trace_mean_load, trace_minmax_vsrc, trace_max_vsrc, trace_mean_vsrc]
            # Layout: main plot and inset axes domains (inset at top-right) with grid on inset
            layout_source = Layout(
                yaxis_title="Feeder F2 Substation Current (A)",
                xaxis_title="Time (h)",
                plot_bgcolor="white",
                paper_bgcolor="white",
                showlegend=true,
                legend=legend_style,
                font=axis_font,
                xaxis=attr(
                    titlefont=axis_font,
                    tickfont=tick_font,
                    tickvals=tickvals,
                    ticktext=ticktext,
                    showline=true, linecolor="black", linewidth=2, mirror=true,
                    showgrid=true, gridcolor="rgba(120,120,120,0.6)", gridwidth=1.2,
                    domain=[0.0, 1.0]
                ),
                yaxis=attr(
                    titlefont=axis_font,
                    tickfont=tick_font,
                    showline=true, linecolor="black", linewidth=2, mirror=true,
                    showgrid=true, gridcolor="rgba(120,120,120,0.6)", gridwidth=1.2,
                    domain=[0.0, 1.0]
                ),
                xaxis2=attr(domain=[inset_x0, inset_x1], anchor="y2", showgrid=true, gridcolor="rgba(120,120,120,0.65)", gridwidth=1.0, zeroline=false, tickvals=tick_idx, ticktext=tick_labels, tickfont=tick_font),
                yaxis2=attr(domain=[inset_y0, inset_y1], anchor="x2", showgrid=true, gridcolor="rgba(120,120,120,0.65)", gridwidth=1.0, zeroline=false, tickfont=tick_font),
                shapes=[
                    Dict(:type => "rect",
                         :xref => "paper", :yref => "paper",
                         :x0 => x0_p_clamped, :x1 => x1_p_clamped,
                         :y0 => clamp(y0_p_c, 0.0, 1.0), :y1 => clamp(y1_p_c, 0.0, 1.0),
                         :line => attr(color="black", width=1),
                         :fillcolor => "rgba(0,0,0,0)"),
                    Dict(:type => "line",
                         :xref => "paper", :yref => "paper",
                         :x0 => xcenter_p, :y0 => ycenter_p,
                         :x1 => inset_x0, :y1 => inset_yc,
                         :line => attr(color="black", width=1, dash="dash"))
                ],
                # no stray annotation text (remove any 'new text' artifacts)
                width=820, height=400
            )

            # combine traces and assign inset traces to x2/y2
            all_traces = vcat(source_traces, inset_traces)
            for i in (length(source_traces)+1):length(all_traces)
                all_traces[i][:xaxis] = "x2"
                all_traces[i][:yaxis] = "y2"
            end

            source_plt = PlotlyJS.Plot(all_traces, layout_source)
            PlotlyJS.savefig(source_plt, "Figures/STATCOM_Isrc_sc$k.pdf")
        end



    # STATCOM current plot (if present)
    if !no_inverter
        phase_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]
        statcom_traces = [
            PlotlyJS.scatter(x=time, y=abs.(inverter_cg[:,1]), name="Phase a", line=attr(color=phase_colors[1], dash="solid")),
            PlotlyJS.scatter(x=time, y=abs.(inverter_cg[:,2]), name="Phase b", line=attr(color=phase_colors[2], dash="solid")),
            PlotlyJS.scatter(x=time, y=abs.(inverter_cg[:,3]), name="Phase c", line=attr(color=phase_colors[3], dash="solid")),
            PlotlyJS.scatter(x=time, y=abs.(inverter_cg[:,4]), name="Neutral", line=attr(color=phase_colors[4], dash="solid")),
            PlotlyJS.scatter(x=time, y=(abs.(inverter_cg[:,1]) .+ abs.(inverter_cg[:,2]) .+ abs.(inverter_cg[:,3])) ./ 3, name="Mean", line=attr(color="black"))
        ]
        layout_statcom = Layout(
            yaxis_title="STATCOM Current (A)",
            xaxis_title="Time (h)",
            plot_bgcolor="white",
            paper_bgcolor="white",
            showlegend=true,
            legend=legend_style,
            font=axis_font,
            xaxis=attr(
                titlefont=axis_font,
                tickfont=tick_font,
                tickvals=tickvals,
                ticktext=ticktext,
                showline=true, linecolor="black", linewidth=2, mirror=true
            ),
            yaxis=attr(
                titlefont=axis_font,
                tickfont=tick_font,
                showline=true, linecolor="black", linewidth=2, mirror=true
            ),
            width=900, height=200
        )
        statcom_plt = PlotlyJS.Plot(statcom_traces, layout_statcom)
        PlotlyJS.savefig(statcom_plt, "Figures/STATCOM_statcom_sc$k.pdf")

        # NOTE: removed individual STATCOM_pdclink_sc$k plot from this function.
        # p_dclink is returned so caller can assemble combined pdclink figure for sc1..sc4.
    else
        statcom_plt = nothing
    end

    # Sequence current plot (negative sequence)
    I_seq_m_load_mn = zeros(3, 0)
    I_seq_m_source_mn = zeros(3, 0)
    for i in 1:n
        _, _, I_seq_m_load = RPMD.get_sequence_components(load_c[i,:])
        _, _, I_seq_m_source = RPMD.get_sequence_components(vsource_cg[i,1:3])
        I_seq_m_load_mn = hcat(I_seq_m_load_mn, I_seq_m_load)
        I_seq_m_source_mn = hcat(I_seq_m_source_mn, I_seq_m_source)
    end
    trace_negseq = PlotlyJS.scatter(x=time, y=I_seq_m_source_mn[3,:], name="Sc. $k", line=attr(color="red"))
    layout_negseq = Layout(
        yaxis_title="Negative Sequence Current (A)",
        xaxis_title="Time (h)",
        plot_bgcolor="white",
        paper_bgcolor="white",
        showlegend=true,
        legend=legend_style,
        font=axis_font,
        xaxis=attr(
            titlefont=axis_font,
            tickfont=tick_font,
            tickvals=tickvals,
            ticktext=ticktext,
            showline=true, linecolor="black", linewidth=2, mirror=true
        ),
        yaxis=attr(
            titlefont=axis_font,
            tickfont=tick_font,
            showline=true, linecolor="black", linewidth=2, mirror=true
        ),
        width=900, height=200
    )
    negseq_plt = Plot([trace_negseq], layout_negseq)
    PlotlyJS.savefig(negseq_plt, "Figures/STATCOM_load_Iseq_sc$k.pdf")

    return (source_plt, statcom_plt, negseq_plt, p_dclink)
end

# Example usage:
Iseq_ldc_plt = nothing
Iseq_plt = nothing
currents_sc_plt = Dict()
for (k, result) in enumerate([result_conv_sc0_mn, result_conv_sc1_mn, result_conv_sc2_mn, result_conv_sc3_mn, result_conv_sc4_mn])
    if k==1 # no inverter => no compensation
        currents_sc_plt[k] = plot_results_plotlyjs(result, k-1, timesteps_plt; Iseq_plt=Iseq_plt, ldc_plt=Iseq_ldc_plt, no_inverter=true)
    else
        currents_sc_plt[k] = plot_results_plotlyjs(result, k-1, timesteps_plt; Iseq_plt=Iseq_plt, ldc_plt=Iseq_ldc_plt)
    end
end

# Build combined pdclink plot for scenarios 1..4 (sc1..sc4) on the same figure
# colors consistent with other inverter plots
pd_colors = ["#0074D9", "#FF4136", "#7BD389", "#C49BF0"]  # more vibrant blue and red
pd_traces = PlotlyJS.GenericTrace{Dict{Symbol,Any}}[]
pd_data = Vector{Vector{Float64}}()
for i in 1:4
    # currents_sc_plt index mapping: 1=>sc0, 2=>sc1, ... so scenario i corresponds to index i+1
    entry = currents_sc_plt[i+1]
    if entry === nothing
        continue
    end
    p = entry[4]  # p_dclink returned as 4th element
    p_vec = collect(Float64.(p))
    push!(pd_data, p_vec)
    t = 1:length(p_vec)
    push!(pd_traces, PlotlyJS.scatter(x=t, y=p_vec, mode="lines", name=all_labels[i], line=attr(color=pd_colors[i], width=2)))
end

# determine y-limits to cover all traces with a small padding
if !isempty(pd_data)
    y_min = minimum([minimum(v) for v in pd_data])
    y_max = maximum([maximum(v) for v in pd_data])
    span = max(1e-6, y_max - y_min)
    pad = 0.05 * span
    y_lo = y_min - pad
    y_hi = y_max + pad
else
    y_lo = 0.0
    y_hi = 1.0
end
layout_pd = Layout(
    yaxis_title="2w DC Link Power (kW)",
    # left y-axis: visible line, grid and ticks
    yaxis=PlotlyJS.attr(
        range=[y_lo, y_hi],
        showgrid=true, gridcolor="rgba(120,120,120,0.25)", gridwidth=1,
        showline=true, linecolor="black", linewidth=2, mirror=false,
        ticks="outside",
        # larger title and tick fonts
        titlefont=PlotlyJS.attr(size=22, family="serif"),
        tickfont=PlotlyJS.attr(size=18, family="serif")
    ),
    # add a mirrored right y-axis to create a box with y axis on both sides
    yaxis2=PlotlyJS.attr(
        overlaying="y",
        side="right",
        showgrid=false,
        showline=true, linecolor="black", linewidth=2,
        ticks="outside",
        titlefont=PlotlyJS.attr(size=22, family="serif"),
        tickfont=PlotlyJS.attr(size=18, family="serif")
    ),
    # make xticks show the same labels defined in `ticktext` and draw top/bottom lines
    xaxis=PlotlyJS.attr(
        title="Time (h)",
        tickmode="array",
        tickvals=tickvals,
        ticktext=ticktext,
        showline=true, linecolor="black", linewidth=2, mirror=true,
        showgrid=true, gridcolor="rgba(120,120,120,0.25)", gridwidth=1,
        titlefont=PlotlyJS.attr(size=22, family="serif"),
        tickfont=PlotlyJS.attr(size=18, family="serif")
    ),
    plot_bgcolor="white",
    paper_bgcolor="white",
    showlegend=true,
    # Increase legend font size and move legend slightly above the figure
    legend=PlotlyJS.attr(orientation="h", x=0.0, y=1.05, yanchor="bottom", font=PlotlyJS.attr(size=18, family="serif")),
    # keep global font a bit smaller than axis titles
    font=PlotlyJS.attr(family="serif", size=16),
    width=900,
    height=520,  # make the plot taller
    # add an explicit vertical line at the far right of the plotting area (paper coords)
    shapes=[
        Dict(
            :type => "line",
            :xref => "paper", :yref => "paper",
            :x0 => 1.0, :x1 => 1.0,
            :y0 => 0.0, :y1 => 1.0,
            :line => PlotlyJS.attr(color="black", width=2)
        )
    ]
)
pdclink_combined = Plot(pd_traces, layout_pd)
PlotlyJS.savefig(pdclink_combined, "Figures/STATCOM_pdclink_sc1_sc4_combined.pdf")


##
load_c_sc1, vsource_cg_sc1, inverter_cg_sc1, p_dclink_sc1 = get_results(result_conv_sc1_mn, timesteps_plt, no_inverter=false)
load_seq1 = zeros(3)
inverter_seq1 = zeros(3)
vsource_seq1 = zeros(3)
for i in timesteps_plt
    _, _, I_seq_m = RPMD.get_sequence_components(load_c_sc1[i,:]);      load_seq1 = hcat(load_seq1, I_seq_m)
    _, _, I_seq_m = RPMD.get_sequence_components(vsource_cg_sc1[i,1:3]);  vsource_seq1 = hcat(vsource_seq1, I_seq_m)
    _, _, I_seq_m = RPMD.get_sequence_components(inverter_cg_sc1[i,1:3]); inverter_seq1 = hcat(inverter_seq1, I_seq_m)
end
statcom_plt1 = Plots.plot(abs.(inverter_cg_sc1[:,1]), color=1, linewidth=1, linestyle=:dash, label="phase a")
Plots.plot!(abs.(inverter_cg_sc1[:,2]), color=2, linewidth=1, linestyle=:dash, label="phase b", xtick=false)
Plots.plot!(abs.(inverter_cg_sc1[:,3]), color=3, linewidth=1, linestyle=:dash, label="phase c", ylim=(0,100))
Plots.plot!(abs.(inverter_cg_sc1[:,4]), color=4, linewidth=1, linestyle=:dash, label="neutral", title="STATCOM Currents (A)", ylabel=all_labels[1])
Plots.plot!((abs.(inverter_cg_sc1[:,1])+abs.(inverter_cg_sc1[:,2])+abs.(inverter_cg_sc1[:,3]))/3, color=:black, linewidth=1, label=false, size=(600, 200))
# ylabel!(L"|I_{statcom}|")

load_c_sc1, vsource_cg_sc1, inverter_cg_sc1, p_dclink_sc1 = get_results(result_conv_sc1_mn, timesteps_plt, no_inverter=false)
load_seq1 = zeros(3)
inverter_seq1 = zeros(3)
vsource_seq1 = zeros(3)
for i in timesteps_plt
    _, _, I_seq_m = RPMD.get_sequence_components(load_c_sc1[i,:]);      load_seq1 = hcat(load_seq1, I_seq_m)
    _, _, I_seq_m = RPMD.get_sequence_components(vsource_cg_sc1[i,1:3]);  vsource_seq1 = hcat(vsource_seq1, I_seq_m)
    _, _, I_seq_m = RPMD.get_sequence_components(inverter_cg_sc1[i,1:3]); inverter_seq1 = hcat(inverter_seq1, I_seq_m)
end
statcom_plt1 = Plots.plot(abs.(inverter_cg_sc1[:,1]), color=1, linewidth=1, linestyle=:dash, label="phase a")
Plots.plot!(abs.(inverter_cg_sc1[:,2]), color=2, linewidth=1, linestyle=:dash, label="phase b", xtick=false)
Plots.plot!(abs.(inverter_cg_sc1[:,3]), color=3, linewidth=1, linestyle=:dash, label="phase c", ylim=(0,100))
Plots.plot!(abs.(inverter_cg_sc1[:,4]), color=4, linewidth=1, linestyle=:dash, label="neutral", title="STATCOM Currents (A)", ylabel=all_labels[1])
Plots.plot!((abs.(inverter_cg_sc1[:,1])+abs.(inverter_cg_sc1[:,2])+abs.(inverter_cg_sc1[:,3]))/3, color=:black, linewidth=1, label=false, size=(600, 200))
# ylabel!(L"|I_{statcom}|")

load_c_sc2, vsource_cg_sc2, inverter_cg_sc2, p_dclink_sc2 = get_results(result_conv_sc2_mn, timesteps_plt, no_inverter=false)
load_seq2 = zeros(3)
inverter_seq2 = zeros(3)
vsource_seq2 = zeros(3)
for i in timesteps_plt
    _, _, I_seq_m = RPMD.get_sequence_components(load_c_sc2[i,:]);      load_seq2 = hcat(load_seq2, I_seq_m)
    _, _, I_seq_m = RPMD.get_sequence_components(vsource_cg_sc2[i,1:3]);  vsource_seq2 = hcat(vsource_seq2, I_seq_m)
    _, _, I_seq_m = RPMD.get_sequence_components(inverter_cg_sc2[i,1:3]); inverter_seq2 = hcat(inverter_seq2, I_seq_m)
end
statcom_plt2 = Plots.plot(abs.(inverter_cg_sc2[:,1]), color=1, linewidth=1, linestyle=:dash, label=false)
Plots.plot!(abs.(inverter_cg_sc2[:,2]), color=2, linewidth=1, linestyle=:dash, label=false, xticks=false)
Plots.plot!(abs.(inverter_cg_sc2[:,3]), color=3, linewidth=1, linestyle=:dash, label=false, ylim=(0,100))
Plots.plot!(abs.(inverter_cg_sc2[:,4]), color=4, linewidth=1, linestyle=:dash, label=false, ylabel=all_labels[2])
Plots.plot!((abs.(inverter_cg_sc2[:,1])+abs.(inverter_cg_sc2[:,2])+abs.(inverter_cg_sc2[:,3]))/3, color=:black, linewidth=1, label=false, size=(600, 200))
# ylabel!(L"|I_{statcom}|")

load_c_sc3, vsource_cg_sc3, inverter_cg_sc3, p_dclink_sc3 = get_results(result_conv_sc3_mn, timesteps_plt, no_inverter=false)
load_seq3 = zeros(3)
vsource_seq3 = zeros(3)
inverter_seq3 = zeros(3)
for i in timesteps_plt
    _, _, I_seq_m = RPMD.get_sequence_components(load_c_sc3[i,:]);      load_seq3 = hcat(load_seq3, I_seq_m)
    _, _, I_seq_m = RPMD.get_sequence_components(vsource_cg_sc3[i,1:3]);  vsource_seq3 = hcat(vsource_seq3, I_seq_m)
    _, _, I_seq_m = RPMD.get_sequence_components(inverter_cg_sc3[i,1:3]); inverter_seq3 = hcat(inverter_seq3, I_seq_m)
end
statcom_plt3 = Plots.plot(abs.(inverter_cg_sc3[:,1]), color=1, linewidth=1, linestyle=:dash, label=false)
Plots.plot!(abs.(inverter_cg_sc3[:,2]), color=2, linewidth=1, linestyle=:dash, label=false, xticks=false)
Plots.plot!(abs.(inverter_cg_sc3[:,3]), color=3, linewidth=1, linestyle=:dash, label=false, ylim=(0,100))
Plots.plot!(abs.(inverter_cg_sc3[:,4]), color=4, linewidth=1, linestyle=:dash, label=false, ylabel=all_labels[3])
Plots.plot!((abs.(inverter_cg_sc3[:,1])+abs.(inverter_cg_sc3[:,2])+abs.(inverter_cg_sc3[:,3]))/3, color=:black, linewidth=1, label=false, size=(600, 200))
# ylabel!(L"|I_{statcom}|")

load_c_sc4, vsource_cg_sc4, inverter_cg_sc4, p_dclink_sc4 = get_results(result_conv_sc4_mn, timesteps_plt, no_inverter=false)
load_seq4 = zeros(3)
vsource_seq4 = zeros(3)
inverter_seq4 = zeros(3)
for i in timesteps_plt
    _, _, I_seq_m = RPMD.get_sequence_components(load_c_sc4[i,:]);        load_seq4 = hcat(load_seq4, I_seq_m)
    _, _, I_seq_m = RPMD.get_sequence_components(vsource_cg_sc4[i,1:3]);  vsource_seq4 = hcat(vsource_seq4, I_seq_m)
    _, _, I_seq_m = RPMD.get_sequence_components(inverter_cg_sc4[i,1:3]); inverter_seq4 = hcat(inverter_seq4, I_seq_m)
end
statcom_plt4 = Plots.plot(abs.(inverter_cg_sc4[:,1]), color=1, linewidth=1, linestyle=:dash, label=false)
Plots.plot!(abs.(inverter_cg_sc4[:,2]), color=2, linewidth=1, linestyle=:dash, label=false, ylim=(0,100))
Plots.plot!(abs.(inverter_cg_sc4[:,3]), color=3, linewidth=1, linestyle=:dash, label=false)
Plots.plot!(abs.(inverter_cg_sc4[:,4]), color=4, linewidth=1, linestyle=:dash, label=false, ylabel=all_labels[4])
Plots.plot!((abs.(inverter_cg_sc4[:,1])+abs.(inverter_cg_sc4[:,2])+abs.(inverter_cg_sc4[:,3]))/3, color=:black, linewidth=1, label=false, size=(600, 200))
# ylabel!(L"|I_{statcom}|")

statcom_plt = Plots.plot(statcom_plt1, statcom_plt2, statcom_plt3, statcom_plt4, layout=(4,1), size=(1000,1000))
Plots.savefig(statcom_plt, "Figures/Scenarios_phase_currents.pdf")


# Define consistent colors for a, b, c, n
phase_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]  # blue, orange, green, purple

# Redefine traces for SC1 (solid lines, legend shown)
trace_sc1_a = PlotlyJS.scatter(y=abs.(inverter_cg_sc1[:,1]), mode="lines", name="Phase a", line=attr(color=phase_colors[1], dash="solid"), showlegend=true)
trace_sc1_b = PlotlyJS.scatter(y=abs.(inverter_cg_sc1[:,2]), mode="lines", name="Phase b", line=attr(color=phase_colors[2], dash="solid"), showlegend=true)
trace_sc1_c = PlotlyJS.scatter(y=abs.(inverter_cg_sc1[:,3]), mode="lines", name="Phase c", line=attr(color=phase_colors[3], dash="solid"), showlegend=true)
trace_sc1_n = PlotlyJS.scatter(y=abs.(inverter_cg_sc1[:,4]), mode="lines", name="Neutral", line=attr(color=phase_colors[4], dash="solid"), showlegend=true)

# Redefine traces for SC2 (solid lines, legend shown)
trace_sc2_a = PlotlyJS.scatter(y=abs.(inverter_cg_sc2[:,1]), mode="lines", name="Phase a", line=attr(color=phase_colors[1], dash="solid"), showlegend=true)
trace_sc2_b = PlotlyJS.scatter(y=abs.(inverter_cg_sc2[:,2]), mode="lines", name="Phase b", line=attr(color=phase_colors[2], dash="solid"), showlegend=true)
trace_sc2_c = PlotlyJS.scatter(y=abs.(inverter_cg_sc2[:,3]), mode="lines", name="Phase c", line=attr(color=phase_colors[3], dash="solid"), showlegend=true)
trace_sc2_n = PlotlyJS.scatter(y=abs.(inverter_cg_sc2[:,4]), mode="lines", name="Neutral", line=attr(color=phase_colors[4], dash="solid"), showlegend=true)

# Redefine traces for SC3 (solid lines, legend hidden)
trace_sc3_a = PlotlyJS.scatter(y=abs.(inverter_cg_sc3[:,1]), mode="lines", name="Phase a", line=attr(color=phase_colors[1], dash="solid"), showlegend=false)
trace_sc3_b = PlotlyJS.scatter(y=abs.(inverter_cg_sc3[:,2]), mode="lines", name="Phase b", line=attr(color=phase_colors[2], dash="solid"), showlegend=false)
trace_sc3_c = PlotlyJS.scatter(y=abs.(inverter_cg_sc3[:,3]), mode="lines", name="Phase c", line=attr(color=phase_colors[3], dash="solid"), showlegend=false)
trace_sc3_n = PlotlyJS.scatter(y=abs.(inverter_cg_sc3[:,4]), mode="lines", name="Neutral", line=attr(color=phase_colors[4], dash="solid"), showlegend=false)

sc1_traces = [trace_sc1_a, trace_sc1_b, trace_sc1_c, trace_sc1_n]
sc2_traces = [trace_sc2_a, trace_sc2_b, trace_sc2_c, trace_sc2_n]
sc3_traces = [trace_sc3_a, trace_sc3_b, trace_sc3_c, trace_sc3_n]
# Define custom tickvals and ticktext for hours


# Compute global y-limits for both subplots
# all_y = vcat(abs.(inverter_cg_sc2), abs.(inverter_cg_sc3))
all_y = vcat(abs.(inverter_cg_sc1), abs.(inverter_cg_sc3))
ymin = minimum(all_y)
ymax = maximum(all_y)

# Define LaTeX-like font
latex_font = attr(
    family="Computer Modern, serif",
    size=18,
    color="black"
)
plt_subplots = PlotlyJS.Plot(
    vcat(sc1_traces, sc3_traces),
    Layout(
        grid=attr(rows=2, columns=1, pattern="independent"),
        width=900,
        height=600,
        plot_bgcolor="white",
        paper_bgcolor="white",
        showlegend=true,
        legend=attr(
            orientation="h",
            x=0.0,
            y=1.12,
            font=attr(size=22, family="serif")
        ),
        font=attr(family="serif", size=22),
        yaxis=attr(
            title="Case 1a)<br>uncst. 2w & neutral",
            range=[ymin, ymax],
            titlefont=attr(size=22, family="serif"),
            tickfont=attr(size=18, family="serif"),
            showline=true, linewidth=0.1, linecolor="grey", mirror=true,
            showgrid=true, gridcolor="rgba(120,120,120,0.45)", gridwidth=1    # darker grid
        ),
        yaxis2=attr(
            title="Case 1c)<br>no 2w & uncst. neutral",
            range=[ymin, ymax],
            titlefont=attr(size=22, family="serif"),
            tickfont=attr(size=18, family="serif"),
            showline=true, linewidth=0.1, linecolor="grey", mirror=true,
            showgrid=true, gridcolor="rgba(120,120,120,0.45)", gridwidth=1    # darker grid
        ),
        xaxis=attr(
            title="Time (h)",
            tickmode="array",
            tickvals=tickvals,
            ticktext=ticktext,
            range=[1, tickvals[end]],    # ensure last tick (00:00) is shown
            titlefont=attr(size=22, family="serif"),
            tickfont=attr(size=18, family="serif"),
            showline=true, linewidth=0.1, linecolor="grey", mirror=true,
            showgrid=true, gridcolor="rgba(120,120,120,0.45)", gridwidth=1    # darker grid
        ),
        xaxis2=attr(
            title="Time (h)",
            tickmode="array",
            tickvals=tickvals,
            ticktext=ticktext,
            range=[1, tickvals[end]],    # ensure last tick (00:00) is shown
            titlefont=attr(size=22, family="serif"),
            tickfont=attr(size=18, family="serif"),
            showline=true, linewidth=0.1, linecolor="grey", mirror=true,
            showgrid=true, gridcolor="rgba(120,120,120,0.45)", gridwidth=1    # darker grid
        )
    )
)
for i in 1:4
    plt_subplots.data[i][:xaxis] = "x"
    plt_subplots.data[i][:yaxis] = "y"
end
for i in 5:8
    plt_subplots.data[i][:xaxis] = "x2"
    plt_subplots.data[i][:yaxis] = "y2"
end
plt_subplots

PlotlyJS.savefig(plt_subplots, "Figures/STATCOM_1a_1c_inverter_currents.pdf")

##
current_pos_seq_plt = Plots.plot(load_seq4[2,2:end], label="load", title="Sequence Current Contribution by STATCOM (A)")
Plots.plot!(inverter_seq1[2,2:end], label="SC. 1", ylabel="Positive", xticks=false)
Plots.plot!(inverter_seq2[2,2:end], label="SC. 2")
Plots.plot!(inverter_seq3[2,2:end], label="SC. 3")
Plots.plot!(inverter_seq4[2,2:end], label="SC. 4")

current_neg_seq_plt = Plots.plot(load_seq3[3,2:end], label="load")
Plots.plot!(inverter_seq1[3,2:end], label="SC. 1", ylabel="Negative", xticks=false)
Plots.plot!(inverter_seq2[3,2:end], label="SC. 2")
Plots.plot!(inverter_seq3[3,2:end], label="SC. 3")
Plots.plot!(inverter_seq4[3,2:end], label="SC. 4", legend=false)

current_zero_seq_plt = Plots.plot(load_seq3[1,2:end], label="load")
Plots.plot!(inverter_seq1[1,2:end], label="SC. 1", ylabel="Zero")
Plots.plot!(inverter_seq2[1,2:end], label="SC. 2")
Plots.plot!(inverter_seq3[1,2:end], label="SC. 3")
Plots.plot!(inverter_seq4[1,2:end], label="SC. 4", legend=false, xticks=(1:72:48*4, string.(1:72:48*4)))

current_seq_plt = Plots.plot(current_pos_seq_plt, current_neg_seq_plt, current_zero_seq_plt, layout=(3,1))
Plots.savefig(current_seq_plt, "Figures/Scenarios_sequence_currents.pdf")




##
using StatsPlots
using PlotlyJS


load_c_sc1, vsource_cg_sc1, inverter_cg_sc1, p_dclink_sc1 = get_results(result_conv_sc1_mn, timesteps, no_inverter=false)
load_c_sc2, vsource_cg_sc2, inverter_cg_sc2, p_dclink_sc2 = get_results(result_conv_sc2_mn, timesteps, no_inverter=false)
load_c_sc3, vsource_cg_sc3, inverter_cg_sc3, p_dclink_sc3 = get_results(result_conv_sc3_mn, timesteps, no_inverter=false)
load_c_sc4, vsource_cg_sc4, inverter_cg_sc4, p_dclink_sc4 = get_results(result_conv_sc4_mn, timesteps, no_inverter=false)
inverter_cga = [abs.(inverter_cg_sc1[:,1]) abs.(inverter_cg_sc2[:,1]) abs.(inverter_cg_sc3[:,1]) abs.(inverter_cg_sc4[:,1])]
inverter_cgb = [abs.(inverter_cg_sc1[:,2]) abs.(inverter_cg_sc2[:,2]) abs.(inverter_cg_sc3[:,2]) abs.(inverter_cg_sc4[:,2])]
inverter_cgc = [abs.(inverter_cg_sc1[:,3]) abs.(inverter_cg_sc2[:,3]) abs.(inverter_cg_sc3[:,3]) abs.(inverter_cg_sc4[:,3])]
inverter_cgn = [abs.(inverter_cg_sc1[:,4]) abs.(inverter_cg_sc2[:,4]) abs.(inverter_cg_sc3[:,4]) abs.(inverter_cg_sc4[:,4])]

load_seq1 = zeros(3)
vsource_seq1 = zeros(3)
inverter_seq1 = zeros(3)
load_seq2 = zeros(3)
vsource_seq2 = zeros(3)
inverter_seq2 = zeros(3)
load_seq3 = zeros(3)
vsource_seq3 = zeros(3)
inverter_seq3 = zeros(3)
load_seq4 = zeros(3)
vsource_seq4 = zeros(3)
inverter_seq4 = zeros(3)
for i in timesteps
    _, _, seq_m = RPMD.get_sequence_components(load_c_sc1[i,:]);        load_seq1 = hcat(load_seq1, seq_m)
    _, _, seq_m = RPMD.get_sequence_components(vsource_cg_sc1[i,1:3]);  vsource_seq1 = hcat(vsource_seq1, seq_m)
    _, _, seq_m = RPMD.get_sequence_components(inverter_cg_sc1[i,1:3]); inverter_seq1 = hcat(inverter_seq1, seq_m)
    _, _, seq_m = RPMD.get_sequence_components(load_c_sc2[i,:]);        load_seq2 = hcat(load_seq2, seq_m)
    _, _, seq_m = RPMD.get_sequence_components(vsource_cg_sc2[i,1:3]);  vsource_seq2 = hcat(vsource_seq2, seq_m)
    _, _, seq_m = RPMD.get_sequence_components(inverter_cg_sc2[i,1:3]); inverter_seq2 = hcat(inverter_seq2, seq_m)
    _, _, seq_m = RPMD.get_sequence_components(load_c_sc3[i,:]);        load_seq3 = hcat(load_seq3, seq_m)
    _, _, seq_m = RPMD.get_sequence_components(vsource_cg_sc3[i,1:3]);  vsource_seq3 = hcat(vsource_seq3, seq_m)
    _, _, seq_m = RPMD.get_sequence_components(inverter_cg_sc3[i,1:3]); inverter_seq3 = hcat(inverter_seq3, seq_m)
    _, _, seq_m = RPMD.get_sequence_components(load_c_sc4[i,:]);        load_seq4 = hcat(load_seq4, seq_m)
    _, _, seq_m = RPMD.get_sequence_components(vsource_cg_sc4[i,1:3]);  vsource_seq4 = hcat(vsource_seq4, seq_m)
    _, _, seq_m = RPMD.get_sequence_components(inverter_cg_sc4[i,1:3]); inverter_seq4 = hcat(inverter_seq4, seq_m)
end


all_labels = ["Case 1a) uncst. 2w & neutral", "Case 1b) uncst. 2w & no neutral", "Case 1c) no 2w & uncst. neutral", "Case 1d) cst. 2w & neutral", "Non-Mitigated F2 Substation Current"]
y_zero = [inverter_seq1[1,2:end] inverter_seq2[1,2:end] inverter_seq3[1,2:end] inverter_seq4[1,2:end] load_seq4[1,2:end]]
y_neg = [inverter_seq1[3,2:end] inverter_seq2[3,2:end] inverter_seq3[3,2:end] inverter_seq4[3,2:end] load_seq4[3,2:end]]

# Use multi-line x-labels with \n and slightly smaller legend font to encourage wrapping into three rows
xlabels = ["Zero sequence (proportional<br>to neutral current)", "Negative sequence (approx.<br>proportional to 2w ripple power)"]
boxplots_Iseqs = [
    box(
        x = repeat(xlabels, 8806),
        y = vec(vcat(vec(y_zero[:,j])', vec(y_neg[:,j])')),
        name = all_labels[j]
    ) for j in 1:5
]
plt_Iseqs = Plot(
    boxplots_Iseqs,
    Layout(
        boxmode = "group",
        yaxis_title = "STATCOM or Substation Current (A)",
        width = 700,
        height = 440,
        plot_bgcolor = "white",
        paper_bgcolor = "white",
        legend = attr(
            orientation = "h",
            x = 0.0,
            y = 1.28,
            xanchor = "left",
            yanchor = "top",
            bgcolor = "rgba(255,255,255,0.7)",
            bordercolor = "black",
            font = attr(size = 16, family = "serif")
        ),
        xaxis = attr(
            showline = true,
            linecolor = "black",
            linewidth = 2,
            mirror = true,
            titlefont = attr(size = 20, family = "serif"),   # increased xlabel font
            tickfont = attr(size = 20, family = "serif")
        ),
        yaxis = attr(
            showline = true,
            linecolor = "black",
            linewidth = 2,
            mirror = true,
            titlefont = attr(size = 20, family = "serif"),   # increased ylabel font
            tickfont = attr(size = 18, family = "serif"),
            showgrid = true,
            gridcolor = "rgba(120,120,120,0.25)",
            gridwidth = 1,
            # Make the zero (zeroline) identical to the horizontal grid lines
            zeroline = true,
            zerolinecolor = "rgba(120,120,120,0.25)",
            zerolinewidth = 1
        ),
        font = attr(family = "serif", size = 18)
    )
)
PlotlyJS.savefig(plt_Iseqs, "Figures/STATCOM_Iseqs_boxplot.pdf")

labels = ["Case 1a) uncst. 2w & neutral", "Case 1b) uncst. 2w & no neutral", "Case 1c) no 2w & uncst. neutral", "Case 1d) cst. 2w & neutral"]
pdclink = [p_dclink_sc1 p_dclink_sc2 p_dclink_sc3 p_dclink_sc4]
boxplots_pdc=[box(x=repeat(["2w and Neutral Scenarios"], 768), y=vec(pdclink[:,j]), name=labels[j]) for j in 1:4]
# boxplots_pdc=[box(x=repeat(["2w and Neutral Scenarios"], 768), y=vec(pdclink[:,j]), name="") for j in 1:4]
plt_pdclink = Plot(
    boxplots_pdc,
    Layout(
        boxmode="group",
        yaxis_title="2w DC Link Power (kW)",
        width=500,
        height=400,
        plot_bgcolor="white",
        paper_bgcolor="white",
        showlegend=true,
        legend=attr(
            orientation="h",
            x=0.0,
            y=1.3,    # shifted higher like the reference layout
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(255,255,255,0.7)",
            bordercolor="black",
            font=attr(size=22, family="serif")
        ),
        xaxis=attr(
            showline=true,
            linecolor="black",
            linewidth=2,
            mirror=true,
            titlefont=attr(size=22, family="serif"),
            tickfont=attr(size=22, family="serif")
        ),
        yaxis=attr(
            showline=true,
            linecolor="black",
            linewidth=2,
            mirror=true,
            titlefont=attr(size=22, family="serif"),
            tickfont=attr(size=18, family="serif")
        ),
        font=attr(family="serif", size=22)
    )
)
PlotlyJS.savefig(plt_pdclink, "Figures/STATCOM_Pdclink_boxplot.pdf")




v_zero = [vsource_seq1[1,2:end] vsource_seq2[1,2:end] vsource_seq3[1,2:end] vsource_seq4[1,2:end]]
v_neg = [vsource_seq1[3,2:end] vsource_seq2[3,2:end] vsource_seq3[3,2:end] vsource_seq4[3,2:end]]
boxplots_Vseqs=[box(x=repeat(["Zero Sequence", "Negative Sequence"], 8806), y=vec(vcat(vec(v_zero[:,j])', vec(v_neg[:,j])')), name=labels[j]) for j in 1:4]
plt_Vseqs = Plot(
    boxplots_Vseqs,
    Layout(
        boxmode="group",
        yaxis_title="Voltage (V)",
        width=500,
        height=400,
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=attr(
            orientation="h",
            x=0.0,
            y=1.3,  # shifted higher
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(255,255,255,0.7)",
            bordercolor="black",
            font=attr(size=22, family="serif")
        ),
        xaxis=attr(
            showline=true,
            linecolor="black",
            linewidth=2,
            mirror=true,
            titlefont=attr(size=22, family="serif"),
            tickfont=attr(size=22, family="serif")
        ),
        yaxis=attr(
            showline=true,
            linecolor="black",
            linewidth=2,
            mirror=true,
            titlefont=attr(size=22, family="serif"),
            tickfont=attr(size=18, family="serif")
        ),
        font=attr(family="serif", size=22)
    )
)
PlotlyJS.savefig(plt_Vseqs, "Figures/STATCOM_Vseqs_boxplot.pdf")


# ldc_plt = plot()
# for (k, result) in enumerate([result_conv_sc0_mn result_conv_sc1_mn result_conv_sc2_mn result_conv_sc3_mn result_conv_sc4_mn])
#     # for (k, result) in enumerate([result_conv_sc2_mn])
#     markershape = :none
#     if k == 1
#         markershape = :circle
#     end

#     branch_cr_1, branch_ci_1, branch_c_1 = get_currents(result, "branch", "1", 1, timesteps)
#     branch_cr_2, branch_ci_2, branch_c_2 = get_currents(result, "branch", "1", 2, timesteps)
#     branch_cr_3, branch_ci_3, branch_c_3 = get_currents(result, "branch", "1", 3, timesteps)
#     branch_cr_4, branch_ci_4, branch_c_4 = get_currents(result, "branch", "1", 4, timesteps)
    
#     #################### plot load duration curves
#     data = maximum(abs.([branch_c_1 branch_c_2 branch_c_3 branch_c_4]), dims=2)
#     plot_load_duration_curve!(data; plt=ldc_plt, label="Sc $(k-1)", markershape=markershape, xlabel="Time Index (10-minutely)", ylabel="Maximum Source Current (A)")
# end
# display(ldc_plt)
# savefig(ldc_plt, "Figures/STATCOM_Isrc_max_LDC_2.pdf")
## Prepare data for PlotlyJS
traces = GenericTrace{Dict{Symbol, Any}}[]
tail_traces = GenericTrace{Dict{Symbol, Any}}[]

# More distinct color palette (Tableau / Plotly distinct colors)
colors = ["#1f77b4", "#ffd700", "#2ca02c", "#b30000", "#9467bd"]

for (k, result) in enumerate([result_conv_sc1_mn, result_conv_sc2_mn, result_conv_sc3_mn, result_conv_sc4_mn, result_conv_sc0_mn])
    branch_cr_1, branch_ci_1, branch_c_1 = get_currents(result, "branch", "1", 1, timesteps)
    branch_cr_2, branch_ci_2, branch_c_2 = get_currents(result, "branch", "1", 2, timesteps)
    branch_cr_3, branch_ci_3, branch_c_3 = get_currents(result, "branch", "1", 3, timesteps)
    branch_cr_4, branch_ci_4, branch_c_4 = get_currents(result, "branch", "1", 4, timesteps)

    # Compute max source current per timestep
    data = maximum(abs.([branch_c_1 branch_c_2 branch_c_3 branch_c_4]), dims=2)
    data_vec = vec(data)

    # Sort for duration curve (descending)
    sorted_data = sort(data_vec, rev=true)
    n = length(sorted_data)
    time_percent = (1:n) ./ n * 100  # Percentage of time

    # Main trace (duration curve) - use consistent distinct color
    push!(traces, PlotlyJS.scatter(x=time_percent, y=sorted_data, mode="lines",
                                   name=all_labels[k],
                                   line=attr(color=colors[k], width=2.2)))

    # Tail trace for zoomed-in section (first 300 sorted points)
    tail_idx = 1:min(300, n)
    # normalized x for inset 0..1
    tail_x = length(tail_idx) > 1 ? collect(0:(1/(length(tail_idx)-1)):1) : [0.0]
    # prepare y for inset, replace values < 0 with NaN (so they are not plotted)
    y_tail = copy(sorted_data[tail_idx])
    y_tail[y_tail .< 0.0] .= NaN

    # store tail traces (no legend)
    push!(tail_traces, PlotlyJS.scatter(x=tail_x, y=y_tail, mode="lines", name="",
                                        line=attr(color=colors[k], width=2.2), showlegend=false))
end

# Main plot
layout_main = Layout(
    width=900,
    height=560,
    plot_bgcolor="white",
    paper_bgcolor="white",
    legend=attr(
        orientation="h",
        x=0.0,
        y=1.22,
        xanchor="left",
        yanchor="top",
        bgcolor="rgba(255,255,255,0.8)",
        bordercolor="black",
        font=attr(size=18, family="serif")
    ),
    xaxis=attr(
        title="Load percentile (%)",
        showline=true, linecolor="black", linewidth=2, mirror=true,
        titlefont=attr(size=20, family="serif"), tickfont=attr(size=16, family="serif"),
        showgrid=true, gridcolor="rgba(120,120,120,0.25)", gridwidth=1
    ),
    yaxis=attr(
        title="Maximum Feeder F2 Substation<br>Phase Current (A)",
        showline=true, linecolor="black", linewidth=2, mirror=true,
        titlefont=attr(size=20, family="serif"), tickfont=attr(size=16, family="serif"),
        showgrid=true, gridcolor="rgba(120,120,120,0.25)", gridwidth=1
    ),
    font=attr(family="serif", size=18)
)
plt_main = PlotlyJS.Plot(traces, layout_main)
PlotlyJS.savefig(plt_main, "Figures/STATCOM_Isrc_max_LDC_2_plotly.pdf")

# Zoomed-in tail plot
layout_tail = Layout(
    width=900,
    height=560,
    plot_bgcolor="white",
    paper_bgcolor="white",
    legend=attr(
        orientation="h",
        x=0.0,
        y=1.22,
        xanchor="left",
        yanchor="top",
        bgcolor="rgba(255,255,255,0.8)",
        bordercolor="black",
        font=attr(size=18, family="serif")
    ),
    xaxis=attr(title="", showline=true, linecolor="black", linewidth=2, mirror=true,
               titlefont=attr(size=18, family="serif"), tickfont=attr(size=16, family="serif"),
               showgrid=true, gridcolor="rgba(120,120,120,0.25)"),
    yaxis=attr(title="", showline=true, linecolor="black", linewidth=2, mirror=true,
               titlefont=attr(size=18, family="serif"), tickfont=attr(size=14, family="serif"),
               showgrid=true, gridcolor="rgba(120,120,120,0.25)"),
    font=attr(family="serif", size=18)
)
plt_tail = Plot(tail_traces, layout_tail)
PlotlyJS.savefig(plt_tail, "Figures/STATCOM_Isrc_max_LDC_2_tail_plotly.pdf")

# Overlay tail_traces as an inset on the top right corner of plt_main
latex_font = attr(family="serif", size=18, color="black")

# derive n from first main trace
n = length(traces[1][:y])
x_rect_percent = 100 * (min(300, n) / n)   # rectangle x-span in percent on main axis

# compute zoom y extents from all tail traces but ignore NaNs
# Force the zoom extents to match the inset figure exactly: 260..390
y_zoom_min = 260.0
y_zoom_max = 390.0

# compute overall y-range for main plot (to map data y to paper coordinates)
all_main_y = reduce(vcat, [t[:y] for t in traces])
y_main_min = minimum(all_main_y)
y_main_max = maximum(all_main_y)
y_span = max(1e-9, y_main_max - y_main_min)

# compute rectangle center in paper coordinates
x_rect_center_data = x_rect_percent / 2.0            # percent units (0..100)
xcenter_paper = x_rect_center_data / 100.0           # domain [0..1]
y_rect_center_data = (y_zoom_min + y_zoom_max) / 2.0
ycenter_paper = (y_rect_center_data - y_main_min) / y_span
ycenter_paper = clamp(ycenter_paper, 0.0, 1.0)

# inset domains
inset_x0, inset_x1 = 0.65, 0.98
inset_y0, inset_y1 = 0.55, 0.98
inset_yc = (inset_y0 + inset_y1) / 2.0

main_xaxis = attr(
    title="Load percentile (%)",
    domain=[0.0, 1.0],
    anchor="y",
    titlefont=latex_font,
    tickfont=latex_font,
    showline=true, linecolor="black", linewidth=2, mirror=true,
    showgrid=true, gridcolor="rgba(120,120,120,0.25)", gridwidth=1
)
main_yaxis = attr(
    title="Maximum Feeder F2 Substation<br>Phase Current (A)",
    domain=[0.0, 1.0],
    anchor="x",
    titlefont=latex_font,
    tickfont=attr(size=14, family="serif"),
    showline=true, linecolor="black", linewidth=2, mirror=true,
    showgrid=true, gridcolor="rgba(120,120,120,0.25)", gridwidth=1
)

# inset x now spans -0.03..1.0 per request, but xticks at 0:0.2:1
inset_xaxis = attr(
    domain=[inset_x0, inset_x1],
    anchor="y2",
    title="",
    range=[-0.03, 1.0],
    tickvals=collect(0.0:0.2:1.0),
    showgrid=true,
    showticklabels=true,
    showlegend=false,
    titlefont=latex_font,
    tickfont=latex_font,
    showline=true, linecolor="black", linewidth=2, mirror=true,
    gridcolor="rgba(120,120,120,0.25)"
)
# set inset y-range to 260..390 as requested
inset_yaxis = attr(
    domain=[inset_y0, inset_y1],
    anchor="x2",
    title="",
    range=[260.0, 390.0],
    showgrid=true,
    showticklabels=true,
    titlefont=latex_font,
    tickfont=attr(size=14, family="serif"),
    showline=true, linecolor="black", linewidth=2, mirror=true,
    gridcolor="rgba(120,120,120,0.25)"
)

for t in tail_traces
    t[:xaxis] = "x2"
    t[:yaxis] = "y2"
    t[:showlegend] = false
end

all_traces = vcat(traces, tail_traces)

# Add rectangle shape on main plot and a connecting dashed line to inset
shapes = [
    Dict(
        :type => "rect",
        :xref => "x", :yref => "y",
        :x0 => 0.0, :x1 => x_rect_percent,
        :y0 => 260.0, :y1 => 390.0,  # Fixed to match inset
        :line => attr(color="black", width=1),
        :fillcolor => "rgba(0,0,0,0)"
    ),
    Dict(
        :type => "line",
        :xref => "paper", :yref => "paper",
        :x0 => clamp(xcenter_paper, 0.02, 0.98), :y0 => clamp(ycenter_paper, 0.02, 0.98),
        :x1 => inset_x0 - 0.005, :y1 => inset_yc,
        :line => attr(color="black", width=1, dash="dash")
    )
]

layout_overlay = Layout(
    xaxis=main_xaxis,
    yaxis=main_yaxis,
    xaxis2=inset_xaxis,
    yaxis2=inset_yaxis,
    width=900,
    height=580,
    plot_bgcolor="white",
    paper_bgcolor="white",
    shapes=shapes,
    showlegend=true,
    legend=attr(
        orientation="h",
        x=0.0,
        y=1.26,
        xanchor="left",
        yanchor="top",
        bgcolor="rgba(255,255,255,0.8)",
        bordercolor="black",
        font=latex_font
    ),
    font=latex_font
)

plt_overlay = Plot(all_traces, layout_overlay)
PlotlyJS.savefig(plt_overlay, "Figures/STATCOM_Isrc_max_LDC_2_overlay_plotly.pdf")


