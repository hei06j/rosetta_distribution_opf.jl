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
highs = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
juniper_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt_solver, "mip_solver" => highs)

objective = "IUF2_inv"
# objective = "loss"

"""
- assymetrical alpha's 
- start with idealised version
- change the testcase
- timeseries analysis
"""

data_path = "./data/inverter_4w_wye_unbalanced_loads.dss"
data_path = "./data/case5_gen_3ph_wye.dss"


data_eng = PMD.parse_file(data_path)
# data_eng["settings"]["sbase_default"] = 1.0
# data_eng["line"]["line1"]["cm_ub"] = [6:-1:3...]
data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)
PMD.add_start_vrvi!(data_math)
PMD.solve_mc_opf(data_math, PMD.IVRENPowerModel, ipopt_solver)

pm_ivr  = PMD.instantiate_mc_model(data_eng, PMD.IVRUPowerModel, PMD.build_mc_opf)
sol_ivr = PMD.optimize_model!(pm_ivr, optimizer=ipopt_solver)


##
function parse_data(data_path)
    data_eng = PMD.parse_file(data_path)
    # data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
    # data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!, PMD.transform_loops!])
    # data_eng["settings"]["sbase_default"] = 1
    # data_eng["voltage_source"]["source"]["rs"] *= 0
    # data_eng["voltage_source"]["source"]["xs"] *= 0
    data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)

    for (i, bus) in data_math["bus"]
        bus["vmin"] = [0.9 * ones(3) ; 0 ]
        bus["vmax"] = [1.1 * ones(3) ; Inf]
    end

    # gen_id = 1
    # gen = data_math["gen"]["$gen_id"]
    # smax = 40
    # pmax = 30
    # gen["pmax"] = pmax/3 * ones(3)
    # gen["pmin"] = zeros(3)
    # gen["qmax"] = sqrt.(smax^2 - pmax^2)/3 * ones(3)
    # gen["qmin"] = -gen["qmax"]

    data_math["gen"]["1"]["cost"] = [10 0]
    data_math["gen"]["2"]["cost"] = [10 0]

    return data_math
end



function get_solutions!(model, results, objective)
    JuMP.optimize!(model)
    @assert(JuMP.termination_status(model) == LOCALLY_SOLVED)
    cost = JuMP.objective_value(model)

    pv_gen_id = [parse(Int,i) for (i,gen) in data_math["gen"] if occursin("pv", gen["name"])][1]
    gen_bus_id = data_math["gen"]["$pv_gen_id"]["gen_bus"]
    
    vr_vals = value.(vr)
    vi_vals = value.(vi)
    v_pv = vr_vals[:,gen_bus_id] + im * vi_vals[:,gen_bus_id]
    v_pv_012 = T * v_pv[1:3]
    vm_pv = abs.(v_pv)
    vm_pv_012 = abs.(v_pv_012)    # va1 = angle.(v_pv).*180/pi
    # va1_pn = va1[1:3] .- va1[4]
    # va1_pn_pp = Array(va1_pn[[1,2,3]]) .- Array(va1_pn[[2,3,1]])
    # va1_pp = Array(va1[[1,2,3]]) .- Array(va1[[2,3,1]])

    
    crg_values = JuMP.value.(crg_bus)
    cig_values = JuMP.value.(cig_bus)
    cg_pv = crg_values[:,pv_gen_id] + im * cig_values[:,pv_gen_id]
    cg_pv_012 = T * cg_pv[1:3]
    cgm_pv = abs.(cg_pv)
    cgm_pv_012 = abs.(cg_pv_012)
    # cga1 = angle.(cg_pv).*180/pi

    cg_src = crg_values[:,2] + im * cig_values[:,2]
    cg_src_012 = T * cg_src[1:3]
    cgm_src = abs.(cg_src)
    cgm_src_012 = abs.(cg_src_012)

    crd_values = JuMP.value.(crd_bus)
    cid_values = JuMP.value.(cid_bus)
    cd = sum(Array(crd_values .+ im * cid_values), dims=2)
    cd_012 = T * cd[1:3]
    # cdm = abs.(cd)
    # cdm_012 = abs.(cd_012)
    
    # pv_branch = (1, 2, 1)
    _, _, pv_branch = RPMD.get_pv_bus_branch(ref)
    c_pv = Array{Float64}(undef, 4, 0)
    c_pv_012 = Array{Float64}(undef, 3, 0)
    for branch in pv_branch
        cr_pv = value.(cr_bus)[:,branch]
        ci_pv = value.(ci_bus)[:,branch]
        c = cr_pv .+ im * ci_pv
        c_pv = [c_pv, c]
        c_pv_012 = [c_pv_012, T*c[1:3]]
    end


    # ref_branch = (2, 3, 1)
    ref_gen, ref_bus, ref_arc, ref_branch = RPMD.get_ref_bus_branch(ref)
    cr_ref = value.(cr_bus)[:,ref_arc]
    ci_ref = value.(ci_bus)[:,ref_arc]
    c_ref = cr_ref .+ im * ci_ref
    c_ref_012 = T * c_ref[1:3]

    
    results[objective] = Dict()
    results[objective]["v_inv"] = v_pv
    results[objective]["v_inv_012"] = v_pv_012
    results[objective]["c_inv"] = cg_pv
    results[objective]["c_inv_012"] = cg_pv_012
    results[objective]["c_source"] = cg_src
    results[objective]["c_source_012"] = cg_src_012
    results[objective]["c_load"] = cd
    results[objective]["c_load_012"] = cd_012
    results[objective]["c_branch1"] = c_pv
    results[objective]["c_branch1_012"] = c_pv_012
    results[objective]["c_branch2"] = c_ref
    results[objective]["c_branch2_012"] = c_ref_012
    results[objective]["sg_inv"] = Array(value.(pg)[:,1] .+ im * value.(qg)[:,1])
    results[objective]["sg_source"] = Array(value.(pg)[:,2] .+ im * value.(qg)[:,2])
    
    return results
end


function plot_phasors(phasor, Imax; labeled=false, I2=[], I0=[])
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


function add_inverter_model!(data_math)
    include("./core/inverter_loss_branch.jl")
    pv_gen_ids = [i for (i, gen) in data_math["gen"] if !occursin("source", gen["name"])]
    for gen_id in pv_gen_ids
        add_inverter_losses!(data_math, gen_id)
    end
end

## 4-leg inverters with multiplexing
multileg = true
multiplexing = true
multiplexing_binary = true

data_math = parse_data(data_path)
add_inverter_model!(data_math)
ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]

m_legs = 8
if multiplexing_binary
    model = JuMP.Model(juniper_solver)
else
    model = JuMP.Model(ipopt_solver)
end
include("./core/variables.jl")
include("./core/constraints.jl")
include("./core/objectives.jl")


results_GFL_4w_mx = Dict()
PMD.add_start_vrvi!(data_math)
get_solutions!(model, results_GFL_4w_mx, objective)

##
@show Bg_mx = value.(bg["1"])
[(results_GFL_4w_mx["IUF2_inv"]["c_load"]) (results_GFL_4w_mx["IUF2_inv"]["c_branch1"]) (results_GFL_4w_mx["IUF2_inv"]["c_branch2"])]
[abs.(results_GFL_4w_mx["IUF2_inv"]["c_inv"]) abs.(results_GFL_4w_mx["IUF2_inv"]["c_source"]) abs.(results_GFL_4w_mx["IUF2_inv"]["c_load"]) abs.(results_GFL_4w_mx["IUF2_inv"]["c_branch1"]) abs.(results_GFL_4w_mx["IUF2_inv"]["c_branch2"])]
[abs.(results_GFL_4w_mx["IUF2_inv"]["c_load_012"]) abs.(results_GFL_4w_mx["IUF2_inv"]["c_branch1_012"]) abs.(results_GFL_4w_mx["IUF2_inv"]["c_branch2_012"])]
[abs.(results_GFL_4w_mx["IUF2_inv"]["c_load"]) abs.(results_GFL_4w_mx["IUF2_inv"]["c_branch1"]).+abs.(results_GFL_4w_mx["IUF2_inv"]["c_branch2"])]


## 4-leg inverters without multiplexing
multileg = true
multiplexing = false
multiplexing_binary = false

data_math = parse_data(data_path)
# add_inverter_model!(data_math)
ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]


# m_legs = 8
model = JuMP.Model(ipopt_solver)
include("./core/variables.jl")
include("./core/constraints.jl")
include("./core/objectives.jl")

results_GFL_4w = Dict()
PMD.add_start_vrvi!(data_math)
get_solutions!(model, results_GFL_4w, objective)

@show Bg_mx = value.(bg["1"])
[(results_GFL_4w["IUF2_inv"]["c_load"]) (results_GFL_4w["IUF2_inv"]["c_branch1"]) (results_GFL_4w["IUF2_inv"]["c_branch2"])]
[abs.(results_GFL_4w["IUF2_inv"]["c_inv"]) abs.(results_GFL_4w["IUF2_inv"]["c_source"]) abs.(results_GFL_4w["IUF2_inv"]["c_load"]) abs.(results_GFL_4w["IUF2_inv"]["c_branch1"]) abs.(results_GFL_4w["IUF2_inv"]["c_branch2"])]
[abs.(results_GFL_4w["IUF2_inv"]["c_load_012"]) abs.(results_GFL_4w["IUF2_inv"]["c_branch1_012"]) abs.(results_GFL_4w["IUF2_inv"]["c_branch2_012"])]
[abs.(results_GFL_4w["IUF2_inv"]["c_load"]) abs.(results_GFL_4w["IUF2_inv"]["c_branch1"]).+abs.(results_GFL_4w["IUF2_inv"]["c_branch2"])]



##
key = [key for (key,value) in results_GFL_4w][1]
round.(abs.(results_GFL_4w[key]["v_inv"]), digits=4)
round.(results_GFL_4w[key]["sg_inv"], digits=3)
round.(results_GFL_4w[key]["c_inv"], digits=3)

round.(abs.(results_GFL_4w_mx[key]["v_inv"]), digits=4)
round.(results_GFL_4w_mx[key]["sg_inv"], digits=3)
round.(results_GFL_4w_mx[key]["c_inv"], digits=3)


using Plots
mkpath("./Figures")

##
key = [key for (key,value) in results_GFL_4w][1]
Imax = maximum(abs.(results_GFL_4w[key]["c_load"]))
I0_4w = round(abs(results_GFL_4w[key]["c_branch1_012"][1]), digits=2)
I2_4w = round(abs(results_GFL_4w[key]["c_branch1_012"][3]), digits=2)
GFL_4w_vuf_c = plot_phasors(results_GFL_4w[key]["c_branch1"], Imax; I2=I2_4w, I0=I0_4w)
Plots.savefig(GFL_4w_vuf_c, "./Figures/GFL_4w_vuf_c.pdf")
GFL_4w_vuf_cg = plot_phasors(results_GFL_4w[key]["c_inv"], Imax)
Plots.savefig(GFL_4w_vuf_cg, "./Figures/GFL_4w_vuf_cg.pdf")
GFL_4w_vuf_cd = plot_phasors(results_GFL_4w[key]["c_load"], Imax)
Plots.savefig(GFL_4w_vuf_cd, "./Figures/GFL_4w_vuf_cd.pdf")
GFL_4w_vuf = Plots.plot(GFL_4w_vuf_c, GFL_4w_vuf_cg, GFL_4w_vuf_cd, layout=(1,3), size=(2200,700))
Plots.savefig(GFL_4w_vuf, "./Figures/GFL_4w_vuf.pdf")
Plots.savefig(GFL_4w_vuf, "./Figures/GFL_4w_vuf.png")


key = [key for (key,value) in results_GFL_4w_mx][1]
# Imax = maximum([abs.(results_GFL_4w_mx[key]["c_branch"])  abs.(results_GFL_4w_mx[key]["c_load"]) abs.(results_GFL_4w_mx[key]["c_inv"])])
Imax = maximum(abs.(results_GFL_4w_mx[key]["c_load"]))
I0_4w_mx = round(abs(results_GFL_4w_mx[key]["c_branch1_012"][1]), digits=2)
I2_4w_mx = round(abs(results_GFL_4w_mx[key]["c_branch1_012"][3]), digits=2)
GFL_4w_mx_vuf_c = plot_phasors(results_GFL_4w_mx[key]["c_branch1"], Imax; I2=I2_4w_mx, I0=I0_4w_mx)
Plots.savefig(GFL_4w_mx_vuf_c, "./Figures/GFL_4w_mx_vuf_c.pdf")
GFL_4w_mx_vuf_cg = plot_phasors(results_GFL_4w_mx[key]["c_inv"], Imax)
Plots.savefig(GFL_4w_mx_vuf_cg, "./Figures/GFL_4w_mx_vuf_cg.pdf")
GFL_4w_mx_vuf_cd = plot_phasors(results_GFL_4w_mx[key]["c_load"], Imax, labeled=true)
Plots.savefig(GFL_4w_mx_vuf_cd, "./Figures/GFL_4w_mx_vuf_cd.pdf")
GFL_4w_mx_vuf = Plots.plot(GFL_4w_mx_vuf_c, GFL_4w_mx_vuf_cg, GFL_4w_mx_vuf_cd, layout=(1,3), size=(2200,700))
Plots.savefig(GFL_4w_mx_vuf, "./Figures/GFL_4w_mx_vuf.pdf")
Plots.savefig(GFL_4w_mx_vuf, "./Figures/GFL_4w_mx_vuf.png")

GFL_4w_3w_vuf = Plots.plot(GFL_4w_vuf, GFL_4w_mx_vuf, layout=(2,1), size=(2200,1400))
Plots.savefig(GFL_4w_3w_vuf, "./Figures/GFL_4w_4w_mx_vuf.pdf")


##
# [results_GFL_4w_mx["loss"]["sg_inv"] results_GFL_4w["loss"]["sg_inv"]]'
# abs.([results_GFL_4w_mx["loss"]["v_inv"] results_GFL_4w["loss"]["v_inv"]])'
# angle.([results_GFL_4w_mx["loss"]["v_inv"] results_GFL_4w["loss"]["v_inv"]])' .* 180/pi

# [results_GFL_4w_mx["loss"]["c_inv"] results_GFL_4w["loss"]["c_inv"]]'
