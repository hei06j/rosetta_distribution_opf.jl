using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP  # bl/array_nl
import LinearAlgebra: diag, diagm
using LaTeXStrings

const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")

objective = "IUF2_inv"

##

function parse_data(data_path)
    data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!, PMD.transform_loops!])
    # RPMD.pv1_correction!(data_eng)
    data_eng["settings"]["sbase_default"] = 1
    data_eng["voltage_source"]["source"]["rs"] *= 0
    data_eng["voltage_source"]["source"]["xs"] *= 0
    data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)

    for (i, bus) in data_math["bus"]
        bus["vmin"] = [0.9 * ones(3) ; 0 ]
        bus["vmax"] = [1.1 * ones(3) ; Inf]
    end

    gen_id = 1
    gen = data_math["gen"]["$gen_id"]
    smax = 40
    pmax = 3
    gen["pmax"] = pmax/3 * ones(3)
    gen["pmin"] = zeros(3)
    gen["qmax"] = sqrt.(smax^2 - pmax^2)/3 * ones(3)
    gen["qmin"] = -gen["qmax"]

    data_math["gen"]["1"]["cost"] = [10 0]
    data_math["gen"]["2"]["cost"] = [10 0]

    return data_math
end

function get_solutions!(model, results)
    JuMP.optimize!(model)
    @assert(JuMP.termination_status(model) == LOCALLY_SOLVED)
    cost = JuMP.objective_value(model)

    gen_id = [parse(Int,i) for (i,gen) in data_math["gen"] if occursin("pv", gen["name"])][1]
    gen_bus_id = data_math["gen"]["$gen_id"]["gen_bus"]
    
    vr_vals = value.(vr)
    vi_vals = value.(vi)
    v1 = vr_vals[:,gen_bus_id] + im * vi_vals[:,gen_bus_id]
    vm1 = abs.(v1)
    # va1 = angle.(v1).*180/pi
    # va1_pn = va1[1:3] .- va1[4]
    # va1_pn_pp = Array(va1_pn[[1,2,3]]) .- Array(va1_pn[[2,3,1]])
    # va1_pp = Array(va1[[1,2,3]]) .- Array(va1[[2,3,1]])
    v1_012 = T * v1[1:3]
    vm1_012 = abs.(v1_012)
    
    crg_values = JuMP.value.(crg_bus)
    cig_values = JuMP.value.(cig_bus)
    cg1 = crg_values[:,gen_id] + im * cig_values[:,gen_id]
    cgm1 = abs.(cg1)
    cg1_012 = T * cg1[1:3]
    cgm1_012 = abs.(cg1_012)
    # cga1 = angle.(cg1).*180/pi
    
    crd_values = JuMP.value.(crd_bus)
    cid_values = JuMP.value.(cid_bus)
    cd1 = sum(Array(crd_values .+ im * cid_values), dims=2)
    # cdm1 = abs.(cd1)
    cd1_012 = T * cd1[1:3]
    cdm1_012 = abs.(cd1_012)
    
    branch_lij = (1, 2, 1)
    cr1 = value.(cr_bus)[:,branch_lij]
    ci1 = value.(ci_bus)[:,branch_lij]
    c1 = cr1 .+ im * ci1
    # cm1 = abs.(c1)
    c1_012 = T * c1[1:3]
    cm1_012 = abs.(c1_012)

    # branch_idx = (1, 2, 1)
    # cr_vals = value.(cr)
    # ci_vals = value.(ci)
    # c = cr_vals[:,branch_idx] .+ im *  ci_vals[:,branch_idx]
    # c012 = T * Array(c[1:3])
    # c2 = abs.(c012)[3]
    
    results[objective] = Dict()
    results[objective]["v_gen"] = v1
    results[objective]["v_gen_012"] = v1_012
    results[objective]["c_gen"] = cg1
    results[objective]["c_gen_012"] = cg1_012
    results[objective]["c_load"] = cd1
    results[objective]["c_load_012"] = cd1_012
    results[objective]["c_branch"] = c1
    results[objective]["c_branch_012"] = c1_012
    results[objective]["sg_gen"] = Array(value.(pg)[:,1] .+ im * value.(qg)[:,1])
    results[objective]["sg_source"] = Array(value.(pg)[:,2] .+ im * value.(qg)[:,2])
    
    return results
end


## 3-leg inverters
data_path = "./data/inverter_3w_wye_unbalanced_loads.dss"
data_math = parse_data(data_path)
# gen["connections"] = gen["connections"][1:3]
# data_math["bus"]["1"]["terminals"] = data_math["bus"]["1"]["terminals"][1:4]
# data_math["bus"]["1"]["grounded"] = data_math["bus"]["1"]["grounded"][1:4]
include("./inverter_loss_branch.jl")
add_inverter_losses(data_math, gen_id, three_wire=true)
ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]

model = JuMP.Model(Ipopt.Optimizer)
include("./variables.jl")
include("./constraints.jl")
include("./objectives.jl")
# include("./VV_VW_controls.jl")

results_GFL_3w = Dict()
get_solutions!(model, results_GFL_3w)


## 4-leg inverters
data_path = "./data/inverter_4w_wye_unbalanced_loads.dss"
data_math = parse_data(data_path)
include("./inverter_loss_branch.jl")
add_inverter_losses(data_math, gen_id)
ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]

model = JuMP.Model(Ipopt.Optimizer)
include("./variables.jl")
include("./constraints.jl")
include("./objectives.jl")
# include("./VV_VW_controls.jl")

results_GFL_4w = Dict()
get_solutions!(model, results_GFL_4w)


###
key = [key for (key,value) in results_GFL_3w][1]
round.(abs.(results_GFL_3w[key]["v_gen"]), digits=4)
round.(results_GFL_3w[key]["sg_gen"], digits=3)
round.(results_GFL_3w[key]["c_gen"], digits=3)

round.(abs.(results_GFL_4w[key]["v_gen"]), digits=4)
round.(results_GFL_4w[key]["sg_gen"], digits=3)
round.(results_GFL_4w[key]["c_gen"], digits=3)


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

using Plots
mkpath("./Figures")

###
key = [key for (key,value) in results_GFL_3w][1]
Imax_3w = maximum([abs.(results_GFL_3w[key]["c_branch"])  abs.(results_GFL_3w[key]["c_load"]) abs.(results_GFL_3w[key]["c_gen"])])
Imax_4w = maximum([abs.(results_GFL_4w[key]["c_branch"])  abs.(results_GFL_4w[key]["c_load"]) abs.(results_GFL_4w[key]["c_gen"])])
I0_3w = round(abs(results_GFL_3w[key]["c_branch_012"][1]), digits=2)
I2_3w = round(abs(results_GFL_3w[key]["c_branch_012"][3]), digits=2)
Imax = maximum([Imax_3w, Imax_4w])
GFL_3w_vuf_c = plot_phasors(results_GFL_3w[key]["c_branch"], Imax; I2=I2_3w, I0=I0_3w)
Plots.savefig(GFL_3w_vuf_c, "./Figures/GFL_3w_vuf_c.pdf")
GFL_3w_vuf_cg = plot_phasors(results_GFL_3w[key]["c_gen"], Imax)
Plots.savefig(GFL_3w_vuf_cg, "./Figures/GFL_3w_vuf_cg.pdf")
GFL_3w_vuf_cd = plot_phasors(results_GFL_3w[key]["c_load"], Imax)
Plots.savefig(GFL_3w_vuf_cd, "./Figures/GFL_3w_vuf_cd.pdf")
GFL_3w_vuf = Plots.plot(GFL_3w_vuf_c, GFL_3w_vuf_cg, GFL_3w_vuf_cd, layout=(1,3), size=(2200,700))
Plots.savefig(GFL_3w_vuf, "./Figures/GFL_3w_vuf.pdf")
Plots.savefig(GFL_3w_vuf, "./Figures/GFL_3w_vuf.png")


key = [key for (key,value) in results_GFL_4w][1]
# Imax = maximum([abs.(results_GFL_4w[key]["c_branch"])  abs.(results_GFL_4w[key]["c_load"]) abs.(results_GFL_4w[key]["c_gen"])])
I0_4w = round(abs(results_GFL_4w[key]["c_branch_012"][1]), digits=2)
I2_4w = round(abs(results_GFL_4w[key]["c_branch_012"][3]), digits=2)
GFL_4w_vuf_c = plot_phasors(results_GFL_4w[key]["c_branch"], Imax; I2=I2_4w, I0=I0_4w)
Plots.savefig(GFL_4w_vuf_c, "./Figures/GFL_4w_vuf_c.pdf")
GFL_4w_vuf_cg = plot_phasors(results_GFL_4w[key]["c_gen"], Imax)
Plots.savefig(GFL_4w_vuf_cg, "./Figures/GFL_4w_vuf_cg.pdf")
GFL_4w_vuf_cd = plot_phasors(results_GFL_4w[key]["c_load"], Imax, labeled=true)
Plots.savefig(GFL_4w_vuf_cd, "./Figures/GFL_4w_vuf_cd.pdf")
GFL_4w_vuf = Plots.plot(GFL_4w_vuf_c, GFL_4w_vuf_cg, GFL_4w_vuf_cd, layout=(1,3), size=(2200,700))
Plots.savefig(GFL_4w_vuf, "./Figures/GFL_4w_vuf.pdf")
Plots.savefig(GFL_4w_vuf, "./Figures/GFL_4w_vuf.png")

GFL_4w_3w_vuf = Plots.plot(GFL_3w_vuf, GFL_4w_vuf, layout=(2,1), size=(2100,1400))
# Plots.savefig(GFL_4w_3w_vuf, "./Figures/GFL_4w_3w_vuf.pdf")


##
[results_GFL_4w["loss"]["sg_gen"] results_GFL_3w["loss"]["sg_gen"]]'
abs.([results_GFL_4w["loss"]["v_gen"] results_GFL_3w["loss"]["v_gen"]])'
angle.([results_GFL_4w["loss"]["v_gen"] results_GFL_3w["loss"]["v_gen"]])' .* 180/pi

[results_GFL_4w["loss"]["c_gen"] results_GFL_3w["loss"]["c_gen"]]'
