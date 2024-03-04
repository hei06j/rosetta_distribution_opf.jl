using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP  # bl/array_nl
import LinearAlgebra: diag, diagm
using Plots
using LaTeXStrings

const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels


##
data_path = "./data/inverter_4w_wye_unbalanced_loads.dss"

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes", "max_iter"=>10000)
data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!])
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
gen["pmax"] = 35/3 * ones(3)
gen["pmin"] = zeros(3)
# gen["pmax"] = 23/3 * ones(3)
# gen["pmin"] = 23/3 * ones(3)
# gen["pg"] = 23/3 * ones(3)
gen["qmax"] = sqrt.(40^2 - 35^2)/3 * ones(3)
gen["qmin"] = -gen["qmax"]


include("./core/inverter_loss_branch.jl")
add_inverter_losses(data_math, gen_id)

data_math["gen"]["1"]["cost"] = [10 0]
data_math["gen"]["2"]["cost"] = [1000 0]

ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]

model =[]
vv_vw_plots = Dict()
control = []
objective = "cost"

##
for control in ["Volt_Watt_pn", "Volt_var_pn", "Volt_Watt_pp", "Volt_var_pp", "Volt_Watt_av", "Volt_var_av"]
    
    model = JuMP.Model(ipopt_solver)
    include("./core/variables.jl")
    include("./core/constraints.jl")
    include("./core/objectives.jl")

    gen_id = 1
    gen = data_math["gen"]["$gen_id"]
    bus_id = gen["gen_bus"]
    vmin = 0.9; vmax = 1.1;
    if 4 in gen["connections"]
        phases = gen["connections"][1:end-1]
        n = gen["connections"][end]
    else
        phases = gen["connections"]
    end

    terminals = data_math["bus"]["$bus_id"]["terminals"]
    vm = JuMP.@variable(model, [t in terminals, i=bus_id], base_name="vm", lower_bound=0)
    JuMP.@constraint(model, [t in terminals], vm[t,bus_id]^2 == vr[t,bus_id]^2 + vi[t,bus_id]^2)

    vmpp = JuMP.@variable(model, [t in 1:3, i=bus_id], base_name="vmpp", lower_bound=0)
    JuMP.@constraint(model, vmpp[1,bus_id]^2 == (vr[1,bus_id]-vr[2,bus_id])^2 + (vi[1,bus_id]-vi[2,bus_id])^2)
    JuMP.@constraint(model, vmpp[2,bus_id]^2 == (vr[2,bus_id]-vr[3,bus_id])^2 + (vi[2,bus_id]-vi[3,bus_id])^2)
    JuMP.@constraint(model, vmpp[3,bus_id]^2 == (vr[3,bus_id]-vr[1,bus_id])^2 + (vi[3,bus_id]-vi[1,bus_id])^2)

    vmpn = JuMP.@variable(model, [t in 1:3, i=bus_id], base_name="vmpn", lower_bound=0)
    JuMP.@constraint(model, [p in phases], vmpn[p,bus_id].^2 == (vr[p,bus_id].-vr[4,bus_id]).^2 + (vi[p,bus_id].-vi[4,bus_id]).^2)
    
    vv_curve(v, qmax) = v <= 0.9 ? qmax : (0.9<=v<=0.95 ? -qmax/0.05*(v-0.95) : (0.95<=v<=1.05 ? 0.0 : (1.05<=v<=1.1 ? -qmax/0.05*(v-1.05) : -qmax)))
    vw_curve(v, pmax) = v <= 0.9 ? pmax : (v >= 1.1 ? 0.0 : -pmax/0.2*(v-1.1))
    JuMP.@operator(model, vv, 2, vv_curve)
    JuMP.@operator(model, vw, 2, vw_curve)

    ### control = Volt_Watt_pn, Volt_var_pn, Volt_Watt_pp, Volt_var_pp, Volt_Watt_av, Volt_var_av
    if control == "Volt_var_pn"
        JuMP.@constraint(model, [i in phases], qg[i,gen_id] == vv(vmpn[i,bus_id], gen["qmax"][i]) )
        # JuMP.@constraint(model, [i in phases], qg[i,gen_id] == vv(vm[i,bus_id]-vm[n,bus_id], gen["qmax"][i]) )

    elseif control == "Volt_Watt_pn"
        # JuMP.@constraint(model, [i in phases], pg[i,gen_id] == vw(vm[i,bus_id], gen["pmax"][i]) )
        # JuMP.@constraint(model, [i in phases], pg[i,gen_id] == vw(vm[i,bus_id]-vm[n,bus_id], gen["pmax"][i]) )
        JuMP.@constraint(model, [i in phases], pg[i,gen_id] == vw(vmpn[i,bus_id], gen["pmax"][i]) )

    elseif control == "Volt_var_pp"
        JuMP.@constraint(model, qg[1,gen_id] - qg[2,gen_id] == vv(vmpp[1,bus_id]/sqrt(3), gen["qmax"][1]) )
        JuMP.@constraint(model, qg[2,gen_id] - qg[3,gen_id] == vv(vmpp[2,bus_id]/sqrt(3), gen["qmax"][2]) )
        JuMP.@constraint(model, qg[3,gen_id] - qg[1,gen_id] == vv(vmpp[3,bus_id]/sqrt(3), gen["qmax"][3]) )

    elseif control == "Volt_Watt_pp"
        JuMP.@constraint(model, pg[1,gen_id] - pg[2,gen_id] == vw(vmpp[1,bus_id]/sqrt(3), gen["pmax"][1]) )
        JuMP.@constraint(model, pg[2,gen_id] - pg[3,gen_id] == vw(vmpp[2,bus_id]/sqrt(3), gen["pmax"][2]) )
        JuMP.@constraint(model, pg[3,gen_id] - pg[1,gen_id] == vw(vmpp[3,bus_id]/sqrt(3), gen["pmax"][3]) )

    elseif control == "Volt_var_av"
        JuMP.@constraint(model, qg[1, gen_id] == vv(sum(vm[i,bus_id] for i in ref[:bus][1]["terminals"][1:3])/3, gen["qmax"][1]) )
        JuMP.@constraint(model, qg[2, gen_id] == qg[1, gen_id])
        JuMP.@constraint(model, qg[3, gen_id] == qg[1, gen_id])

    elseif control == "Volt_Watt_av"
        JuMP.@constraint(model, pg[1, gen_id] == vw(sum(vm[i,bus_id] for i in ref[:bus][1]["terminals"][1:3])/3, gen["pmax"][1]) )
        JuMP.@constraint(model, pg[2, gen_id] == pg[1, gen_id])
        JuMP.@constraint(model, pg[3, gen_id] == pg[1, gen_id])
    end

    ###
    JuMP.optimize!(model)
    cost = JuMP.objective_value(model)

    ##
    gen_id = [parse(Int,i) for (i,gen) in data_math["gen"] if occursin("pv", gen["name"])][1]
    bus_id = data_math["gen"]["$gen_id"]["gen_bus"]
    pg_values = JuMP.value.(pg)
    qg_values = JuMP.value.(qg)

    ##
    plt = Plots.plot()
    vm_range = 0.85:0.01:1.15
    if control in ["Volt_var_pn", "Volt_var_pp", "Volt_var_av"]
        if control == "Volt_var_pn"
            vmpn_vals = Array(value.(vmpn))
            # plt_Volt_var_pn = plot(vm_range, vv_curve.(vm_range,data_math["gen"]["1"]["qmax"][1]), label=false, guidefontsize=10)
            plt = plot(vm_range, vv_curve.(vm_range,data_math["gen"]["1"]["qmax"][1]), label=false, guidefontsize=10)
            plot!([vmpn_vals[1]], [(qg_values[1,gen_id])], seriestype=:scatter, label=L"(q_{an}, U_{an})", xticks=(0.85:0.05:1.15, 0.85:0.05:1.15), xrotation = 45)
            plot!([vmpn_vals[2]], [(qg_values[2,gen_id])], seriestype=:scatter, label=L"(q_{bn}, U_{bn})")
            plot!([vmpn_vals[3]], [(qg_values[3,gen_id])], seriestype=:scatter, label=L"(q_{cn}, U_{cn})")
            # xlabel!("Voltage Magnitudes (V pu)")
            ylabel!("Reactive Power (kvar)")
            
        elseif control == "Volt_var_pp"
            vmpp_vals = Array(value.(vmpp))
            # plt_Volt_var_pp = plot(vm_range, vv_curve.(vm_range,data_math["gen"]["1"]["qmax"][1]), label=false, guidefontsize=10)
            plt = plot(vm_range, vv_curve.(vm_range,data_math["gen"]["1"]["qmax"][1]), label=false, guidefontsize=10)
            plot!([vmpp_vals[1]/sqrt(3)], [qg_values[1,gen_id]-qg_values[2,gen_id]], seriestype=:scatter, label=L"(q_{ab}, U_{ab})", xticks=(0.85:0.05:1.15, 0.85:0.05:1.15), xrotation = 45)
            plot!([vmpp_vals[2]/sqrt(3)], [qg_values[2,gen_id]-qg_values[3,gen_id]], seriestype=:scatter, label=L"(q_{bc}, U_{bc})")
            plot!([vmpp_vals[3]/sqrt(3)], [qg_values[3,gen_id]-qg_values[1,gen_id]], seriestype=:scatter, label=L"(q_{ca}, U_{ca})")
            # xlabel!("Voltage Magnitudes (V pu)")
            # ylabel!("Reactive Power (kvar)")

        elseif control == "Volt_var_av"
            vm_vals = Array(value.(vm))
            # plt_Volt_var_av = plot(vm_range, vv_curve.(vm_range,data_math["gen"]["1"]["qmax"][1]), label=false, guidefontsize=10)
            plt = plot(vm_range, vv_curve.(vm_range,data_math["gen"]["1"]["qmax"][1]), label=false, guidefontsize=10)
            plot!([sum(vm_vals[phases])/3], [qg_values[1,gen_id]], seriestype=:scatter, label=L"(q_p, (U_a+U_b+U_c)/3)", legend=:bottomleft, xticks=(0.85:0.05:1.15, 0.85:0.05:1.15), xrotation = 45)
            # xlabel!("Voltage Magnitudes (V pu)")
            # ylabel!("Reactive Power (kvar)")
            
        end

    elseif control in ["Volt_Watt_pn", "Volt_Watt_pp", "Volt_Watt_av"]
        if control == "Volt_Watt_pn"
            vmpn_vals = Array(value.(vmpn))
            # plt_Volt_Watt_pn = plot(vm_range, vw_curve.(vm_range,data_math["gen"]["1"]["pmax"][1]), label=false, guidefontsize=10)
            plt = plot(vm_range, vw_curve.(vm_range,data_math["gen"]["1"]["pmax"][1]), label=false, guidefontsize=10)
            plot!([vmpn_vals[1]], [(pg_values[1,gen_id])], seriestype=:scatter, label=L"(p_{an}, U_{an})", xticks=(0.85:0.05:1.15, 0.85:0.05:1.15), xrotation = 45)
            plot!([vmpn_vals[2]], [(pg_values[2,gen_id])], seriestype=:scatter, label=L"(p_{bn}, U_{bn})")
            plot!([vmpn_vals[3]], [(pg_values[3,gen_id])], seriestype=:scatter, label=L"(p_{cn}, U_{cn})")
            xlabel!("Voltage (V pu)")
            ylabel!("Active Power (kW)")

        elseif control == "Volt_Watt_pp"
            vmpp_vals = Array(value.(vmpp))
            # plt_Volt_Watt_pp = plot(vm_range, vw_curve.(vm_range,data_math["gen"]["1"]["pmax"][1]), label=false, guidefontsize=10)
            plt = plot(vm_range, vw_curve.(vm_range,data_math["gen"]["1"]["pmax"][1]), label=false, guidefontsize=10)
            plot!([vmpp_vals[1]/sqrt(3)], [pg_values[1,gen_id]-pg_values[2,gen_id]], seriestype=:scatter, label=L"(p_{ab}, U_{ab})", xticks=(0.85:0.05:1.15, 0.85:0.05:1.15), xrotation = 45)
            plot!([vmpp_vals[2]/sqrt(3)], [pg_values[2,gen_id]-pg_values[3,gen_id]], seriestype=:scatter, label=L"(p_{bc}, U_{bc})")
            plot!([vmpp_vals[3]/sqrt(3)], [pg_values[3,gen_id]-pg_values[1,gen_id]], seriestype=:scatter, label=L"(p_{ca}, U_{ca})")
            xlabel!("Voltage (V pu)")
            # ylabel!("Active Power (kW)")

        elseif control == "Volt_Watt_av"
            vm_vals = Array(value.(vm))
            # plt_Volt_Watt_av = plot(vm_range, vw_curve.(vm_range, data_math["gen"]["1"]["pmax"][1]), label=false, guidefontsize=10)
            plt = plot(vm_range, vw_curve.(vm_range, data_math["gen"]["1"]["pmax"][1]), label=false, guidefontsize=10)
            plot!([sum(vm_vals[phases])/3], [pg_values[1,gen_id]], seriestype=:scatter, label=L"(p_p, (U_a+U_b+U_c)/3)", legend=:bottomleft, xticks=(0.85:0.05:1.15, 0.85:0.05:1.15), xrotation = 45)
            xlabel!("Voltage (V pu)")
            # ylabel!("Active Power (kW)")
        end
    end

    vv_vw_plots[control] = plt

end

##

using Plots
mkpath("./Figures")

vv_vw_plots = Plots.plot(vv_vw_plots["Volt_var_pn"], vv_vw_plots["Volt_var_pp"], vv_vw_plots["Volt_var_av"], vv_vw_plots["Volt_Watt_pn"], vv_vw_plots["Volt_Watt_pp"], vv_vw_plots["Volt_Watt_av"], layout=(2,3))
Plots.savefig(vv_vw_plots, "./Figures/vv_vw_plots.pdf")