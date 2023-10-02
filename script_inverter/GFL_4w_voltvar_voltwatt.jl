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


include("./inverter_loss_branch.jl")
add_inverter_losses(data_math, gen_id)

data_math["gen"]["1"]["cost"] = [10 0]
data_math["gen"]["2"]["cost"] = [1000 0]

ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]

###
vv_vw_plots = Dict()

for control in ["Volt_Watt_pn", "Volt_var_pn", "Volt_Watt_pp", "Volt_var_pp", "Volt_Watt_av", "Volt_var_av"]

# control = "Volt_Watt_av"

    model = JuMP.Model(ipopt_solver)

    # function IVR_EN(model, ref)
    terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref[:bus])
    v_start = [exp.(im.*collect(0:-1:-2)*2/3*pi) ; 0]

    n_ph = 4
    vr = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vr", start = real(v_start)[t], lower_bound = -ref[:bus][i]["vmax"][t], upper_bound = ref[:bus][i]["vmax"][t] ) for i in keys(ref[:bus]))
    vr = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vr[i].axes[1] ? vr[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
    vi = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vi", start = imag(v_start)[t], lower_bound = -ref[:bus][i]["vmax"][t], upper_bound = ref[:bus][i]["vmax"][t] ) for i in keys(ref[:bus]))
    vi = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vi[i].axes[1] ? vi[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))

    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref[:branch])
    conds = Dict(l => branch["f_connections"] for (l,branch) in ref[:branch])
    cr = Dict((l,i,j) => JuMP.@variable(model, [c in conds[l]], base_name="cr_$((l,i,j))") for (l,i,j) in ref[:arcs_branch]) # , lower_bound = -ref[:branch][l]["c_rating_a"], upper_bound = ref[:branch][l]["c_rating_a"]
    cr = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in cr[(l,i,j)].axes[1] ? cr[(l,i,j)][c] : 0.0 for c in 1:n_ph, (l,i,j) in ref[:arcs_branch]]), 1:n_ph, ref[:arcs_branch])
    ci = Dict((l,i,j) => JuMP.@variable(model, [c in conds[l]], base_name="ci_$((l,i,j))") for (l,i,j) in ref[:arcs_branch]) # , lower_bound = -ref[:branch][l]["c_rating_a"], upper_bound = ref[:branch][l]["c_rating_a"]
    ci = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in ci[(l,i,j)].axes[1] ? ci[(l,i,j)][c] : 0.0 for c in 1:n_ph, (l,i,j) in ref[:arcs_branch]]), 1:n_ph, ref[:arcs_branch])

    csr = Dict(l => JuMP.@variable(model, [c in conds[l]], base_name="csr_$l") for (l,i,j) in ref[:arcs_branch])
    csi = Dict(l => JuMP.@variable(model, [c in conds[l]], base_name="csi_$l") for (l,i,j) in ref[:arcs_branch])
    csr = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in csr[l].axes[1] ? csr[l][c] : 0.0 for c in 1:n_ph, l in keys(ref[:branch])]), 1:n_ph, keys(ref[:branch]))
    csi = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in csi[l].axes[1] ? csi[l][c] : 0.0 for c in 1:n_ph, l in keys(ref[:branch])]), 1:n_ph, keys(ref[:branch]))
    cr_bus = Dict{Tuple{Int,Int,Int}, Any}()
    ci_bus = Dict{Tuple{Int,Int,Int}, Any}()

    int_dim = Dict(i => RPMD._infer_int_dim_unit(gen, false) for (i,gen) in ref[:gen])
    crg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="crg_$i") for i in keys(ref[:gen]))
    cig = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cig_$i") for i in keys(ref[:gen]))
    crg = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? crg[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
    cig = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? cig[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))

    pg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="pg_$i") for i in keys(ref[:gen]))
    qg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="qg_$i") for i in keys(ref[:gen]))
    pg = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? pg[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
    qg = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? qg[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
    crg_bus = Dict{Int, Any}()
    cig_bus = Dict{Int, Any}()

    int_dim = Dict(i => RPMD._infer_int_dim_unit(load, false) for (i,load) in ref[:load])
    connections = Dict(i => load["connections"][1:int_dim[i]] for (i, load) in ref[:load])
    crd = Dict(i => JuMP.@variable(model, [c in connections[i]], base_name="crd_$i") for i in keys(ref[:load]))
    cid = Dict(i => JuMP.@variable(model, [c in connections[i]], base_name="cid_$i") for i in keys(ref[:load]))
    crd = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in crd[i].axes[1] ? crd[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:load])]), 1:n_ph, keys(ref[:load]))
    cid = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in cid[i].axes[1] ? cid[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:load])]), 1:n_ph, keys(ref[:load]))
    crd_bus = Dict{Int, Any}()
    cid_bus = Dict{Int, Any}()


    ###
    for (i, bus) in ref[:bus]
        terminals = bus["terminals"]
        grounded = bus["grounded"]
        
        if i in keys(ref[:ref_buses])
            vref = bus["vm"][bus["terminals"]] .* exp.(im*bus["va"][bus["terminals"]])
            vrefre = real.(vref)
            vrefim = imag.(vref)
            JuMP.@constraint(model, vr[:,i] .== vrefre)
            JuMP.@constraint(model, vi[:,i] .== vrefim)
        end

        ungrounded_terminals = [t for (idx,t) in enumerate(terminals) if !grounded[idx]]
        nonzero_vmin_terminals = ungrounded_terminals[bus["vmin"][ungrounded_terminals] .> 0]
        JuMP.@constraint(model, Vector{JuMP.AffExpr}(vr[nonzero_vmin_terminals,i]).^2 .+ Vector{JuMP.AffExpr}(vi[nonzero_vmin_terminals,i]).^2 .>= bus["vmin"][nonzero_vmin_terminals].^2)

        nonInf_vmax_terminals = ungrounded_terminals[bus["vmax"][ungrounded_terminals] .< Inf]
        JuMP.@constraint(model, Vector{JuMP.AffExpr}(vr[nonInf_vmax_terminals,i]).^2 .+ Vector{JuMP.AffExpr}(vi[nonInf_vmax_terminals,i]).^2 .<= bus["vmax"][nonInf_vmax_terminals].^2)
    end


    for (id, generator) in ref[:gen]
        nphases = RPMD._infer_int_dim_unit(generator, false)
        bus_id = generator["gen_bus"]
        bus = ref[:bus][bus_id]
        configuration = generator["configuration"]
        connections = generator["connections"]

        N = length(connections)
        pmin = get(generator, "pmin", fill(-Inf, N))
        pmax = get(generator, "pmax", fill( Inf, N))
        qmin = get(generator, "qmin", fill(-Inf, N))
        qmax = get(generator, "qmax", fill( Inf, N))

        # constraint_mc_generator_current(pm, id)
        if configuration==PMD.WYE || length(pmin)==1 || nphases==1
            phases = connections[1:end-1]
            n = connections[end]

            crg_bus[id] = JuMP.Containers.DenseAxisArray([crg[phases,id]..., -sum(crg[phases,id])], connections)
            cig_bus[id] = JuMP.Containers.DenseAxisArray([cig[phases,id]..., -sum(cig[phases,id])], connections)

            JuMP.@constraint(model,  pg[phases,id] .==  (vr[phases,bus_id] .- vr[n,bus_id]) .* crg[phases,id] .+ (vi[phases,bus_id] .- vi[n,bus_id]) .* cig[phases,id])
            nonInf_pmin_pmax = [idx for (idx,t) in enumerate(pmin) if pmin[idx].>-Inf || pmax[idx].<Inf]
            JuMP.@constraint(model, pmin[nonInf_pmin_pmax] .<= pg[nonInf_pmin_pmax,id] .<= pmax[nonInf_pmin_pmax] )
            # JuMP.@constraint(model, pmin .<= pg[phases,id] .<= pmax )

            nonInf_qmin_qmax = [idx for (idx,t) in enumerate(qmin) if qmin[idx].>-Inf || qmax[idx].<Inf]
            JuMP.@constraint(model, qg[nonInf_qmin_qmax,id] .== -(vr[nonInf_qmin_qmax,bus_id] .- vr[n,bus_id]) .* cig[nonInf_qmin_qmax,id] .+ (vi[nonInf_qmin_qmax,bus_id] .- vi[n,bus_id]) .* crg[nonInf_qmin_qmax,id])
            JuMP.@constraint(model, qmin[nonInf_qmin_qmax] .<= qg[nonInf_qmin_qmax,id] .<= qmax[nonInf_qmin_qmax] )
            # JuMP.@constraint(model, qmin .<= qg[phases,id] .<= qmax )

            # for (idx, p) in enumerate(phases)
            #     if pmin[idx]>-Inf
            #         JuMP.@constraint(model, pmin[idx] .<= pg[p,id])
            #     end
            #     if pmax[idx]< Inf
            #         JuMP.@constraint(model, pmax[idx] .>= pg[p,id])
            #     end
            #     if qmin[idx]>-Inf || qmax[idx]< Inf
            #         JuMP.@constraint(model, qg[p,id] .== -(vr[p,bus_id] .- vr[n,bus_id]) .* cig[p,id] .+ (vi[p,bus_id] .- vi[n,bus_id]) .* crg[p,id])
            #         if qmin[idx]>-Inf
            #             JuMP.@constraint(model, qmin[idx] .<= qg[p,id])
            #         end
            #         if qmax[idx]< Inf
            #             JuMP.@constraint(model, qmax[idx] .>= qg[p,id])
            #         end
            #     end
            # end

        else ## configuration==PMD.DELTA

            Md = PMD._get_delta_transformation_matrix(length(connections))
            crg_bus[id] = JuMP.Containers.DenseAxisArray(Md'*Vector{JuMP.AffExpr}(crg[connections,id]), connections)
            cig_bus[id] = JuMP.Containers.DenseAxisArray(Md'*Vector{JuMP.AffExpr}(cig[connections,id]), connections)

            connections2 = [connections[2:end]..., connections[1]]
            vrg = Vector{JuMP.AffExpr}(vr[connections,bus_id]) .- Vector{JuMP.AffExpr}(vr[connections2,bus_id])
            vig = Vector{JuMP.AffExpr}(vi[connections,bus_id]) .- Vector{JuMP.AffExpr}(vi[connections2,bus_id])

            JuMP.@constraint(model, pg[connections,id] .==  vrg .* crg[connections,id] .+ vig .* cig[connections,id])
            JuMP.@constraint(model, qg[connections,id] .== -vrg .* cig[connections,id] .+ vig .* crg[connections,id])

            JuMP.@constraint(model, pmin[connections] .<= pg[connections,id] .<= pmax[connections])
            JuMP.@constraint(model, qmin[connections] .<= qg[connections,id] .<= qmax[connections])
        end
    end
    crg_bus = JuMP.Containers.DenseAxisArray([t in crg_bus[i].axes[1] ? crg_bus[i][t] : 0 for t in 1:n_ph, i in keys(ref[:gen])], 1:n_ph, keys(ref[:gen]))
    cig_bus = JuMP.Containers.DenseAxisArray([t in cig_bus[i].axes[1] ? cig_bus[i][t] : 0 for t in 1:n_ph, i in keys(ref[:gen])], 1:n_ph, keys(ref[:gen]))


    for (i, branch) in ref[:branch]
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_idx = (i, f_bus, t_bus)
        t_idx = (i, t_bus, f_bus)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        r = branch["br_r"]
        x = branch["br_x"]
        f_connections =  branch["f_connections"]
        t_connections =  branch["t_connections"]
        c_rating = branch["c_rating_a"]

        vr_fr = [vr[idx,f_bus] for (idx,v) in enumerate(vr[:,f_bus])]
        vi_fr = [vi[idx,f_bus] for (idx,v) in enumerate(vi[:,f_bus])]
        vr_to = [vr[idx,t_bus] for (idx,v) in enumerate(vr[:,t_bus])]
        vi_to = [vi[idx,t_bus] for (idx,v) in enumerate(vi[:,t_bus])]

        cr_fr = [cr[idx,f_idx] for (idx,c) in enumerate(cr[:,f_idx])]
        ci_fr = [ci[idx,f_idx] for (idx,c) in enumerate(ci[:,f_idx])]
        cr_to = [cr[idx,t_idx] for (idx,c) in enumerate(cr[:,t_idx])]
        ci_to = [ci[idx,t_idx] for (idx,c) in enumerate(ci[:,t_idx])]

        csr_fr = [csr[idx,f_idx[1]] for (idx,c) in enumerate(csr[:,f_idx[1]])]
        csi_fr = [csi[idx,f_idx[1]] for (idx,c) in enumerate(csi[:,f_idx[1]])]
        csr_to = -[csr[idx,t_idx[1]] for (idx,c) in enumerate(csr[:,t_idx[1]])]
        csi_to = -[csi[idx,t_idx[1]] for (idx,c) in enumerate(csi[:,t_idx[1]])]

        ### constraint_mc_current_from
        JuMP.@constraint(model, cr_fr .== csr_fr .+ g_fr * vr_fr .- b_fr * vi_fr)
        JuMP.@constraint(model, ci_fr .== csi_fr .+ g_fr * vi_fr .+ b_fr * vr_fr)
        cr_bus[f_idx] = JuMP.Containers.DenseAxisArray(cr_fr, f_connections)
        ci_bus[f_idx] = JuMP.Containers.DenseAxisArray(ci_fr, f_connections)

        ### constraint_mc_current_to
        JuMP.@constraint(model, cr_to .== csr_to .+ g_to * vr_to .- b_to * vi_to)
        JuMP.@constraint(model, ci_to .== csi_to .+ g_to * vi_to .+ b_to * vr_to)
        cr_bus[t_idx] = JuMP.Containers.DenseAxisArray(cr_to, t_connections)
        ci_bus[t_idx] = JuMP.Containers.DenseAxisArray(ci_to, t_connections)

        ### constraint_mc_bus_voltage_drop
        JuMP.@constraint(model, vr_to .== vr_fr .- r*csr_fr .+ x*csi_fr)
        JuMP.@constraint(model, vi_to .== vi_fr .- r*csi_fr .- x*csr_fr)

        ### constraint_mc_branch_current_limit
        cnds_finite_rating = [c for (c,r) in enumerate(c_rating) if r<Inf]
        JuMP.@constraint(model, cr_fr[cnds_finite_rating].^2 .+ ci_fr[cnds_finite_rating].^2 .<= c_rating[cnds_finite_rating].^2)
        JuMP.@constraint(model, cr_to[cnds_finite_rating].^2 .+ ci_to[cnds_finite_rating].^2 .<= c_rating[cnds_finite_rating].^2)

        ### constraint_mc_thermal_limit
        if haskey(branch, "rate_a") && any(branch["rate_a"] .< Inf)
            ### constraint_mc_thermal_limit_from
            pf_idx = JuMP.@expression(model,  vr_fr .* cr_fr .+ vi_fr .* ci_fr)
            qf_idx = JuMP.@expression(model, -vr_fr .* ci_fr .+ vi_fr .* cr_fr)
            JuMP.@constraint(model, pf_idx.^2 .+ qf_idx.^2 .<= branch["rate_a"].^2)

            ### constraint_mc_thermal_limit_to
            pt_idx = JuMP.@expression(model,  vr_to .* cr_to .+ vi_to .* ci_to)
            qt_idx = JuMP.@expression(model, -vr_to .* ci_to .+ vi_to .* cr_to)
            JuMP.@constraint(model, pt_idx.^2 .+ qt_idx.^2 .<= branch["rate_a"].^2)
        end

    end
    cr_bus = JuMP.Containers.DenseAxisArray([cr_bus[(l,i,j)][t] for t in 1:n_ph , (l,i,j) in ref[:arcs_branch]], 1:n_ph, ref[:arcs_branch])
    ci_bus = JuMP.Containers.DenseAxisArray([ci_bus[(l,i,j)][t] for t in 1:n_ph , (l,i,j) in ref[:arcs_branch]], 1:n_ph, ref[:arcs_branch])


    for (id, load) in ref[:load]
        bus_id = load["load_bus"]
        bus = ref[:bus][bus_id]
        configuration = load["configuration"]
        connections = load["connections"]
        load_model = load["model"]
        a, alpha, b, beta = PMD._load_expmodel_params(load, bus)

        int_dim = RPMD._infer_int_dim_unit(load, false)
        if configuration==PMD.WYE || int_dim==1
            phases = connections[1:end-1]
            n = connections[end]

            vr_pn = Vector{JuMP.AffExpr}(vr[phases,bus_id] .- vr[n,bus_id])
            vi_pn = Vector{JuMP.AffExpr}(vi[phases,bus_id] .- vi[n,bus_id])

            if load_model==PMD.POWER
                pd = a
                qd = b
            elseif load_model==PMD.IMPEDANCE
                pd = a .* (vr_pn.^2 .+ vi_pn.^2)
                qd = b .* (vr_pn.^2 .+ vi_pn.^2)
            elseif load_model==PMD.CURRENT
                pd = JuMP.@variable(model, [c in 1:int_dim])
                qd = JuMP.@variable(model, [c in 1:int_dim])
                JuMP.@constraint(model, pd.^2 .== a.^2 .* (vr_pn.^2 .+ vi_pn.^2))
                JuMP.@constraint(model, sign.(a).*pd .>= 0)
                JuMP.@constraint(model, qd.^2 .== b.^2 .* (vr_pn.^2 .+ vi_pn.^2))
                JuMP.@constraint(model, sign.(b).*qd .>= 0)
            else
                error("Load model $model for load $id is not supported by this formulation.")
            end

            JuMP.@constraint(model, Vector{JuMP.AffExpr}(crd[phases,id]) .== 
                    a .* vr_pn .* (vr_pn.^2 .+ vi_pn.^2).^(alpha/2 .-1)
                .+ b .* vi_pn .* (vr_pn.^2 .+ vi_pn.^2).^(beta/2 .-1))
            JuMP.@constraint(model, Vector{JuMP.AffExpr}(cid[phases,id]) .== 
                    a .* vi_pn .* (vr_pn.^2 .+ vi_pn.^2).^(alpha/2 .-1)
                .- b .* vr_pn .* (vr_pn.^2 .+ vi_pn.^2).^(beta/2 .-1))
            # JuMP.@constraint(model, pd .==  vr_pn .* Vector{JuMP.AffExpr}(crd[phases,id]) .+ vi_pn .* Vector{JuMP.AffExpr}(cid[p,id]))
            # JuMP.@constraint(model, qd .== -vr_pn .* Vector{JuMP.AffExpr}(cid[p,id]) .+ vi_pn .* Vector{JuMP.AffExpr}(crd[p,id]))
            
            crd_bus[id] = JuMP.Containers.DenseAxisArray([Vector{JuMP.AffExpr}(crd[phases,id])..., -sum(Vector{JuMP.AffExpr}(crd[phases,id]))], connections)
            cid_bus[id] = JuMP.Containers.DenseAxisArray([Vector{JuMP.AffExpr}(cid[phases,id])..., -sum(Vector{JuMP.AffExpr}(cid[phases,id]))], connections)
            
        else
            phases = connections
            phases_next = [connections[2:end]..., connections[1]]
            P = length(connections)
            idxs = 1:P
            idxs_prev = [idxs[end], idxs[1:end-1]...]
            
            vrd = Vector{JuMP.AffExpr}(vr[phases,bus_id]) .- Vector{JuMP.AffExpr}(vr[phases_next,bus_id])
            vid = Vector{JuMP.AffExpr}(vi[phases,bus_id]) .- Vector{JuMP.AffExpr}(vi[phases_next,bus_id])

            if load_model==PMD.POWER
                pd = a
                qd = b
            elseif load_model==PMD.IMPEDANCE
                pd = a .* (vrd.^2 .+ vid.^2)
                qd = b .* (vrd.^2 .+ vid.^2)
            elseif load_model==PMD.CURRENT
                pd = JuMP.@variable(model, [c in 1:int_dim])
                qd = JuMP.@variable(model, [c in 1:int_dim])
                JuMP.@constraint(model, pd[id].^2 .== a.^2 .* (vrd.^2 .+ vid.^2))
                JuMP.@constraint(model, sign.(a).*pd[id] .>= 0)
                JuMP.@constraint(model, qd[id].^2 .== b.^2 .* (vrd.^2 .+ vid.^2))
                JuMP.@constraint(model, sign.(b).*qd[id] .>= 0)
            else
                error("Load model $model for load $id is not supported by this formulation.")
            end

            JuMP.@constraint(model, Vector{JuMP.AffExpr}(crd[idxs,id]) .== 
                        a[idxs] .* vrd[idxs,bus_id] .* (vrd[idxs,bus_id].^2 .+ vid[idxs,bus_id].^2).^(alpha[idxs]/2 .-1) .+ 
                        b[idxs] .* vid[idxs,bus_id] .* (vrd[idxs,bus_id].^2 .+ vid[idxs,bus_id].^2).^(beta[idxs]/2 .-1))
            JuMP.@constraint(model, Vector{JuMP.AffExpr}(cid[idxs,id]) .== 
                        a[idxs] .* vid[idxs,bus_id] .* (vrd[idxs,bus_id].^2 .+ vid[idxs,bus_id].^2).^(alpha[idxs]/2 .-1) .- 
                        b[idxs] .* vrd[idxs,bus_id] .* (vrd[idxs,bus_id].^2 .+ vid[idxs,bus_id].^2).^(beta[idxs]/2 .-1))
            # JuMP.@constraint(model, pd[idxs] .==  vrd[idxs,bus_id] .* Vector{JuMP.AffExpr}(crd[idxs,id]) .+ vid[idxs,bus_id] .* Vector{JuMP.AffExpr}(cid[idxs,id]))
            # JuMP.@constraint(model, qd[idxs] .== -vrd[idxs,bus_id] .* Vector{JuMP.AffExpr}(cid[idxs,id]) .+ vid[idxs,bus_id] .* Vector{JuMP.AffExpr}(crd[idxs,id]))

            crd_bus[id] = JuMP.Containers.DenseAxisArray(Vector{JuMP.AffExpr}(crd[idxs,id]).-Vector{JuMP.AffExpr}(crd[idxs_prev,id]), connections)
            cid_bus[id] = JuMP.Containers.DenseAxisArray(Vector{JuMP.AffExpr}(cid[idxs,id]).-Vector{JuMP.AffExpr}(cid[idxs_prev,id]), connections)
        end
        
    end
    crd_bus = JuMP.Containers.DenseAxisArray([t in crd_bus[i].axes[1] ? crd_bus[i][t] : 0 for t in 1:n_ph, i in keys(ref[:load])], 1:n_ph, keys(ref[:load]))
    cid_bus = JuMP.Containers.DenseAxisArray([t in cid_bus[i].axes[1] ? cid_bus[i][t] : 0 for t in 1:n_ph, i in keys(ref[:load])], 1:n_ph, keys(ref[:load]))


    for (i, bus) in ref[:bus]
        ### constraint_mc_current_balance
        bus_arcs = ref[:bus_arcs_conns_branch][i]
        bus_arcs_sw = ref[:bus_arcs_conns_switch][i]
        bus_arcs_trans = ref[:bus_arcs_conns_transformer][i]
        bus_gens = ref[:bus_conns_gen][i]
        bus_storage = ref[:bus_conns_storage][i]
        bus_loads = ref[:bus_conns_load][i]
        bus_shunts = ref[:bus_conns_shunt][i]

        terminals = bus["terminals"]
        grounded = bus["grounded"]

        Gt, Bt = RPMD._build_bus_shunt_matrices(ref, terminals, bus_shunts)

        ungrounded_terminals = [t for (idx,t) in enumerate(terminals) if !grounded[idx]]

        JuMP.@constraint(model, sum(Vector{JuMP.AffExpr}(cr_bus[ungrounded_terminals,a]) for (a, conns) in bus_arcs) .==
                                sum(Vector{JuMP.AffExpr}(crg_bus[ungrounded_terminals,g]) for (g, conns) in bus_gens) .- 
                                sum(Vector{JuMP.AffExpr}(crd_bus[ungrounded_terminals,d]) for (d, conns) in bus_loads) .-
                                Gt[ungrounded_terminals,ungrounded_terminals] * Vector{JuMP.AffExpr}(vr[ungrounded_terminals,i]) .- 
                                    Bt[ungrounded_terminals,ungrounded_terminals] * Vector{JuMP.AffExpr}(vi[ungrounded_terminals,i])
                                )
        
        JuMP.@constraint(model, sum(Vector{JuMP.AffExpr}(ci_bus[ungrounded_terminals,a]) for (a, conns) in bus_arcs) .==
                                sum(Vector{JuMP.AffExpr}(cig_bus[ungrounded_terminals,g]) for (g, conns) in bus_gens) .- 
                                sum(Vector{JuMP.AffExpr}(cid_bus[ungrounded_terminals,d]) for (d, conns) in bus_loads) .-
                                Gt[ungrounded_terminals,ungrounded_terminals] * Vector{JuMP.AffExpr}(vi[ungrounded_terminals,i]) .+
                                    Bt[ungrounded_terminals,ungrounded_terminals] * Vector{JuMP.AffExpr}(vr[ungrounded_terminals,i])
                                )
        
    end


    # alpha = exp(im*2/3*pi)
    # T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
    # Tre = real.(T)
    # Tim = imag.(T)

    ### objectives: cost, loss, VUF, VUF2, PVUR, LVUR, IUF, IUF2, PIUR, PPUR (x), PQUR (x)
    objective = "cost"

    if objective == "cost"
        JuMP.@objective(model, Min, sum(gen["cost"][1]*sum(pg[:,i]) + gen["cost"][2] for (i,gen) in ref[:gen]))

    elseif objective == "loss"
        branch_id = [i for (i, branch) in ref[:branch] if occursin("internal", branch["name"])][1]
        br_r = ref[:branch][branch_id]["br_r"]
        br_x = ref[:branch][branch_id]["br_x"]
        JuMP.@objective(model, Min, sum(sqrt.(br_r.^2 .+ br_x.^2) * Array(csr[:,branch_id].^2 .+ csi[:,branch_id].^2) ))

    elseif objective in ["VUF" "VUF2" "PVUR" "LVUR"]
        # terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref[:bus])
        # terminals = Dict(i => collect(1:3) for (i, bus) in ref[:bus])
        # vm = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm_$i", lower_bound = 0 ) for i in keys(ref[:bus]))
        # vm = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm[i].axes[1] ? vm[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
        # vr_012 = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vr_012") for i in keys(ref[:bus]))
        # vr_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vr_012[i].axes[1] ? vr_012[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
        # vi_012 = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vi_012") for i in keys(ref[:bus]))
        # vi_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vi_012[i].axes[1] ? vi_012[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
        # vm_012 = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm_012") for i in keys(ref[:bus]))
        # vm_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm_012[i].axes[1] ? vm_012[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))

        vm = JuMP.@variable(model, [t in 1:3], base_name="vm", lower_bound = 0 )
        vr_012 = JuMP.@variable(model, [t in 1:3], base_name="vr_012")
        vi_012 = JuMP.@variable(model, [t in 1:3], base_name="vi_012")
        vm_012 = JuMP.@variable(model, [t in 1:3], base_name="vm_012", lower_bound = 0)

        # for (i, bus) in ref[:bus]
            i = 1
            bus = ref[:bus][i]
            # JuMP.@constraint(model, vm[:,i].^2 .== vr[:,i].^2 .+ vi[:,i].^2)
            terminals = ref[:bus][i]["terminals"][1:end-1]
            JuMP.@constraint(model, vm[terminals,i].^2 .== vr[terminals,i].^2 .+ vi[terminals,i].^2)
            JuMP.@constraint(model, vr_012[terminals,i] .== Tre * Array(vr[terminals,i]) .- Tim * Array(vi[terminals,i]))
            JuMP.@constraint(model, vi_012[terminals,i] .== Tre * Array(vi[terminals,i]) .+ Tim * Array(vr[terminals,i]))
            JuMP.@constraint(model, vm_012[terminals,i].^2 .== vr_012[terminals,i].^2 .+ vi_012[terminals,i].^2)
        # end

        if objective == "VUF"
            JuMP.@objective(model, Min, vm_012[3,1] / vm_012[2,1])

        elseif objective == "VUF2"
            JuMP.@objective(model, Min, vm_012[3,1])

        elseif objective == "PVUR"
            phase_voltage = JuMP.@variable(model, [i=1], base_name="phase_voltage_$i")
            # phase_voltage = JuMP.@variable(model, [i in keys(ref[:bus])], base_name="phase_voltage_$i")
            # for (i, bus) in ref[:bus]
                i = 1
                bus = ref[:bus][i]
                terminals = ref[:bus][i]["terminals"][1:end-1]
                JuMP.@constraint(model, [t in terminals], phase_voltage[i] >= vm[t,i] - sum(vm[terminals,i])/3)
            # end
            JuMP.@objective(model, Min, phase_voltage[1] / (sum(vm[terminals,1])/3) )

        elseif objective == "LVUR"
            # line_voltage = JuMP.@variable(model, [i in keys(ref[:bus])], base_name="line_voltage_$i")
            # terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref[:bus])
            # vm_ll = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm_ll_$i") for i in keys(ref[:bus]))
            # vm_ll = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm_ll[i].axes[1] ? vm_ll[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
            line_voltage = JuMP.@variable(model, base_name="line_voltage")
            vm_ll = JuMP.@variable(model, [t in 1:3], base_name="vm_ll", lower_bound=0)
            vr_ll = JuMP.@variable(model, [t in 1:3], base_name="vr_ll")
            vi_ll = JuMP.@variable(model, [t in 1:3], base_name="vi_ll")
            # for (i, bus) in ref[:bus]
                # terminals = ref[:bus][i]["terminals"][1:end-1]
                # JuMP.@constraint(model, vm_ll[terminals,i] .== Array(vm[terminals,i]) .- Array(vm[terminals2,i]))
                # JuMP.@constraint(model, [t in terminals], line_voltage[i] >= vm_ll[t,i] - sum(vm_ll[terminals,i])/3)
                i = 1
                bus = ref[:bus][i]
                terminals = collect(1:3)
                terminals2 = [terminals[2:end]..., terminals[1]]
                # JuMP.@constraint(model, vm_ll[terminals] .== vm[terminals] .- vm[terminals2])
                JuMP.@constraint(model, vr_ll[terminals] .== Array(vr[terminals,i]) .- Array(vr[terminals2,i]))
                JuMP.@constraint(model, vi_ll[terminals] .== Array(vi[terminals,i]) .- Array(vi[terminals2,i]))
                JuMP.@constraint(model, vm_ll[terminals].^2 .== vr_ll[terminals].^2 .+ vi_ll[terminals].^2)

                JuMP.@constraint(model, line_voltage .>= vm_ll[terminals] .- sum(vm_ll[terminals])/3)
            # end
            JuMP.@objective(model, Min, line_voltage / (sum(vm_ll[terminals,1])/3) )
        end


    elseif objective in ["IUF" "IUF2" "PIUR"]
        int_dim = Dict(i => RPMD._infer_int_dim_unit(gen, !(4 in gen["connections"])) for (i,gen) in ref[:gen])
        cmg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cmg_$i", lower_bound=0) for i in keys(ref[:gen]))
        cmg = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? cmg[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
        crg_012 = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="crg_012_$i") for i in keys(ref[:gen]))
        crg_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? crg_012[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
        cig_012 = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cig_012_$i") for i in keys(ref[:gen]))
        cig_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? cig_012[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
        cmg_012 = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cmg_012_$i", lower_bound=0) for i in keys(ref[:gen]))
        cmg_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? cmg_012[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))

        for (i, gen) in ref[:gen]
            if 4 in gen["connections"]
                phases = ref[:gen][i]["connections"][1:end-1]
            else
                phases = ref[:gen][i]["connections"]
            end
            JuMP.@constraint(model, cmg[phases,i].^2 .== crg_bus[phases,i].^2 .+ crg_bus[phases,i].^2)
            JuMP.@constraint(model, crg_012[phases,i] .== Tre * Array(crg_bus[phases,i]) .- Tim * Array(cig_bus[phases,i]))
            JuMP.@constraint(model, cig_012[phases,i] .== Tre * Array(cig_bus[phases,i]) .+ Tim * Array(crg_bus[phases,i]))
            JuMP.@constraint(model, cmg_012[phases,i].^2 .== crg_012[phases,i].^2 .+ cig_012[phases,i].^2)
        end

        if objective == "IUF"
            JuMP.@objective(model, Min, cmg_012[3,1] / cmg_012[2,1])

        elseif objective == "IUF2"
            JuMP.@objective(model, Min, cmg_012[3,1])
            
        elseif objective == "PIUR"
            phase_current = JuMP.@variable(model, [i in keys(ref[:gen])], base_name="phase_current_$i")
            # for (i, gen) in ref[:gen]
                i = 1
                gen = ref[:gen][i]
                connections = ref[:gen][i]["connections"][1:end-1]
                JuMP.@constraint(model, [t in connections], phase_current[i] >= cmg[t,i] - sum(cmg[connections,i])/3)
            # end
            JuMP.@objective(model, Min, phase_current[1] / (sum(cmg[connections,1])/3) )
        end


    elseif objective == "PPUR"
        phase_p = JuMP.@variable(model, base_name="phase_p", lower_bound=0)
        # phase_p = JuMP.@variable(model, [i in keys(ref[:gen])], base_name="phase_p_$i")
        # for (i, gen) in ref[:gen]
            i = 1
            gen = ref[:gen][i]
            connections = gen["connections"][1:end-1]
            JuMP.@constraint(model, phase_p .>= pg[connections,i] .- sum(pg[connections,i])/3)
            # JuMP.@constraint(model, [t in connections], phase_p[i] >= pg[t,i] - sum(pg[connections,i])/3)
        # end
        # JuMP.@objective(model, Min, phase_p[1] / (sum(pg[connections,1])/3) )
        JuMP.@objective(model, Min, phase_p / (sum(pg[connections,1])/3) )


    elseif objective == "PQUR"
        phase_q = JuMP.@variable(model, [i=1], base_name="phase_q_$i")
        # phase_q = JuMP.@variable(model, [i in keys(ref[:gen])], base_name="phase_q_$i")
        # for (i, gen) in ref[:gen]
            i = 1
            gen = ref[:gen][i]
            connections = ref[:gen][i]["connections"][1:end-1]
            JuMP.@constraint(model, [c in connections], phase_q[i] >= qg[connections,i] - sum(qg[connections,i])/3)
        # end
        JuMP.@objective(model, Min, phase_q[1] / (sum(qg[connections,1])/3) )

    end
    
    ### control = Volt_Watt_pn, Volt_var_pn, Volt_Watt_pp, Volt_var_pp, Volt_Watt_av, Volt_var_av
    #### "Volt_Watt_pn" does not work with vm[phases]-vm[n], it is perhaps due to only vm[phases] in [0.9,1.1], and no bound on vm_pn
    # include("./VV_VW_controls.jl")

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
    # JuMP.@constraint(model, vmpn[2,bus_id]^2 == (vr[2,bus_id]-vr[4,bus_id])^2 + (vi[2,bus_id]-vi[4,bus_id])^2)
    # JuMP.@constraint(model, vmpn[3,bus_id]^2 == (vr[3,bus_id]-vr[4,bus_id])^2 + (vi[3,bus_id]-vi[4,bus_id])^2)

    # vv_curve(v, qmax) = v <= 0.95 ? qmax : (v >= 1.05 ? -qmax : -2*qmax/0.1*(v-1))
    # vw_curve(v, pmax) = v <= 0.95 ? pmax : (v >= 1.05 ? 0.0 : -pmax/0.1*(v-1.05))
    vv_curve(v, qmax) = v <= 0.9 ? qmax : (v >= 1.1 ? -qmax : -2*qmax/0.2*(v-1))
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
    vm_range = 0.8:0.01:1.2
    if control in ["Volt_var_pn", "Volt_var_pp", "Volt_var_av"]
        if control == "Volt_var_pn"
            vmpn_vals = Array(value.(vmpn))
            # plt_Volt_var_pn = plot(vm_range, vv_curve.(vm_range,data_math["gen"]["1"]["qmax"][1]), label=false, guidefontsize=10)
            plt = plot(vm_range, vv_curve.(vm_range,data_math["gen"]["1"]["qmax"][1]), label=false, guidefontsize=10)
            plot!([vmpn_vals[1]], [(qg_values[1,gen_id])], seriestype=:scatter, label=L"(q_{an}, U_{an})")
            plot!([vmpn_vals[2]], [(qg_values[2,gen_id])], seriestype=:scatter, label=L"(q_{bn}, U_{bn})")
            plot!([vmpn_vals[3]], [(qg_values[3,gen_id])], seriestype=:scatter, label=L"(q_{cn}, U_{cn})")
            # xlabel!("Voltage Magnitudes (V pu)")
            ylabel!("Reactive Power (kvar)")
            
        elseif control == "Volt_var_pp"
            vmpp_vals = Array(value.(vmpp))
            # plt_Volt_var_pp = plot(vm_range, vv_curve.(vm_range,data_math["gen"]["1"]["qmax"][1]), label=false, guidefontsize=10)
            plt = plot(vm_range, vv_curve.(vm_range,data_math["gen"]["1"]["qmax"][1]), label=false, guidefontsize=10)
            plot!([vmpp_vals[1]/sqrt(3)], [qg_values[1,gen_id]-qg_values[2,gen_id]], seriestype=:scatter, label=L"(q_{ab}, U_{ab})")
            plot!([vmpp_vals[2]/sqrt(3)], [qg_values[2,gen_id]-qg_values[3,gen_id]], seriestype=:scatter, label=L"(q_{bc}, U_{bc})")
            plot!([vmpp_vals[3]/sqrt(3)], [qg_values[3,gen_id]-qg_values[1,gen_id]], seriestype=:scatter, label=L"(q_{ca}, U_{ca})")
            # xlabel!("Voltage Magnitudes (V pu)")
            # ylabel!("Reactive Power (kvar)")

        elseif control == "Volt_var_av"
            vm_vals = Array(value.(vm))
            # plt_Volt_var_av = plot(vm_range, vv_curve.(vm_range,data_math["gen"]["1"]["qmax"][1]), label=false, guidefontsize=10)
            plt = plot(vm_range, vv_curve.(vm_range,data_math["gen"]["1"]["qmax"][1]), label=false, guidefontsize=10)
            plot!([sum(vm_vals[phases])/3], [qg_values[1,gen_id]], seriestype=:scatter, label=L"(q_p, (U_a+U_b+U_c)/3)", legend=:bottomleft)
            # xlabel!("Voltage Magnitudes (V pu)")
            # ylabel!("Reactive Power (kvar)")
            
        end

    elseif control in ["Volt_Watt_pn", "Volt_Watt_pp", "Volt_Watt_av"]
        if control == "Volt_Watt_pn"
            vmpn_vals = Array(value.(vmpn))
            # plt_Volt_Watt_pn = plot(vm_range, vw_curve.(vm_range,data_math["gen"]["1"]["pmax"][1]), label=false, guidefontsize=10)
            plt = plot(vm_range, vw_curve.(vm_range,data_math["gen"]["1"]["pmax"][1]), label=false, guidefontsize=10)
            plot!([vmpn_vals[1]], [(pg_values[1,gen_id])], seriestype=:scatter, label=L"(p_{an}, U_{an})")
            plot!([vmpn_vals[2]], [(pg_values[2,gen_id])], seriestype=:scatter, label=L"(p_{bn}, U_{bn})")
            plot!([vmpn_vals[3]], [(pg_values[3,gen_id])], seriestype=:scatter, label=L"(p_{cn}, U_{cn})")
            xlabel!("Voltage (V pu)")
            ylabel!("Active Power (kW)")

        elseif control == "Volt_Watt_pp"
            vmpp_vals = Array(value.(vmpp))
            # plt_Volt_Watt_pp = plot(vm_range, vw_curve.(vm_range,data_math["gen"]["1"]["pmax"][1]), label=false, guidefontsize=10)
            plt = plot(vm_range, vw_curve.(vm_range,data_math["gen"]["1"]["pmax"][1]), label=false, guidefontsize=10)
            plot!([vmpp_vals[1]/sqrt(3)], [pg_values[1,gen_id]-pg_values[2,gen_id]], seriestype=:scatter, label=L"(p_{ab}, U_{ab})")
            plot!([vmpp_vals[2]/sqrt(3)], [pg_values[2,gen_id]-pg_values[3,gen_id]], seriestype=:scatter, label=L"(p_{bc}, U_{bc})")
            plot!([vmpp_vals[3]/sqrt(3)], [pg_values[3,gen_id]-pg_values[1,gen_id]], seriestype=:scatter, label=L"(p_{ca}, U_{ca})")
            xlabel!("Voltage (V pu)")
            # ylabel!("Active Power (kW)")

        elseif control == "Volt_Watt_av"
            vm_vals = Array(value.(vm))
            # plt_Volt_Watt_av = plot(vm_range, vw_curve.(vm_range, data_math["gen"]["1"]["pmax"][1]), label=false, guidefontsize=10)
            plt = plot(vm_range, vw_curve.(vm_range, data_math["gen"]["1"]["pmax"][1]), label=false, guidefontsize=10)
            plot!([sum(vm_vals[phases])/3], [pg_values[1,gen_id]], seriestype=:scatter, label=L"(p_p, (U_a+U_b+U_c)/3)", legend=:bottomleft)
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

##

# vv_vw_plots = Plots.plot(plt_Volt_var_pn, plt_Volt_Watt_pn, plt_Volt_var_pp, plt_Volt_Watt_pp, plt_Volt_var_av, plt_Volt_Watt_av, layout=(3,2))
vv_vw_plots = Plots.plot(plt_Volt_var_pn, plt_Volt_var_av, plt_Volt_Watt_pn, plt_Volt_Watt_av, layout=(2,2))
Plots.savefig(vv_vw_plots, "./Figures/vv_vw_plots.pdf")


# vv_vw_plots["Volt_var_pp"]
##
alpha = exp(im*2/3*pi)
T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
Tre = real.(T)
Tim = imag.(T)

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


results_GFL_4w[objective] = Dict()
results_GFL_4w[objective]["vm_gen_012"] = vm1_012
results_GFL_4w[objective]["cm_gen_012"] = cgm1_012
results_GFL_4w[objective]["cm_load_012"] = cdm1_012
results_GFL_4w[objective]["cm_branch_012"] = cm1_012

# end