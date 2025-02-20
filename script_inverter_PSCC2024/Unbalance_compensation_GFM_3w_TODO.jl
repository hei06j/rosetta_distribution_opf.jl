using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP  # bl/array_nl
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels
import LinearAlgebra: diag, diagm


###
data_path = "./data/inverter_GFM_4w_wye_unbalanced_loads.dss"

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!])
# RPMD.pv1_correction!(data_eng)
data_eng["settings"]["sbase_default"] = 1
vbase = data_eng["settings"]["vbases_default"]["b1"]
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)
data_math_original = deepcopy(data_math)
data_math["bus"]["1"]["vmin"] = [0.9 * ones(3) ; 0 ]
data_math["bus"]["1"]["vmax"] = [1.1 * ones(3) ; Inf]

### - reomve the voltage source gen
### - remove the voltage source bus 
### - remove the voltage source branch that connects to the load bus

### - add the inverter internal bus: connect the inverter to this bus
###   - terminals = [1,2,3,4]
###   - grounded = [0,0,0,1]
### - add the inverter internal branch: f_bus and t_bus are new and old inverter buses

voltage_source_id = "2"
voltage_source_bus = data_math["gen"][voltage_source_id]["gen_bus"]
voltage_source_branch = [i for (i,branch) in data_math["branch"] if branch["f_bus"]==voltage_source_bus || branch["t_bus"]==voltage_source_bus][1]
delete!(data_math["gen"], "$voltage_source_id")
delete!(data_math["bus"], "$voltage_source_bus")
delete!(data_math["branch"], voltage_source_branch)
# data_math["bus"]["$voltage_source_bus"]["bus_type"] = 1
# data_math["bus"]["$voltage_source_bus"]["grounded"] = Bool[0, 0, 0, 0]


gen_id = 1
gen = data_math["gen"]["$gen_id"]
# gen["pmax"] = 23/3 * ones(3)
# gen["pmin"] = 23/3 * ones(3)
# gen["pmin"] = zeros(3)
# gen["qmax"] = sqrt.(25^2 - 23^2)/3 * ones(3)
# gen["qmin"] = -gen["qmax"]
# data_math["gen"]["1"]["cost"] = [0 0]
# data_math["gen"]["2"]["cost"] = [1000 0]

Rf = 0.015
Lf = 0.42E-3
Cf = 0.33E-9
old_gen_bus = gen["gen_bus"]
new_gen_bus = length(data_math["bus"]) + 1
# data_math["bus"]["$gen_bus"]["bus_type"] = 1
data_math["bus"]["$new_gen_bus"] = deepcopy(data_math["bus"]["$old_gen_bus"])
data_math["bus"]["$new_gen_bus"]["bus_i"] = new_gen_bus
data_math["bus"]["$new_gen_bus"]["index"] = new_gen_bus
data_math["bus"]["$new_gen_bus"]["name"] = "GFM_bus"
data_math["bus"]["$new_gen_bus"]["type"] = "GFM"
data_math["bus"]["$new_gen_bus"]["bus_type"] = 2
data_math["bus"]["$new_gen_bus"]["vm"] = [1. ; 1 ; 1 ; 0]
# data_math["bus"]["$new_gen_bus"]["va] = [0 ; -2. ; 1 ; 0]
data_math["bus"]["$new_gen_bus"]["grounded"] = Bool[0, 0, 0, 1]
data_math["bus"]["$old_gen_bus"]["bus_type"] = 1
gen["gen_bus"] = new_gen_bus

new_branch_id = length(data_math["branch"]) + 1
data_math["branch"]["$new_branch_id"] = deepcopy(data_math["branch"]["$(new_branch_id-1)"])
data_math["branch"]["$new_branch_id"]["index"] = new_branch_id
data_math["branch"]["$new_branch_id"]["name"] = "GFM_internal_z"
zbase = 230.94^2 / 1000
data_math["branch"]["$new_branch_id"]["br_r"] = diagm(Rf/zbase * ones(4))
data_math["branch"]["$new_branch_id"]["br_x"] = diagm(Lf*2*pi*50/zbase * ones(4))
data_math["branch"]["$new_branch_id"]["b_to"] = diagm(Cf*2*pi*50*zbase * ones(4))
data_math["branch"]["$new_branch_id"]["f_bus"] = new_gen_bus
data_math["branch"]["$new_branch_id"]["t_bus"] = old_gen_bus

ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]

##

results_GFM_4w = Dict()

# objective = "VUF2"
for objective in ["cost", "VUF", "VUF2", "PVUR", "LVUR", "IUF", "IUF2", "PIUR"]

    model = JuMP.Model(Ipopt.Optimizer)

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
    alpha = exp(im*2/3*pi)
    T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
    Tre = real.(T)
    Tim = imag.(T)

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

        # v_neg_seq_real = JuMP.@expression(model, Tre[3,:]' * Array(vr[1:3,i]) - Tim[3,:]' * Array(vi[1:3,i]))
        # v_neg_seq_imag = JuMP.@expression(model, Tre[3,:]' * Array(vi[1:3,i]) + Tim[3,:]' * Array(vr[1:3,i]))
        # JuMP.@constraint(model, v_neg_seq_real^2 + v_neg_seq_imag^2 <= 0.02^2)

        if haskey(bus, "type")
            if bus["type"] == "GFM"
                phases = terminals[1:end-1]
                n = terminals[end]
                # vm = JuMP.@variable(model, [t in phases], base_name="vm_$i", lower_bound=0)
                # JuMP.@constraint(model, [t in phases], vm[t]^2 == vr[t,i]^2 + vi[t,i]^2)
                JuMP.@constraint(model, [t in phases], vr[t,i]^2 + vi[t,i]^2 == bus["vm"][t]^2)
                # va_pn = JuMP.@variable(model, [t in terminals], base_name="va_pn_$i")
                # JuMP.@constraint(model, vm[phases] .- vm[n] .== bus["vm"][phases])
                JuMP.@constraint(model, (vr[1,i] * vi[2,i] - vi[1,i] * vr[2,i]) / (vr[1,i] * vr[2,i] + vi[1,i] * vi[2,i]) == tan(2*pi/3))
                JuMP.@constraint(model, (vr[2,i] * vi[3,i] - vi[2,i] * vr[3,i]) / (vr[2,i] * vr[3,i] + vi[2,i] * vi[3,i]) == tan(2*pi/3))
                # JuMP.@constraint(model, (vr[3,i] * vi[1,i] - vi[3,i] * vr[1,i]) / (vr[3,i] * vr[1,i] + vi[3,i] * vi[1,i]) == tan(2*pi/3))
                # JuMP.@constraint(model, vr[n,i] == 0.0)
                # JuMP.@constraint(model, vi[1,i] == 0.0)
            end
        end
        
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


    ### objectives: cost, VUF, VUF2, PVUR, LVUR, IUF, IUF2, PIUR, PPUR (x), PQUR (x)
    if objective == "cost"
        JuMP.@objective(model, Min, sum(gen["cost"][1]*sum(pg[:,i]) + gen["cost"][2] for (i,gen) in ref[:gen]))

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

    JuMP.optimize!(model)
    @assert(JuMP.termination_status(model) == LOCALLY_SOLVED)
    cost = JuMP.objective_value(model)

    ##
    gen_id = [parse(Int,i) for (i,gen) in data_math["gen"] if occursin("pv", gen["name"])][1]
    gen_bus_id = data_math["gen"]["$gen_id"]["gen_bus"]

    vr_vals = value.(vr)
    vi_vals = value.(vi)
    v1 = vr_vals[:,gen_bus_id] + im * vi_vals[:,gen_bus_id]
    vm1 = abs.(v1)
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

    branch_lij = (1, 1, 2)
    cr1 = value.(cr_bus)[:,branch_lij]
    ci1 = value.(ci_bus)[:,branch_lij]
    c1 = cr1 .+ im * ci1
    # cm1 = abs.(c1)
    c1_012 = T * c1[1:3]
    cm1_012 = abs.(c1_012)

    @show Array(cd1) .- Array(cg1)
    # cd1_012 .- cg1_012
    # [cd1_012  cg1_012]
    # cdm1_012 .- cgm1_012

    results_GFM_4w[objective] = Dict()
    results_GFM_4w[objective]["vm_gen_012"] = v1_012
    results_GFM_4w[objective]["cm_gen_012"] = cg1_012
    results_GFM_4w[objective]["cm_load_012"] = cd1_012
    results_GFM_4w[objective]["cm_branch_012"] = c1_012

end

##
pg_vals = JuMP.value.(pg)
qg_vals = JuMP.value.(qg)
vr_vals = JuMP.value.(vr)
vi_vals = JuMP.value.(vi)

vm = abs.(vr_vals .+ im*vi_vals)
va = angle.(vr_vals .+ im*vi_vals) .* 180/pi 


##
crg_values = JuMP.value.(crg_bus)
cig_values = JuMP.value.(cig_bus)
cmg_012_values = JuMP.value.(cmg_012)

vm_012_values = JuMP.value.(vm_012)

alpha = exp(im*2/3*pi)
T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
Tre = real.(T)
Tim = imag.(T)
val_cg = crg_values[1:3,1] .+ im .* cig_values[1:3,1]
val_cg_012 = T * Array(val_cg)
val_cmg_012 = abs.(val_cg_012)

# cmg_vv_mincost_balLoad = [0.0819961164734981
#  7.699458638290301
#  0.12776140632900412]

# cmg_vv_minI2_balLoad = [0.08199611647350194
# 7.699458638290301
# 0.12776140632900793]

# cmg_vv_mincost_unbalLoad = 0.8775286636572217
#  7.6564277512939665
#  1.1362190436422837]

# cmg_vv_minI2_unbalLoad = [0.8775286636572188
#  7.6564277512939665
#  1.1362190436422803

# cmg_NOvv_mincost_unbalLoad = 1.5255515209463981
#  7.726651234561468
#  1.3611366662907889

# cmg_NOvv_minI2_unbalLoad = [0.7097603618426963
#  7.8068514320686875
#  2.168646298832432e-6