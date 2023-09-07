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

# data_path = "./data/test_gen_3ph_wye.dss"         # works
# data_path = "./data/test_gen_3ph_delta.dss"
# data_path = "./data/test_gen_1ph_wye.dss"         # works
# data_path = "./data/test_gen_1ph_delta.dss"
# data_path = "./data/test_load_1ph_delta_cp.dss" 
# data_path = "./data/test_load_1ph_wye_cp.dss"     # works
# data_path = "./data/test_load_3ph_delta_cp.dss"
data_path = "./data/test_load_3ph_wye_cp.dss"      # works
# data_path = "./data/test_load_3ph_delta_ci.dss"
# data_path = "./data/test_load_3ph_wye_ci.dss"     # works
# data_path = "./data/test_load_3ph_delta_cz.dss"
# data_path = "./data/test_load_3ph_wye_cz.dss"     # works


ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!])
RPMD.pv1_correction!(data_eng)
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)

ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]

##
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
crd = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="crd_$i") for i in keys(ref[:load]))
cid = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cid_$i") for i in keys(ref[:load]))
crd = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? crd[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:load])]), 1:n_ph, keys(ref[:load]))
cid = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? cid[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:load])]), 1:n_ph, keys(ref[:load]))
crd_bus = Dict{Int, Any}()
cid_bus = Dict{Int, Any}()


##
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
    
    # ungrounded_terminals = terminals[(!).(grounded)]
    # for (idx,t) in enumerate(terminals)
    #     if !grounded[idx]
    #         if bus["vmax"][idx] < Inf
    #             JuMP.@constraint(model, vr[t,i]^2+vi[t,i]^2 <= bus["vmax"][idx]^2)
    #         end
    #         if bus["vmin"][idx] > 0.0
    #             JuMP.@constraint(model, vr[t,i]^2+vi[t,i]^2 >= bus["vmin"][idx]^2)
    #         end
    #     end
    # end

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
        # JuMP.@constraint(model, pmin .<= pg[phases,id] .<= pmax )
        nonInf_pmin_pmax = [idx for (idx,t) in enumerate(pmin) if pmin[idx].>-Inf || pmax[idx].<Inf]
        JuMP.@constraint(model, pmin[nonInf_pmin_pmax] .<= pg[nonInf_pmin_pmax,id] .<= pmax[nonInf_pmin_pmax] )


        nonInf_qmin_qmax = [idx for (idx,t) in enumerate(qmin) if qmin[idx].>-Inf || qmax[idx].<Inf]
        JuMP.@constraint(model,  qg[nonInf_qmin_qmax,id] .== -(vr[nonInf_qmin_qmax,bus_id] .- vr[n,bus_id]) .* cig[nonInf_qmin_qmax,id] .+ (vi[nonInf_qmin_qmax,bus_id] .- vi[n,bus_id]) .* crg[nonInf_qmin_qmax,id])
        JuMP.@constraint(model, qmin[nonInf_qmin_qmax] .<= qg[nonInf_qmin_qmax,id] .<= qmax[nonInf_qmin_qmax] )
        JuMP.@constraint(model, qmin .<= qg[phases,id] .<= qmax )

        for (idx, p) in enumerate(phases)
            if pmin[idx]>-Inf
                JuMP.@constraint(model, pmin[idx] .<= pg[p,id])
            end
            if pmax[idx]< Inf
                JuMP.@constraint(model, pmax[idx] .>= pg[p,id])
            end
            if qmin[idx]>-Inf || qmax[idx]< Inf
                JuMP.@constraint(model, qg[p,id] .== -(vr[p,bus_id] .- vr[n,bus_id]) .* cig[p,id] .+ (vi[p,bus_id] .- vi[n,bus_id]) .* crg[p,id])
                if qmin[idx]>-Inf
                    JuMP.@constraint(model, qmin[idx] .<= qg[p,id])
                end
                if qmax[idx]< Inf
                    JuMP.@constraint(model, qmax[idx] .>= qg[p,id])
                end
            end
        end

    else ## configuration==PMD.DELTA

        Md = PMD._get_delta_transformation_matrix(length(connections))
        # crg_bus[id] = _merge_bus_flows(model, Md'*crg[id], connections)
        # cig_bus[id] = _merge_bus_flows(model, Md'*cig[id], connections)

        crg_bus[id] = JuMP.Containers.DenseAxisArray(Md'*crg[id], connections)
        cig_bus[id] = JuMP.Containers.DenseAxisArray(Md'*cig[id], connections)

        nph = length(pmin)
        vrg = Dict()
        vig = Dict()
        for (idx,c,d) in zip(1:nph, connections, [connections[2:end]..., connections[1]])
            vrg[idx] = JuMP.@NLexpression(model, vr[bus_id][c]-vr[bus_id][d])
            vig[idx] = JuMP.@NLexpression(model, vi[bus_id][c]-vi[bus_id][d])
        end
        pg_expr = Vector{JuMP.NonlinearExpression}([])
        qg_expr = Vector{JuMP.NonlinearExpression}([])
        for idx in 1:nph
            JuMP.@NLconstraint(model, pg[id][idx] ==  vrg[idx]*crg[id][idx]+vig[idx]*cig[id][idx])
            JuMP.@NLconstraint(model, qg[id][idx] == -vrg[idx]*cig[id][idx]+vig[idx]*crg[id][idx])
        end
        JuMP.@NLconstraint(model, [i in 1:nph], pmin[i] <= pg[id][i])
        JuMP.@NLconstraint(model, [i in 1:nph], pmax[i] >= pg[id][i])
        JuMP.@NLconstraint(model, [i in 1:nph], qmin[i] <= qg[id][i])
        JuMP.@NLconstraint(model, [i in 1:nph], qmax[i] >= qg[id][i])
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
    # cr_bus[f_idx] = _merge_bus_flows(model, cr_fr, f_connections)
    # ci_bus[f_idx] = _merge_bus_flows(model, ci_fr, f_connections)
    cr_bus[f_idx] = JuMP.Containers.DenseAxisArray(cr_fr, f_connections)
    ci_bus[f_idx] = JuMP.Containers.DenseAxisArray(ci_fr, f_connections)

    ### constraint_mc_current_to
    JuMP.@constraint(model, cr_to .== csr_to .+ g_to * vr_to .- b_to * vi_to)
    JuMP.@constraint(model, ci_to .== csi_to .+ g_to * vi_to .+ b_to * vr_to)
    # cr_bus[t_idx] = _merge_bus_flows(model, cr_to, t_connections)
    # ci_bus[t_idx] = _merge_bus_flows(model, ci_to, t_connections)
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
        p = connections[1:end-1]
        n = connections[end]

        vr_pn = Vector{JuMP.AffExpr}(vr[p,bus_id] .- vr[n,bus_id])
        vi_pn = Vector{JuMP.AffExpr}(vi[p,bus_id] .- vi[n,bus_id])

        # crd_vec = [crd[idx,id] for (idx,c) in enumerate(crd[:,id])]
        # cid_vec = [cid[idx,id] for (idx,c) in enumerate(cid[:,id])]

        # vr_pn = [vr_pn[idx] for (idx,v) in enumerate(vr_pn)]
        # vi_pn = [vi_pn[idx] for (idx,v) in enumerate(vi_pn)]

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

        JuMP.@constraint(model, Vector{JuMP.AffExpr}(crd[p,id]) .== 
                a .* vr_pn .* (vr_pn.^2 .+ vi_pn.^2).^(alpha/2 .-1)
             .+ b .* vi_pn .* (vr_pn.^2 .+ vi_pn.^2).^(beta/2 .-1))
        JuMP.@constraint(model, Vector{JuMP.AffExpr}(cid[p,id]) .== 
                a .* vi_pn .* (vr_pn.^2 .+ vi_pn.^2).^(alpha/2 .-1)
             .- b .* vr_pn .* (vr_pn.^2 .+ vi_pn.^2).^(beta/2 .-1))
        # JuMP.@constraint(model, pd .==  vr_pn .* Vector{JuMP.AffExpr}(crd[p,id]) .+ vi_pn .* Vector{JuMP.AffExpr}(cid[p,id]))
        # JuMP.@constraint(model, qd .== -vr_pn .* Vector{JuMP.AffExpr}(cid[p,id]) .+ vi_pn .* Vector{JuMP.AffExpr}(crd[p,id]))
        crd_bus[id] = JuMP.Containers.DenseAxisArray([Vector{JuMP.AffExpr}(crd[p,id])..., -sum(Vector{JuMP.AffExpr}(crd[p,id]))], connections)
        cid_bus[id] = JuMP.Containers.DenseAxisArray([Vector{JuMP.AffExpr}(cid[p,id])..., -sum(Vector{JuMP.AffExpr}(cid[p,id]))], connections)

        # JuMP.@constraint(model, crd_vec[p] .== 
        #         a .* vr_pn .* (vr_pn.^2 .+ vi_pn.^2).^(alpha/2 .-1)
        #      .+ b .* vi_pn .* (vr_pn.^2 .+ vi_pn.^2).^(beta/2 .-1))
        # JuMP.@constraint(model, cid_vec[p] .== 
        #         a .* vi_pn .* (vr_pn.^2 .+ vi_pn.^2).^(alpha/2 .-1)
        #      .- b .* vr_pn .* (vr_pn.^2 .+ vi_pn.^2).^(beta/2 .-1))
        # JuMP.@constraint(model, pd .==  vr_pn .* crd_vec[p] .+ vi_pn .* cid_vec[p])
        # JuMP.@constraint(model, qd .== -vr_pn .* cid_vec[p] .+ vi_pn .* crd_vec[p])
        # crd_bus[id] = JuMP.Containers.DenseAxisArray([crd_vec[p]..., -sum(crd_vec[p])], connections)
        # cid_bus[id] = JuMP.Containers.DenseAxisArray([cid_vec[p]..., -sum(cid_vec[p])], connections)

        # # crd_bus[id] = JuMP.Containers.DenseAxisArray([t in crd_bus[id].axes[1] ? crd_bus[id][t] : 0 for t in 1:n_ph], 1:n_ph)
        # # cid_bus[id] = JuMP.Containers.DenseAxisArray([t in cid_bus[id].axes[1] ? cid_bus[id][t] : 0 for t in 1:n_ph], 1:n_ph)
        
    else
        phases = connections
        phases_next = [connections[2:end]..., connections[1]]
        P = length(connections)
        idxs = 1:P
        idxs_prev = [idxs[end], idxs[1:end-1]...]
        
        # Md = PMD._get_delta_transformation_matrix(length(connections))
        # vrd = Md*[vr[p] for p in phases]
        # vid = Md*[vi[p] for p in phases]

        vrd = [vr[bus_id][c]-vr[bus_id][d] for (c,d) in zip(phases,phases_next)]
        vid = [vi[bus_id][c]-vi[bus_id][d] for (c,d) in zip(phases,phases_next)]

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

        JuMP.@NLconstraint(model, [i in 1:P], crd[id][i] == a[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1) + b[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1))
        JuMP.@NLconstraint(model, [i in 1:P], cid[id][i] == a[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1) - b[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1))
        JuMP.@constraint(model, [i in 1:P], pd[i] ==  vrd[i]*crd[id][i] + vid[i]*cid[id][i])
        JuMP.@constraint(model, [i in 1:P], qd[i] == -vrd[i]*cid[id][i] + vid[i]*crd[id][i])

        crd_bus_expr = JuMP.@NLexpression(model, [i in 1:P], crd[id][i]-crd[id][idxs_prev[i]])
        cid_bus_expr = JuMP.@NLexpression(model, [i in 1:P], cid[id][i]-cid[id][idxs_prev[i]])
        crd_bus[id] = _merge_bus_flows(model, crd_bus_expr, connections)
        cid_bus[id] = _merge_bus_flows(model, cid_bus_expr, connections)
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


# "gen connections adaptation of min fuel cost polynomial linquad objective"
JuMP.@objective(model, Min, sum(gen["cost"][1]*sum(pg[:,i]) + gen["cost"][2] for (i,gen) in ref[:gen]))


JuMP.optimize!(model)
cost = JuMP.objective_value(model)


# ##
# model_variables = JuMP.num_variables(model)

# # for consistency with other solvers, skip the variable bounds in the constraint count
# non_nl_constraints = sum(JuMP.num_constraints(model, ft, st) for (ft, st) in JuMP.list_of_constraint_types(model) if ft != JuMP.VariableRef)
# model_constraints = JuMP.num_nonlinear_constraints(model) + non_nl_constraints

# model_build_time = time() - time_model_start

# time_solve_start = time()

# JuMP.optimize!(model)
# cost = JuMP.objective_value(model)
# feasible = (JuMP.termination_status(model) == JuMP.LOCALLY_SOLVED)

# solve_time = time() - time_solve_start
# total_time = time() - time_data_start

##
PMD.add_start_vrvi!(data_math)
pm = PMD.instantiate_mc_model(data_math, PMD.IVRENPowerModel, PMD.build_mc_opf);
res = PMD.optimize_model!(pm, optimizer=ipopt_solver)
# println(pm.model)
