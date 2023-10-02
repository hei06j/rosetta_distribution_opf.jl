using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import OpenDSSDirect
import InfrastructureModels
using Ipopt
using JuMP  # bl/array_nl
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels
const ODSS = OpenDSSDirect


##
ODSS.dss("Redirect ./data/test_gen_3ph_wye.dss")
vm_noControl = ODSS.Circuit.AllBusVolts()
cm_noControl = ODSS.PDElements.AllCurrentsAllCurrents()[5:8]
vm1_noControl = abs.(vm_noControl)[1:3]
vm2_noControl = abs.(vm_noControl)[4:7]
pg_noControl = ODSS.PVsystems.kW()
qg_noControl = ODSS.PVsystems.kvar()

ODSS.dss("Redirect ./data/test_gen_3w_wye_InvControl.dss")
vm_Control = ODSS.Circuit.AllBusVolts()
cm_Control = ODSS.PDElements.AllCurrentsAllCurrents()[5:8]
vm1_Control = abs.(vm_Control)[1:3]
vm2_Control = abs.(vm_Control)[4:7]
pg_Control = ODSS.PVsystems.kW()
qg_Control = ODSS.PVsystems.kvar()

# plot([vm1_noControl[1:3] vm2_noControl[1:3]]', marker=:circle, linestyle=:dash, color=[:blue :red :green], label=["v_a" "v_b" "v_c"])
# plot!([vm1_Control[1:3] vm2_Control[1:3]]', marker=:circle, linestyle=:solid, color=[:blue :red :green], label=["v_a^c" "v_b^c" "v_c^c"])

[pg_noControl pg_Control]
[qg_noControl qg_Control]
qg_Control
round.(vm2_Control./230.94, digits=4)

sum(qg_values[:,gen_id])
round.(abs.(JuMP.value.(vm))[:,1], digits=4)
# qmax = sqrt.((25)^2 - (23)^2)/3 * ones(3)
# # sum(-2*qmax/0.2 .* ((vm2_noControl[1:3] .-  vm2_noControl[4])/230.94.- 1))
# sum(-2*qmax/0.2 .* ((vm2_Control[1:3])/230.94.- 1))

##
data_path = "./data/inverter_3w_wye_unbalanced_loads.dss"

ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
data_eng = PMD.parse_file(data_path, transformations=[PMD.remove_all_bounds!, PMD.transform_loops!])
RPMD.pv1_correction!(data_eng)
data_eng["settings"]["sbase_default"] = 1
vbase = data_eng["settings"]["vbases_default"]["b1"]
data_eng["voltage_source"]["source"]["rs"] *= 0
data_eng["voltage_source"]["source"]["xs"] *= 0
data_math = PMD.transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)

data_math["bus"]["1"]["vmin"] = [0.9 * ones(3) ; 0 ]
data_math["bus"]["1"]["vmax"] = [1.1 * ones(3) ; Inf]

gen_id = 1
gen = data_math["gen"]["$gen_id"]
smax = 40
pmax = 35
gen["pmax"] = pmax/3 * ones(3)
gen["pmin"] = pmax/3 * ones(3)
gen["pg"] = pmax/3 * ones(3)
gen["qmax"] = sqrt.(smax^2 - pmax^2)/3 * ones(3)
gen["qmin"] = -gen["qmax"]

gen["connections"] = gen["connections"][1:3]
data_math["bus"]["1"]["terminals"] = data_math["bus"]["1"]["terminals"][1:4]
data_math["bus"]["1"]["grounded"] = data_math["bus"]["1"]["grounded"][1:4]

gen["cost"] = [0 0]
data_math["gen"]["2"]["cost"] = [1000 0]
##
ref = IM.build_ref(data_math, PMD.ref_add_core!, PMD._pmd_global_keys, PMD.pmd_it_name)[:it][:pmd][:nw][0]

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

# int_dim = Dict(i => RPMD._infer_int_dim_unit(gen, false) for (i,gen) in ref[:gen])
int_dim = Dict(i => RPMD._infer_int_dim_unit(gen, !(4 in gen["connections"])) for (i,gen) in ref[:gen])
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
    ### if neutral is in generator connections, then false, otherwise it is true
    explicit_neutral = 4 in generator["connections"]
    nphases = RPMD._infer_int_dim_unit(generator, !explicit_neutral)
    # nphases = RPMD._infer_int_dim_unit(generator, false)
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
        if explicit_neutral
            phases = connections[1:end-1]
            n = connections[end]
        else
            phases = connections
            n = 4
        end
            crg_bus[id] = JuMP.Containers.DenseAxisArray([crg[phases,id]..., -sum(crg[phases,id])], connections)
            cig_bus[id] = JuMP.Containers.DenseAxisArray([cig[phases,id]..., -sum(cig[phases,id])], connections)
            # JuMP.@constraint(model,  pg[phases,id] .==  (vr[phases,bus_id] .- vr[n,bus_id]) .* crg[phases,id] .+ (vi[phases,bus_id] .- vi[n,bus_id]) .* cig[phases,id])
            JuMP.@constraint(model,  pg[phases,id] .==  vr[phases,bus_id] .* crg[phases,id] .+ vi[phases,bus_id] .* cig[phases,id])
            
            nonInf_pmin_pmax = [idx for (idx,t) in enumerate(pmin) if pmin[idx].>-Inf || pmax[idx].<Inf]
            JuMP.@constraint(model, pmin[nonInf_pmin_pmax] .<= pg[nonInf_pmin_pmax,id] .<= pmax[nonInf_pmin_pmax] )
    
            nonInf_qmin_qmax = [idx for (idx,t) in enumerate(qmin) if qmin[idx].>-Inf || qmax[idx].<Inf]
            # JuMP.@constraint(model, qg[nonInf_qmin_qmax,id] .== -(vr[nonInf_qmin_qmax,bus_id] .- vr[n,bus_id]) .* cig[nonInf_qmin_qmax,id] .+ (vi[nonInf_qmin_qmax,bus_id] .- vi[n,bus_id]) .* crg[nonInf_qmin_qmax,id])
            JuMP.@constraint(model, qg[nonInf_qmin_qmax,id] .== -vr[nonInf_qmin_qmax,bus_id] .* cig[nonInf_qmin_qmax,id] .+ vi[nonInf_qmin_qmax,bus_id] .* crg[nonInf_qmin_qmax,id])
            JuMP.@constraint(model, qmin[nonInf_qmin_qmax] .<= qg[nonInf_qmin_qmax,id] .<= qmax[nonInf_qmin_qmax] )
    
        # else
        #     phases = connections
        #     crg_bus[id] = crg[:,id]
        #     cig_bus[id] = cig[:,id]
        #     JuMP.@constraint(model,  pg[phases,id] .==  vr[phases,bus_id] .* crg[phases,id] .+ vi[phases,bus_id] .* cig[phases,id])

        #     nonInf_pmin_pmax = [idx for (idx,t) in enumerate(pmin) if pmin[idx].>-Inf || pmax[idx].<Inf]
        #     JuMP.@constraint(model, pmin[nonInf_pmin_pmax] .<= pg[nonInf_pmin_pmax,id] .<= pmax[nonInf_pmin_pmax] )
        #     # JuMP.@constraint(model, pmin .<= pg[phases,id] .<= pmax )

        #     nonInf_qmin_qmax = [idx for (idx,t) in enumerate(qmin) if qmin[idx].>-Inf || qmax[idx].<Inf]
        #     JuMP.@constraint(model, qg[nonInf_qmin_qmax,id] .== -vr[nonInf_qmin_qmax,bus_id] .* cig[nonInf_qmin_qmax,id] .+ vi[nonInf_qmin_qmax,bus_id] .* crg[nonInf_qmin_qmax,id])
        #     JuMP.@constraint(model, qmin[nonInf_qmin_qmax] .<= qg[nonInf_qmin_qmax,id] .<= qmax[nonInf_qmin_qmax] )
        #     # JuMP.@constraint(model, qmin .<= qg[phases,id] .<= qmax )
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

# "gen connections adaptation of min fuel cost polynomial linquad objective"
JuMP.@objective(model, Min, sum(gen["cost"][1]*sum(pg[:,i]) + gen["cost"][2] for (i,gen) in ref[:gen]))


## Inverter control Volt-var
gen_id = 1
gen = data_math["gen"]["$gen_id"]
bus_id = gen["gen_bus"]
vmin = 0.9; vmax = 1.1;
phases = gen["connections"]
# n = gen["connections"][end]

terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref[:bus])
vm = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm", lower_bound = 0 ) for i in keys(ref[:bus]))
vm = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm[i].axes[1] ? vm[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
for (i, bus) in ref[:bus]
    JuMP.@constraint(model, [p in ref[:bus][i]["terminals"]], vm[p,i]^2 == vr[p,i]^2 + vi[p,i]^2)
end

vv_curve(v, qmax) = v <= 0.95 ? qmax : (v >= 1.05 ? -qmax : -2*qmax/0.1*(v-1))
JuMP.@operator(model, vv, 2, vv_curve)

# JuMP.@constraint(model, [i in phases], qg[i,gen_id] == vv(vm_pv[i,bus_id], gen["qmax"][i]) )
JuMP.@constraint(model, qg[1, gen_id] == vv(sum(vm[i,bus_id] for i in ref[:bus][1]["terminals"][1:3])/3, gen["qmax"][1]) )
JuMP.@constraint(model, qg[2, gen_id] == qg[1, gen_id])
JuMP.@constraint(model, qg[3, gen_id] == qg[1, gen_id])

##
JuMP.optimize!(model)
cost = JuMP.objective_value(model)
v = abs.(JuMP.value.(vm)) * 230.94
pg_values = JuMP.value.(pg)
qg_values = JuMP.value.(qg)
vbase
sum(pg_values[:,gen_id])
sum(qg_values[:,gen_id])
round.(abs.(JuMP.value.(vm))[:,1], digits=4)

# vm_range = 0.8:0.01:1.2
# plot(vm_range, vv_curve.(vm_range,gen["qmax"][1]), label="VV curve")
# # plot!(vm_range, vv_curve.(vm_range,gen["qmax"][2]))
# # plot!(vm_range, vv_curve.(vm_range,gen["qmax"][3]))
# plot!([sum(v[phases,bus_id])/vbase/1000/3], [qg_values[1,gen_id]], seriestype=:scatter, label="VV point")

# vm_range = 0.8:0.01:1.2
# vw_curve(v, pmax) = v <= 0.95 ? pmax : (v >= 1.05 ? 0 : -pmax/0.1*(v-1.05))
# plot(vm_range, vw_curve.(vm_range,3), label="VW curve")
