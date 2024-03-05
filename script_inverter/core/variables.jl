terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref[:bus])
v_start = [exp.(im.*collect(0:-1:-2)*2/3*pi) ; 0]

n_ph = 4

# """ GFL may not need vm variables """
vr = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vr", start = real(v_start)[t], lower_bound = -ref[:bus][i]["vmax"][t], upper_bound = ref[:bus][i]["vmax"][t] ) for i in keys(ref[:bus]))
vr = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vr[i].axes[1] ? vr[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
vi = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vi", start = imag(v_start)[t], lower_bound = -ref[:bus][i]["vmax"][t], upper_bound = ref[:bus][i]["vmax"][t] ) for i in keys(ref[:bus]))
vi = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vi[i].axes[1] ? vi[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
vm = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm", start = 1, lower_bound = 0 ) for i in keys(ref[:bus]))
vm = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm[i].axes[1] ? vm[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))

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
connections = Dict(i => gen["connections"][1:int_dim[i]] for (i, gen) in ref[:gen])

# crg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="crg_$i") for i in keys(ref[:gen]))
# cig = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cig_$i") for i in keys(ref[:gen]))
# crg = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? crg[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
# cig = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? cig[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
crg = Dict(i => JuMP.@variable(model, [c in connections[i]], base_name="crg_$i") for i in keys(ref[:gen]))
cig = Dict(i => JuMP.@variable(model, [c in connections[i]], base_name="cig_$i") for i in keys(ref[:gen]))
crg = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in crg[i].axes[1] ? crg[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
cig = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in cig[i].axes[1] ? cig[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))

# pg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="pg_$i") for i in keys(ref[:gen]))
# qg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="qg_$i") for i in keys(ref[:gen]))
# pg = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? pg[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
# qg = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? qg[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
pg = Dict(i => JuMP.@variable(model, [c in connections[i]], base_name="pg_$i") for i in keys(ref[:gen]))
qg = Dict(i => JuMP.@variable(model, [c in connections[i]], base_name="qg_$i") for i in keys(ref[:gen]))
pg = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in pg[i].axes[1] ? pg[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
qg = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in qg[i].axes[1] ? qg[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
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