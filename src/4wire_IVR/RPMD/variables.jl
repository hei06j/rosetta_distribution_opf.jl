terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref[:bus])
v_start = [exp.(im.*collect(0:-1:-2)*2/3*pi) ; 0]
vr = Dict(i => JuMP.@variable(model, [t in terminals[i]], base_name="vr_$(i)", start = real(v_start)[t], lower_bound = -ref[:bus][i]["vmax"][t], upper_bound = ref[:bus][i]["vmax"][t] ) for i in keys(ref[:bus]))
vi = Dict(i => JuMP.@variable(model, [t in terminals[i]], base_name="vi_$(i)", start = imag(v_start)[t], lower_bound = -ref[:bus][i]["vmax"][t], upper_bound = ref[:bus][i]["vmax"][t] ) for i in keys(ref[:bus]))
# perfectly grounded terminals get a constant 0 instead of a variable
vr = Dict(i=>JuMP.Containers.DenseAxisArray(Vector{JuMP.AffExpr}([t in vr[i].axes[1] ? vr[i][t] : 0.0 for t in bus["terminals"]]), bus["terminals"]) for (i, bus) in ref[:bus])
vi = Dict(i=>JuMP.Containers.DenseAxisArray(Vector{JuMP.AffExpr}([t in vi[i].axes[1] ? vi[i][t] : 0.0 for t in bus["terminals"]]), bus["terminals"]) for (i, bus) in ref[:bus])

nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref[:branch])
cr = Dict((l,i,j) => JuMP.@variable(model, [c in 1:nconds[l]], base_name="cr_$((l,i,j))", start = _PMD.comp_start_value(ref[:branch][l], "cr_start", c, 0.0) ) for (l,i,j) in ref[:arcs_branch]) # , lower_bound = -ref[:branch][l]["c_rating_a"], upper_bound = ref[:branch][l]["c_rating_a"]
ci = Dict((l,i,j) => JuMP.@variable(model, [c in 1:nconds[l]], base_name="ci_$((l,i,j))", start = _PMD.comp_start_value(ref[:branch][l], "ci_start", c, 0.0) ) for (l,i,j) in ref[:arcs_branch]) 
csr = Dict(l => JuMP.@variable(model, [c in 1:nconds[l]], base_name="csr_$(l)", start = _PMD.comp_start_value(ref[:branch][l], "csr_start", c, 0.0)) for (l,i,j) in ref[:arcs_branch])
csi = Dict(l => JuMP.@variable(model, [c in 1:nconds[l]], base_name="csi_$(l)", start = _PMD.comp_start_value(ref[:branch][l], "csi_start", c, 0.0)) for (l,i,j) in ref[:arcs_branch])
cr_bus = Dict{Tuple{Int,Int,Int}, Any}()
ci_bus = Dict{Tuple{Int,Int,Int}, Any}()

int_dim = Dict(i => _RPMD._infer_int_dim_unit(gen, false) for (i,gen) in ref[:gen])
crg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="crg_$(i)", start = _PMD.comp_start_value(ref[:gen][i], "crg_start", c, 0.0)) for i in keys(ref[:gen]))
cig = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cig_$(i)", start = _PMD.comp_start_value(ref[:gen][i], "cig_start", c, 0.0)) for i in keys(ref[:gen]))
pg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="pg_$(i)", start = _PMD.comp_start_value(ref[:gen][i], "pg_start", c, 0.0)) for i in keys(ref[:gen]))
qg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="qg_$(i)", start = _PMD.comp_start_value(ref[:gen][i], "qg_start", c, 0.0)) for i in keys(ref[:gen]))
crg_bus = Dict{Int, Any}()
cig_bus = Dict{Int, Any}()
# # store active and reactive power expressions for use in objective + post processing
# var(pm, nw)[:pg] = Dict{Int, Any}()
# var(pm, nw)[:qg] = Dict{Int, Any}()

int_dim = Dict(l => _RPMD._infer_int_dim_transformer(trans, false) for (l,trans) in ref[:transformer])
crt = Dict((l,i,j) => JuMP.@variable(model, [c in 1:int_dim[l]], base_name="crt_$((l,i,j))", start = _PMD.comp_start_value(ref[:transformer][l], "crt_start", c, 0.0) ) for (l,i,j) in ref[:arcs_transformer])
cit = Dict((l,i,j) => JuMP.@variable(model, [c in 1:int_dim[l]], base_name="cit_$((l,i,j))", start = _PMD.comp_start_value(ref[:transformer][l], "cit_start", c, 0.0) ) for (l,i,j) in ref[:arcs_transformer])
pt = Dict((l,i,j) => JuMP.@variable(model, [c in 1:int_dim[l]], base_name="pt_$((l,i,j))", start = _PMD.comp_start_value(ref[:transformer][l], "pt_start", c, 0.0) ) for (l,i,j) in ref[:arcs_transformer])
qt = Dict((l,i,j) => JuMP.@variable(model, [c in 1:int_dim[l]], base_name="qt_$((l,i,j))", start = _PMD.comp_start_value(ref[:transformer][l], "qt_start", c, 0.0) ) for (l,i,j) in ref[:arcs_transformer])
crt_bus = Dict{Tuple{Int,Int,Int}, Any}()
cit_bus = Dict{Tuple{Int,Int,Int}, Any}()

nconds = Dict(l => length(switch["f_connections"]) for (l,switch) in ref[:switch])
crsw = Dict((l,i,j) => JuMP.@variable(model, [c in 1:nconds[l]], base_name="crsw_$((l,i,j))", start = _PMD.comp_start_value(ref[:switch][l], "crsw_start", c, 0.0)) for (l,i,j) in ref[:arcs_switch]) # , lower_bound = -ref[:switch][l]["current_rating"], upper_bound = ref[:switch][l]["current_rating"]
cisw = Dict((l,i,j) => JuMP.@variable(model, [c in 1:nconds[l]], base_name="cisw_$((l,i,j))", start = _PMD.comp_start_value(ref[:switch][l], "cisw_start", c, 0.0)) for (l,i,j) in ref[:arcs_switch]) # , lower_bound = -ref[:switch][l]["current_rating"], upper_bound = ref[:switch][l]["current_rating"]
crsw_bus = Dict{Tuple{Int,Int,Int}, Any}()
cisw_bus = Dict{Tuple{Int,Int,Int}, Any}()

int_dim = Dict(i => _RPMD._infer_int_dim_unit(load, false) for (i,load) in ref[:load])
crd = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="crd_$(i)") for i in keys(ref[:load]))
cid = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cid_$(i)") for i in keys(ref[:load]))
# pd = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="pd_$(i)") for i in keys(ref[:load]))
# qd = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="qd_$(i)") for i in keys(ref[:load]))
crd_bus = Dict{Int, Any}()
cid_bus = Dict{Int, Any}()
# var(pm, nw)[:pd] = Dict{Int, Any}()
# var(pm, nw)[:qd] = Dict{Int, Any}()
# var(pm, nw)[:pd_bus] = Dict{Int, Any}()
# var(pm, nw)[:qd_bus] = Dict{Int, Any}()