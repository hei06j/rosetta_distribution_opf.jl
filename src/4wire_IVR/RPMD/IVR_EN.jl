"The solar parsing is a bit off, this method corrects for that."
function pv1_correction!(data_eng)
    if haskey(data_eng, "solar")
        pv1 = data_eng["solar"]["pv1"]
        pv1["pg_lb"] = pv1["pg_ub"] = pv1["pg"]
        pv1["qg_lb"] = pv1["qg_ub"] = pv1["qg"]
    end
end

function IVR_EN(model, ref)
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

    int_dim = Dict(i => _infer_int_dim_unit(gen, false) for (i,gen) in ref[:gen])
    crg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="crg_$(i)", start = _PMD.comp_start_value(ref[:gen][i], "crg_start", c, 0.0)) for i in keys(ref[:gen]))
    cig = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cig_$(i)", start = _PMD.comp_start_value(ref[:gen][i], "cig_start", c, 0.0)) for i in keys(ref[:gen]))
    pg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="pg_$(i)", start = _PMD.comp_start_value(ref[:gen][i], "pg_start", c, 0.0)) for i in keys(ref[:gen]))
    qg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="qg_$(i)", start = _PMD.comp_start_value(ref[:gen][i], "qg_start", c, 0.0)) for i in keys(ref[:gen]))
    crg_bus = Dict{Int, Any}()
    cig_bus = Dict{Int, Any}()
    # # store active and reactive power expressions for use in objective + post processing
    # var(pm, nw)[:pg] = Dict{Int, Any}()
    # var(pm, nw)[:qg] = Dict{Int, Any}()

    int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in ref[:transformer])
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
    
    int_dim = Dict(i => _infer_int_dim_unit(load, false) for (i,load) in ref[:load])
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

    # to_variable(x::AffExpr) = first(keys(x.terms))
    # v_start = exp.(im.*collect(0:-1:-2)*2/3*pi)
    # for i in keys(ref[:bus])
    #     for j in 1:3
    #         JuMP.set_start_value(to_variable(vr[i][j]), real(v_start[j]))
    #         JuMP.set_start_value(to_variable(vi[i][j]), imag(v_start[j]))
    #     end
    #     if vr[i][4] != 0
    #         JuMP.set_start_value(to_variable(vr[i][4]), 0)
    #         JuMP.set_start_value(to_variable(vi[i][4]), 0)
    #     end
    # end

    for (i, bus) in ref[:bus]
        terminals = bus["terminals"]
        grounded = bus["grounded"]

        if i in keys(ref[:ref_buses])
            # if haskey(bus, "va") && !haskey(bus, "vm")
            #     #constraint_mc_theta_ref(pm, nw, id, bus["va"], terminals, grounded)
            # elseif haskey(bus, "vm") && !haskey(bus, "va")
            #     #constraint_mc_voltage_magnitude_fixed(pm, nw, id, bus["vm"], terminals, grounded)
            # elseif haskey(bus, "vm") && haskey(bus, "va")
                idxs = findall((!).(grounded))
                JuMP.@constraint(model, [k in idxs], vr[i][terminals[k]] == bus["vm"][k]*cos(bus["va"][k]))
                JuMP.@constraint(model, [k in idxs], vi[i][terminals[k]] == bus["vm"][k]*sin(bus["va"][k]))
            # end
        end

        ungrounded_terminals = terminals[(!).(grounded)]
        for (idx,t) in enumerate(terminals)
            if !grounded[idx]
                if bus["vmax"][idx] < Inf
                    JuMP.@constraint(model, vr[i][t]^2+vi[i][t]^2 <= bus["vmax"][idx]^2)
                end
                if bus["vmin"][idx] > 0.0
                    JuMP.@constraint(model, vr[i][t]^2+vi[i][t]^2 >= bus["vmin"][idx]^2)
                end
            end
        end

        # for (a,b,lb) in bus["vm_pair_lb"]
        #     if lb > 0.0
        #         JuMP.@constraint(model, (vr[i][a]-vr[i][b])^2 + (vi[i][a]-vi[i][b])^2 >= lb^2)
        #     end
        # end

        # for (a,b,ub) in bus["vm_pair_ub"]
        #     if ub < Inf
        #         JuMP.@constraint(model, (vr[i][a]-vr[i][b])^2 + (vi[i][a]-vi[i][b])^2 <= ub^2)
        #     end
        # end
        
    end


    for (id, generator) in ref[:gen]
        nphases = _infer_int_dim_unit(generator, false)
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
        if configuration==_PMD.WYE || length(pmin)==1 || nphases==1
            phases = connections[1:end-1]
            n      = connections[end]

            crg_bus[id] = _merge_bus_flows(model, [crg[id]..., -sum(crg[id])], connections)
            cig_bus[id] = _merge_bus_flows(model, [cig[id]..., -sum(cig[id])], connections)

            for (idx, p) in enumerate(phases)
                JuMP.@NLconstraint(model,  pg[id][idx] ==  (vr[bus_id][p]-vr[bus_id][n])*crg[id][idx]+(vi[bus_id][p]-vi[bus_id][n])*cig[id][idx])
                JuMP.@NLconstraint(model,  qg[id][idx] == -(vr[bus_id][p]-vr[bus_id][n])*cig[id][idx]+(vi[bus_id][p]-vi[bus_id][n])*crg[id][idx])
            end
            for (idx, p) in enumerate(phases)
                if pmin[idx]>-Inf
                    JuMP.@constraint(model, pmin[idx] .<= (vr[bus_id][p]-vr[bus_id][n])*crg[id][idx]  + (vi[bus_id][p]-vi[bus_id][n])*cig[id][idx])
                end
                if pmax[idx]< Inf
                    JuMP.@constraint(model, pmax[idx] .>= (vr[bus_id][p]-vr[bus_id][n])*crg[id][idx]  + (vi[bus_id][p]-vi[bus_id][n])*cig[id][idx])
                end
                if qmin[idx]>-Inf
                    JuMP.@constraint(model, qmin[idx] .<= (vi[bus_id][p]-vi[bus_id][n])*crg[id][idx]  - (vr[bus_id][p]-vr[bus_id][n])*cig[id][idx])
                end
                if qmax[idx]< Inf
                    JuMP.@constraint(model, qmax[idx] .>= (vi[bus_id][p]-vi[bus_id][n])*crg[id][idx]  - (vr[bus_id][p]-vr[bus_id][n])*cig[id][idx])
                end
            end

        else ## configuration==_PMD.DELTA

            Md =_PMD._get_delta_transformation_matrix(length(connections))
            crg_bus[id] = _merge_bus_flows(model, Md'*crg[id], connections)
            cig_bus[id] = _merge_bus_flows(model, Md'*cig[id], connections)
            _merge_bus_flows(model, Md'*crg[id], connections)
            _merge_bus_flows(model, Md'*cig[id], connections)

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

            # pg[id] = JuMP.Containers.DenseAxisArray(pg_expr, connections)
            # qg[id] = JuMP.Containers.DenseAxisArray(qg_expr, connections)
        end
    end


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

        vr_fr = [vr[f_bus][c] for c in f_connections]
        vi_fr = [vi[f_bus][c] for c in f_connections]
        vr_to = [vr[t_bus][c] for c in t_connections]
        vi_to = [vi[t_bus][c] for c in t_connections]

        cr_fr = cr[f_idx]
        ci_fr = ci[f_idx]
        cr_to = cr[t_idx]
        ci_to = ci[t_idx]

        csr_fr = csr[f_idx[1]]
        csi_fr = csi[f_idx[1]]
        csr_to = -csr[t_idx[1]]
        csi_to = -csi[t_idx[1]]

        ### constraint_mc_current_from
        JuMP.@constraint(model, cr_fr .== csr_fr + g_fr*vr_fr - b_fr*vi_fr)
        JuMP.@constraint(model, ci_fr .== csi_fr + g_fr*vi_fr + b_fr*vr_fr)
        cr_bus[f_idx] = _merge_bus_flows(model, cr_fr, f_connections)
        ci_bus[f_idx] = _merge_bus_flows(model, ci_fr, f_connections)

        ### constraint_mc_current_to
        JuMP.@constraint(model, cr_to .== csr_to + g_to*vr_to - b_to*vi_to)
        JuMP.@constraint(model, ci_to .== csi_to + g_to*vi_to + b_to*vr_to)
        cr_bus[t_idx] = _merge_bus_flows(model, cr_to, t_connections)
        ci_bus[t_idx] = _merge_bus_flows(model, ci_to, t_connections)

        ### constraint_mc_bus_voltage_drop
        JuMP.@constraint(model, vr_to .== vr_fr - r*csr_fr + x*csi_fr)
        JuMP.@constraint(model, vi_to .== vi_fr - r*csi_fr - x*csr_fr)

        ### constraint_mc_branch_current_limit
        cnds_finite_rating = [c for (c,r) in enumerate(c_rating) if r<Inf]
        JuMP.@constraint(model, [c in cnds_finite_rating], cr_fr[c]^2+ci_fr[c]^2 <= c_rating[c]^2)
        JuMP.@constraint(model, [c in cnds_finite_rating], cr_to[c]^2+ci_to[c]^2 <= c_rating[c]^2)

        ### constraint_mc_thermal_limit
        if haskey(branch, "rate_a") && any(branch["rate_a"] .< Inf)
            for idx in 1:length(rate_a)
                if branch["rate_a"][idx]<Inf
                    ### constraint_mc_thermal_limit_from
                    pf_idx = JuMP.@NLexpression(model,  vr_fr[idx]*cr_fr[idx] + vi_fr[idx]*ci_fr[idx])
                    qf_idx = JuMP.@NLexpression(model, -vr_fr[idx]*ci_fr[idx] + vi_fr[idx]*cr_fr[idx])
                    JuMP.@NLconstraint(model, pf_idx^2 + qf_idx^2 <= branch["rate_a"][idx]^2)

                    ### constraint_mc_thermal_limit_to
                    pt_idx = JuMP.@NLexpression(model,  vr_to[idx]*cr_to[idx] + vi_to[idx]*ci_to[idx])
                    qt_idx = JuMP.@NLexpression(model, -vr_to[idx]*ci_to[idx] + vi_to[idx]*cr_to[idx])
                    JuMP.@NLconstraint(model, pt_idx^2 + qt_idx^2 <= branch["rate_a"][idx]^2)
                end
            end
        end

    end


    for (id, load) in ref[:load]
        bus_id = load["load_bus"]
        bus = ref[:bus][bus_id]
        configuration = load["configuration"]
        connections = load["connections"]
        load_model = load["model"]
        a, alpha, b, beta = _PMD._load_expmodel_params(load, bus)

        int_dim = _infer_int_dim_unit(load, false)
        if configuration==_PMD.WYE || int_dim==1
            phases = connections[1:end-1]
            n      = connections[end]

            vr_pn = [vr[bus_id][p]-vr[bus_id][n] for p in phases]
            vi_pn = [vi[bus_id][p]-vi[bus_id][n] for p in phases]

            if load_model==_PMD.POWER
                pd = a
                qd = b
            elseif load_model==_PMD.IMPEDANCE
                pd = a .* (vr_pn.^2 .+ vi_pn.^2)
                qd = b .* (vr_pn.^2 .+ vi_pn.^2)
            elseif load_model==_PMD.CURRENT
                pd = JuMP.@variable(model, [c in 1:int_dim])
                qd = JuMP.@variable(model, [c in 1:int_dim])
                JuMP.@constraint(model, pd.^2 .== a.^2 .* (vr_pn.^2 .+ vi_pn.^2))
                JuMP.@constraint(model, sign.(a).*pd .>= 0)
                JuMP.@constraint(model, qd.^2 .== b.^2 .* (vr_pn.^2 .+ vi_pn.^2))
                JuMP.@constraint(model, sign.(b).*qd .>= 0)
            else
                error("Load model $model for load $id is not supported by this formulation.")
            end


            for (idx, p) in enumerate(phases)
                JuMP.@NLconstraint(model, crd[id][idx] == 
                    a[idx] * vr_pn[idx] * (vr_pn[idx]^2+vi_pn[idx]^2)^(alpha[idx]/2-1)
                    +b[idx] * vi_pn[idx] * (vr_pn[idx]^2+vi_pn[idx]^2)^(beta[idx]/2 -1))
                JuMP.@NLconstraint(model, cid[id][idx] == 
                    a[idx] * vi_pn[idx] * (vr_pn[idx]^2+vi_pn[idx]^2)^(alpha[idx]/2-1)
                    -b[idx] * vr_pn[idx] * (vr_pn[idx]^2+vi_pn[idx]^2)^(beta[idx]/2 -1))

                JuMP.@NLconstraint(model, pd[idx] ==  vr_pn[idx]*crd[id][idx] + vi_pn[idx]*cid[id][idx])
                JuMP.@NLconstraint(model, qd[idx] == -vr_pn[idx]*cid[id][idx] + vi_pn[idx]*crd[id][idx])
            end
            crd_bus_n = JuMP.@NLexpression(model, -sum(crd[id][i] for i in 1:length(phases)))
            cid_bus_n = JuMP.@NLexpression(model, -sum(cid[id][i] for i in 1:length(phases)))
            crd_bus[id] = _merge_bus_flows(model, [crd[id]..., crd_bus_n], connections)
            cid_bus[id] = _merge_bus_flows(model, [cid[id]..., cid_bus_n], connections)


            
        else
            phases = connections
            phases_next = [connections[2:end]..., connections[1]]
            P = length(connections)
            idxs = 1:P
            idxs_prev = [idxs[end], idxs[1:end-1]...]
            
            # Md = _PMD._get_delta_transformation_matrix(length(connections))
            # vrd = Md*[vr[p] for p in phases]
            # vid = Md*[vi[p] for p in phases]

            vrd = [vr[bus_id][c]-vr[bus_id][d] for (c,d) in zip(phases,phases_next)]
            vid = [vi[bus_id][c]-vi[bus_id][d] for (c,d) in zip(phases,phases_next)]

            if load_model==_PMD.POWER
                pd = a
                qd = b
            elseif load_model==_PMD.IMPEDANCE
                pd = a .* (vrd.^2 .+ vid.^2)
                qd = b .* (vrd.^2 .+ vid.^2)
            elseif load_model==_PMD.CURRENT
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


    fix_taps = true
    for (i, transformer) in ref[:transformer]
        f_bus = transformer["f_bus"]
        t_bus = transformer["t_bus"]
        f_idx = (i, f_bus, t_bus)
        t_idx = (i, t_bus, f_bus)
        configuration = transformer["configuration"]
        f_connections = transformer["f_connections"]
        t_connections = transformer["t_connections"]
        tm_set = transformer["tm_set"]
        tm_fixed = fix_taps ? ones(Bool, length(tm_set)) : transformer["tm_fix"]
        tm_scale = _PMD.calculate_tm_scale(transformer, ref[:bus][f_bus], ref[:bus][t_bus])
        pol = transformer["polarity"]
        sm_ub = transformer["sm_ub"]

        vr_fr = vr[f_bus]
        vi_fr = vi[f_bus]
        vr_to = vr[t_bus]
        vi_to = vi[t_bus]

        crt_fr = crt[f_idx]
        cit_fr = cit[f_idx]
        crt_to = crt[t_idx]
        cit_to = cit[t_idx]

        if configuration == _PMD.WYE  || length(crt_fr)==1
            ### constraint_mc_transformer_voltage_yy
            vr_fr_P = [vr_fr[c] for c in f_connections[1:end-1]]
            vi_fr_P = [vi_fr[c] for c in f_connections[1:end-1]]
            vr_fr_n = vr_fr[f_connections[end]]
            vi_fr_n = vi_fr[f_connections[end]]
            vr_to_P = [vr_to[c] for c in t_connections[1:end-1]]
            vi_to_P = [vi_to[c] for c in t_connections[1:end-1]]
            vr_to_n = vr_to[t_connections[end]]
            vi_to_n = vi_to[t_connections[end]]
            # construct tm as a parameter or scaled variable depending on whether it is fixed or not
            tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
            scale = (tm_scale*pol).*tm_set
            JuMP.@constraint(model, (vr_fr_P.-vr_fr_n) .== scale.*(vr_to_P.-vr_to_n))
            JuMP.@constraint(model, (vi_fr_P.-vi_fr_n) .== scale.*(vi_to_P.-vi_to_n))

            ### constraint_mc_transformer_current_yy
            # construct tm as a parameter or scaled variable depending on whether it is fixed or not
            tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
            scale = (tm_scale*pol).*tm_set

            JuMP.@constraint(model, scale.*crt_fr .+ crt_to .== 0)
            JuMP.@constraint(model, scale.*cit_fr .+ cit_to .== 0)
            crt_bus[f_idx] = _merge_bus_flows(model, [crt_fr..., -sum(crt_fr)], f_connections)
            cit_bus[f_idx] = _merge_bus_flows(model, [cit_fr..., -sum(cit_fr)], f_connections)
            crt_bus[t_idx] = _merge_bus_flows(model, [crt_to..., -sum(crt_to)], t_connections)
            cit_bus[t_idx] = _merge_bus_flows(model, [cit_to..., -sum(cit_to)], t_connections)

        elseif configuration == _PMD.DELTA
            ### constraint_mc_transformer_voltage_dy
            vr_fr_P = [vr_fr[c] for c in f_connections]
            vi_fr_P = [vi_fr[c] for c in f_connections]
            vr_to_P = [vr_to[c] for c in t_connections[1:end-1]]
            vi_to_P = [vi_to[c] for c in t_connections[1:end-1]]
            vr_to_n = vr_to[t_connections[end]]
            vi_to_n = vi_to[t_connections[end]]
            # construct tm as a parameter or scaled variable depending on whether it is fixed or not
            tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
            scale = (tm_scale*pol).*tm_set
            n_phases = length(tm)
            Md = _PMD._get_delta_transformation_matrix(n_phases)
            JuMP.@constraint(model, Md*vr_fr_P .== scale.*(vr_to_P .- vr_to_n))
            JuMP.@constraint(model, Md*vi_fr_P .== scale.*(vi_to_P .- vi_to_n))

            ### constraint_mc_transformer_current_dy
            # construct tm as a parameter or scaled variable depending on whether it is fixed or not
            tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
            scale = (tm_scale*pol).*tm_set
            n_phases = length(tm)
            Md = _PMD._get_delta_transformation_matrix(n_phases)
            JuMP.@constraint(model, scale.*crt_fr .+ crt_to .== 0)
            JuMP.@constraint(model, scale.*cit_fr .+ cit_to .== 0)
            crt_bus[f_idx] = _merge_bus_flows(model, Md'*crt_fr, f_connections)
            cit_bus[f_idx] = _merge_bus_flows(model, Md'*cit_fr, f_connections)
            crt_bus[t_idx] = _merge_bus_flows(model, [crt_to..., -sum(crt_to)], t_connections)
            cit_bus[t_idx] = _merge_bus_flows(model, [cit_to..., -sum(cit_to)], t_connections)

        elseif configuration == "zig-zag"
            error("Zig-zag not yet supported.")
        end

        ### constraint_mc_transformer_thermal_limit
        if configuration==_PMD.WYE || length(crt_fr)==1
            P_fr = f_connections[1:end-1]
            n_fr = f_connections[end]
            vrt_fr = [vr_fr[p]-vr_fr[n_fr] for p in P_fr]
            vit_fr = [vi_fr[p]-vi_fr[n_fr] for p in P_fr]
        elseif configuration==_PMD.DELTA && length(crt_fr)==3
            M = _PMD._get_delta_transformation_matrix(3)
            vrt_fr = M*[vr_to[p] for p in f_connections]
            vit_fr = M*[vi_to[p] for p in f_connections]
        else
            error("The configuration $config of dimension $(length(crt)) is not supported.")
        end
        P_to = t_connections[1:end-1]
        n_to = t_connections[end]
        vrt_to = [vr_to[p]-vr_to[n_to] for p in P_to]
        vit_to = [vi_to[p]-vi_to[n_to] for p in P_to]

        idxs = [1:length(vrt_fr)...]
        pt_fr = JuMP.@NLexpression(model, [i in idxs],  vrt_fr[i]*crt_fr[i] + vit_fr[i]*cit_fr[i])
        qt_fr = JuMP.@NLexpression(model, [i in idxs], -vrt_fr[i]*cit_fr[i] + vit_fr[i]*crt_fr[i])
        pt_to = JuMP.@NLexpression(model, [i in idxs],  vrt_to[i]*crt_to[i] + vit_to[i]*cit_to[i])
        qt_to = JuMP.@NLexpression(model, [i in idxs], -vrt_to[i]*cit_to[i] + vit_to[i]*crt_to[i])

        if sm_ub<Inf
            JuMP.@NLconstraint(model, sum(pt_fr[i] for i in idxs)^2 + sum(qt_fr[i] for i in idxs)^2 <= sm_ub^2)
            JuMP.@NLconstraint(model, sum(pt_to[i] for i in idxs)^2 + sum(qt_to[i] for i in idxs)^2 <= sm_ub^2)
        end

    end


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

        cr    = cr_bus #get(cr_bus,   Dict()); _check_var_keys(cr, bus_arcs, "real current", "branch")
        ci    = ci_bus #get(ci_bus,   Dict()); _check_var_keys(ci, bus_arcs, "imaginary current", "branch")
        crd   = crd_bus #get(crd_bus,  Dict()); _check_var_keys(crd, bus_loads, "real current", "load")
        cid   = cid_bus #get(cid_bus,  Dict()); _check_var_keys(cid, bus_loads, "imaginary current", "load")
        crg   = crg_bus #get(crg_bus,  Dict()); _check_var_keys(crg, bus_gens, "real current", "generator")
        cig   = cig_bus #get(cig_bus,  Dict()); _check_var_keys(cig, bus_gens, "imaginary current", "generator")
        # crs   = get(crs_bus,  Dict()); _check_var_keys(crs, bus_storage, "real currentr", "storage")
        # cis   = get(cis_bus,  Dict()); _check_var_keys(cis, bus_storage, "imaginary current", "storage")
        # crsw  = get(crsw_bus, Dict()); _check_var_keys(crsw, bus_arcs_sw, "real current", "switch")
        # cisw  = get(cisw_bus, Dict()); _check_var_keys(cisw, bus_arcs_sw, "imaginary current", "switch")
        crt   = crt_bus #get(crt_bus,  Dict()); _check_var_keys(crt, bus_arcs_trans, "real current", "transformer")
        cit   = cit_bus #get(cit_bus,  Dict()); _check_var_keys(cit, bus_arcs_trans, "imaginary current", "transformer")

        Gt, Bt = _build_bus_shunt_matrices(ref, terminals, bus_shunts)

        ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

        for (idx, t) in ungrounded_terminals
            @smart_constraint(model,  [cr, crd, crg, crt, vr[i]],
                                        # [cr, crd, crg, crs, crsw, crt, vr],
                                        sum(cr[a][t] for (a, conns) in bus_arcs if t in conns)
                                        # + sum(crsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                        + sum(crt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                        ==
                                        sum(crg[g][t]         for (g, conns) in bus_gens if t in conns)
                                        # - sum(crs[s][t]         for (s, conns) in bus_storage if t in conns)
                                        - sum(crd[d][t]         for (d, conns) in bus_loads if t in conns)
                                        - sum( Gt[idx,jdx]*vr[i][u] -Bt[idx,jdx]*vi[i][u] for (jdx,u) in ungrounded_terminals) # shunts
                                        )
            @smart_constraint(model, [ci, cid, cig, cit, vi[i]],
                                        # [ci, cid, cig, cis, cisw, cit, vi],
                                        sum(ci[a][t] for (a, conns) in bus_arcs if t in conns)
                                        # + sum(cisw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                        + sum(cit[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                        ==
                                        sum(cig[g][t]         for (g, conns) in bus_gens if t in conns)
                                        # - sum(cis[s][t]         for (s, conns) in bus_storage if t in conns)
                                        - sum(cid[d][t]         for (d, conns) in bus_loads if t in conns)
                                        - sum( Gt[idx,jdx]*vi[i][u] +Bt[idx,jdx]*vr[i][u] for (jdx,u) in ungrounded_terminals) # shunts
                                        )
        end
        
    end


    # "gen connections adaptation of min fuel cost polynomial linquad objective"
    JuMP.@objective(model, Min, sum( gen["cost"][1]*sum(pg[i]) + gen["cost"][2] for (i, gen) in ref[:gen] ) )

end
