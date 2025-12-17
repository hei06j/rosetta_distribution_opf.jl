# BRANCH - Constraints

"""
	function constraint_mc_current_from(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		f_bus::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		g_sh_fr::Matrix{<:Real},
		b_sh_fr::Matrix{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
defines how current distributes over series and shunt impedances of a pi-model branch.

```
cr_fr == csr_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr
ci_fr == csi_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr
```
"""
function constraint_mc_current_from(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, f_bus::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real}; report::Bool=true)
    vr_fr = [PMD.var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr = [PMD.var(pm, nw, :vi, f_bus)[c] for c in f_connections]

    cr_fr =  PMD.var(pm, nw, :cr, f_idx)
    ci_fr =  PMD.var(pm, nw, :ci, f_idx)

    csr_fr =  PMD.var(pm, nw, :csr, f_idx[1])
    csi_fr =  PMD.var(pm, nw, :csi, f_idx[1])

    # JuMP.@constraint(pm.model, cr_fr .== csr_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr)
    # JuMP.@constraint(pm.model, ci_fr .== csi_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr)

    PMD.var(pm, nw, :cr_bus)[f_idx] = cr_bus_fr = PMD._merge_bus_flows(pm, cr_fr, f_connections)
    PMD.var(pm, nw, :ci_bus)[f_idx] = ci_bus_fr = PMD._merge_bus_flows(pm, ci_fr, f_connections)

    if report
        PMD.sol(pm, nw, :branch, f_idx[1])[:pf] =  cr_fr.*vr_fr .+ ci_fr.*vi_fr
        PMD.sol(pm, nw, :branch, f_idx[1])[:qf] = -cr_fr.*vi_fr .+ ci_fr.*vr_fr
    end

end


"""
	function constraint_mc_current_to(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		t_bus,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int},
		g_sh_to::Matrix{<:Real},
		b_sh_to::Matrix{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
defines how current distributes over series and shunt impedances of a pi-model branch.

```
cr_to == csr_to + g_sh_to*vr_to - b_sh_to*vi_to
ci_to == csi_to + g_sh_to*vi_to + b_sh_to*vr_to
```
"""
function constraint_mc_current_to(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, t_bus, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, g_sh_to::Matrix{<:Real}, b_sh_to::Matrix{<:Real}; report::Bool=true)
    vr_to = [PMD.var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to = [PMD.var(pm, nw, :vi, t_bus)[c] for c in t_connections]

    cr_to = PMD.var(pm, nw, :cr, t_idx)
    ci_to = PMD.var(pm, nw, :ci, t_idx)

    csr_to = -PMD.var(pm, nw, :csr, f_idx[1])
    csi_to = -PMD.var(pm, nw, :csi, f_idx[1])

    # JuMP.@constraint(pm.model, cr_to .== csr_to + g_sh_to*vr_to - b_sh_to*vi_to)
    # JuMP.@constraint(pm.model, ci_to .== csi_to + g_sh_to*vi_to + b_sh_to*vr_to)

    PMD.var(pm, nw, :cr_bus)[t_idx] = cr_bus_to = PMD._merge_bus_flows(pm, cr_to, t_connections)
    PMD.var(pm, nw, :ci_bus)[t_idx] = ci_bus_to = PMD._merge_bus_flows(pm, ci_to, t_connections)

    if report
        PMD.sol(pm, nw, :branch, t_idx[1])[:pt] =  cr_to.*vr_to .+ ci_to.*vi_to
        PMD.sol(pm, nw, :branch, t_idx[1])[:qt] = -cr_to.*vi_to .+ ci_to.*vr_to
    end
end


"""
	function constraint_mc_branch_current_limit(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector,
		t_connections::Vector,
		c_rating::Vector{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
imposes a bound on the current magnitude per conductor
at both ends of the branch (total current, i.e. including shunt contributions).

```
cr_fr^2 + ci_fr^2 <= c_rating^2
cr_to^2 + ci_to^2 <= c_rating^2
```
"""
function constraint_mc_branch_current_limit(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector, t_connections::Vector, c_rating::Vector{<:Real}; report::Bool=true)
    cr_fr = PMD.var(pm, nw, :cr, f_idx)
    ci_fr = PMD.var(pm, nw, :ci, f_idx)
    cr_to = PMD.var(pm, nw, :cr, t_idx)
    ci_to = PMD.var(pm, nw, :ci, t_idx)

    cnds_finite_rating = [c for (c,r) in enumerate(c_rating) if r<Inf]
    # JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_fr[c]^2+ci_fr[c]^2 <= c_rating[c]^2)
    # JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_to[c]^2+ci_to[c]^2 <= c_rating[c]^2)
    # JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_fr[c]^2+ci_fr[c]^2 <= sum(c_rating)^2)
    # JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_to[c]^2+ci_to[c]^2 <= sum(c_rating)^2)
    JuMP.@constraint(pm.model, sum(cr_fr[c]^2+ci_fr[c]^2 for c in cnds_finite_rating)  <= sum(c_rating)^2)
    JuMP.@constraint(pm.model, sum(cr_to[c]^2+ci_to[c]^2 for c in cnds_finite_rating)  <= sum(c_rating)^2)
end


"""
	function constraint_mc_thermal_limit_from(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		rate_a::Vector{<:Real}
	)

For IVR models with explicit neutrals,
imposes a bound on the from-side line power magnitude.
"""
function constraint_mc_thermal_limit_from(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})
    vr_fr = [PMD.var(pm, nw, :vr, f_idx[1])[t] for t in f_connections]
    vi_fr = [PMD.var(pm, nw, :vi, f_idx[1])[t] for t in f_connections]
    cr_fr = PMD.var(pm, nw, :cr, f_idx)
    ci_fr = PMD.var(pm, nw, :ci, f_idx)

    for idx in 1:length(rate_a)
        if rate_a[idx]<Inf
            pf_idx = JuMP.@expression(pm.model,  vr_fr[idx]*cr_fr[idx] + vi_fr[idx]*ci_fr[idx])
            qf_idx = JuMP.@expression(pm.model, -vr_fr[idx]*ci_fr[idx] + vi_fr[idx]*cr_fr[idx])

            # JuMP.@constraint(pm.model, pf_idx^2 + qf_idx^2 <= rate_a[idx]^2)
            JuMP.@constraint(pm.model, pf_idx^2 + qf_idx^2 <= sum(rate_a)^2)
        end
    end
end


"""
	function constraint_mc_thermal_limit_to(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		t_idx::Tuple{Int,Int,Int},
		t_connections::Vector{Int},
		rate_a::Vector{<:Real}
	)

For IVR models with explicit neutrals,
imposes a bound on the to-side line power magnitude.
"""
function constraint_mc_thermal_limit_to(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})
    vr_to = [PMD.var(pm, nw, :vr, t_idx[1])[t] for t in t_connections]
    vi_to = [PMD.var(pm, nw, :vi, t_idx[1])[t] for t in t_connections]
    cr_to = PMD.var(pm, nw, :cr, t_idx)
    ci_to = PMD.var(pm, nw, :ci, t_idx)

    for idx in 1:length(rate_a)
        if rate_a[idx]<Inf
            pt_idx = JuMP.@expression(pm.model,  vr_to[idx]*cr_to[idx] + vi_to[idx]*ci_to[idx])
            qt_idx = JuMP.@expression(pm.model, -vr_to[idx]*ci_to[idx] + vi_to[idx]*cr_to[idx])

            # JuMP.@constraint(pm.model, pt_idx^2 + qt_idx^2 <= rate_a[idx]^2)
            JuMP.@constraint(pm.model, pt_idx^2 + qt_idx^2 <= sum(rate_a)^2)
        end
    end
end




# GENERATOR - Constraints - Non-linear

"""
	function constraint_mc_generator_power_wye(
		pm::AbstractNLExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
		pmin::Vector{<:Real},
		pmax::Vector{<:Real},
		qmin::Vector{<:Real},
		qmax::Vector{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
creates non-linear expressions for the generator power `:pd` and `:qd`
of wye-connected generators as a function of voltage and current
"""
function constraint_mc_generator_power_wye(pm::PMD.AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true)
    vr = PMD.var(pm, nw, :vr, bus_id)
    vi = PMD.var(pm, nw, :vi, bus_id)
    crg = PMD.var(pm, nw, :crg, id)
    cig = PMD.var(pm, nw, :cig, id)

    phases = connections[1:end-1]
    n      = connections[end]

    pg = JuMP.NonlinearExpr[]
    qg = JuMP.NonlinearExpr[]

    for (idx, p) in enumerate(phases)
        push!(pg, JuMP.@expression(pm.model,   (vr[p]-vr[n])*crg[idx]  + (vi[p]-vi[n])*cig[idx]))
        push!(qg, JuMP.@expression(pm.model, - (vr[p]-vr[n])*cig[idx]  + (vi[p]-vi[n])*crg[idx]))
    end

    JuMP.@constraint(pm.model, sum(pmin) <= sum(pg))
    JuMP.@constraint(pm.model, sum(pmax) >= sum(pg))
    JuMP.@constraint(pm.model, sum(qmin) <= sum(qg))
    JuMP.@constraint(pm.model, sum(qmax) >= sum(qg))

    # for (idx, p) in enumerate(phases)
    #     if pmin[idx]>-Inf
    #         JuMP.@constraint(pm.model, pmin[idx] .<= (vr[p]-vr[n])*crg[idx]  + (vi[p]-vi[n])*cig[idx])
    #     end
    #     if pmax[idx]< Inf
    #         JuMP.@constraint(pm.model, pmax[idx] .>= (vr[p]-vr[n])*crg[idx]  + (vi[p]-vi[n])*cig[idx])
    #     end
    #     if qmin[idx]>-Inf
    #         JuMP.@constraint(pm.model, qmin[idx] .<= (vi[p]-vi[n])*crg[idx]  - (vr[p]-vr[n])*cig[idx])
    #     end
    #     if qmax[idx]< Inf
    #         JuMP.@constraint(pm.model, qmax[idx] .>= (vi[p]-vi[n])*crg[idx]  - (vr[p]-vr[n])*cig[idx])
    #     end
    # end

    PMD.var(pm, nw, :pg)[id] = pg
    PMD.var(pm, nw, :qg)[id] = qg

    if report
        PMD.sol(pm, nw, :gen, id)[:pg] = pg
        PMD.sol(pm, nw, :gen, id)[:qg] = qg
    end
end


"""
	function constraint_mc_generator_power_delta(
		pm::AbstractNLExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
		pmin::Vector{<:Real},
		pmax::Vector{<:Real},
		qmin::Vector{<:Real},
		qmax::Vector{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
creates non-linear expressions for the generator power `:pd` and `:qd`
of delta-connected generators as a function of voltage and current
"""
### TODO update to have pg and qg limits, similar to wye case above
function constraint_mc_generator_power_delta(pm::PMD.AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true)
    vr = PMD.var(pm, nw, :vr, bus_id)
    vi = PMD.var(pm, nw, :vi, bus_id)
    crg = PMD.var(pm, nw, :crg, id)
    cig = PMD.var(pm, nw, :cig, id)

    nph = length(pmin)

    vrg = Dict()
    vig = Dict()
    for (idx,c,d) in zip(1:nph, connections, [connections[2:end]..., connections[1]])
        vrg[idx] = JuMP.@expression(pm.model, vr[c]-vr[d])
        vig[idx] = JuMP.@expression(pm.model, vi[c]-vi[d])
    end

    pg = JuMP.NonlinearExpr[]
    qg = JuMP.NonlinearExpr[]
    
    for idx in 1:nph
        push!(pg, JuMP.@expression(pm.model,  vrg[idx]*crg[idx]+vig[idx]*cig[idx]))
        push!(qg, JuMP.@expression(pm.model, -vrg[idx]*cig[idx]+vig[idx]*crg[idx]))
    end

    JuMP.@constraint(pm.model, [i in 1:nph], pmin[i] <= pg[i])
    JuMP.@constraint(pm.model, [i in 1:nph], pmax[i] >= pg[i])
    JuMP.@constraint(pm.model, [i in 1:nph], qmin[i] <= qg[i])
    JuMP.@constraint(pm.model, [i in 1:nph], qmax[i] >= qg[i])

    PMD.var(pm, nw, :pg)[id] = JuMP.Containers.DenseAxisArray(pg, connections)
    PMD.var(pm, nw, :qg)[id] = JuMP.Containers.DenseAxisArray(qg, connections)

    if report
        PMD.sol(pm, nw, :gen, id)[:pg] = pg
        PMD.sol(pm, nw, :gen, id)[:qg] = qg
    end
end


"""
	function constraint_mc_generator_current_limit(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		connections::Vector{Int};
		report::Bool=true,
		bounded::Bool=true
	)

For IVR models with explicit neutrals,
creates expressions for the terminal current flows `:crg_bus` and `:cig_bus` of wye-connected generators
"""
function constraint_mc_generator_current_limit(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, id::Int, connections::Vector{Int}, c_rating::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
    # crg = PMD.var(pm, nw, :crg, id)
    # cig = PMD.var(pm, nw, :cig, id)

    crg_bus = PMD.var(pm, nw, :crg_bus)[id]
    cig_bus = PMD.var(pm, nw, :cig_bus)[id]
    
    # PMD.var(pm, nw, :crg_bus)[id] = crg_bus = PMD._merge_bus_flows(pm, [crg..., -sum(crg)], connections)
    # PMD.var(pm, nw, :cig_bus)[id] = cig_bus = PMD._merge_bus_flows(pm, [cig..., -sum(cig)], connections)
    # JuMP.@constraint(pm.model, sum(PMD.var(pm, nw, :crg_bus)[id]) == 0)
    # JuMP.@constraint(pm.model, sum(PMD.var(pm, nw, :cig_bus)[id]) == 0)

    @assert length(c_rating) == length(crg_bus)
    cnds_finite_nonzero_rating = [c for (c,r) in enumerate(c_rating) if (r<Inf && r!==0)]
    cnds_zero_rating = [c for (c,r) in enumerate(c_rating) if r==0]

    # JuMP.@constraint(pm.model, [c in cnds_finite_rating], crg[c]^2+cig[c]^2 <= c_rating[c]^2)
    JuMP.@constraint(pm.model, [c in cnds_finite_nonzero_rating], crg_bus[c]^2+cig_bus[c]^2 <= c_rating[c]^2)
    JuMP.@constraint(pm.model, [c in cnds_zero_rating], crg_bus[c] == 0)
    JuMP.@constraint(pm.model, [c in cnds_zero_rating], cig_bus[c] == 0)
end


"""
	function constraint_mc_generator_current_wye(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		connections::Vector{Int};
		report::Bool=true,
		bounded::Bool=true
	)

For IVR models with explicit neutrals,
creates expressions for the terminal current flows `:crg_bus` and `:cig_bus` of wye-connected generators
"""
function constraint_mc_generator_current_wye(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crg = PMD.var(pm, nw, :crg, id)
    cig = PMD.var(pm, nw, :cig, id)
    PMD.var(pm, nw, :crg_bus)[id] = PMD._merge_bus_flows(pm, [crg..., -sum(crg)], connections)
    PMD.var(pm, nw, :cig_bus)[id] = PMD._merge_bus_flows(pm, [cig..., -sum(cig)], connections)

    JuMP.@constraint(pm.model, sum(PMD.var(pm, nw, :crg_bus)[id]) == 0)
    JuMP.@constraint(pm.model, sum(PMD.var(pm, nw, :cig_bus)[id]) == 0)
end


"""
	function constraint_mc_generator_current_delta(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		connections::Vector{Int};
		report::Bool=true,
		bounded::Bool=true
	)

For IVR models with explicit neutrals,
creates expressions for the terminal current flows `:crg_bus` and `:cig_bus` of delta-connected generators
"""
function constraint_mc_generator_current_delta(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crg = PMD.var(pm, nw, :crg, id)
    cig = PMD.var(pm, nw, :cig, id)
    Md = _PMD.get_delta_transformation_matrix(length(connections))
    PMD.var(pm, nw, :crg_bus)[id] = PMD._merge_bus_flows(pm, Md'*crg, connections)
    PMD.var(pm, nw, :cig_bus)[id] = PMD._merge_bus_flows(pm, Md'*cig, connections)
end


# function constraint_mc_load_current_wye_magnitude(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
function constraint_mc_load_current_wye_magnitude(pm::PMD.AbstractExplicitNeutralIVRModel, id::Int; nw::Int=PMD.nw_id_default, report::Bool=true)
    load = PMD.ref(pm, nw, :load, id)
    bus = PMD.ref(pm, nw,:bus, load["load_bus"])
    load_current = load["cm"]
    connections = load["connections"]
    phases = connections[1:end-1]

    ccmd = JuMP.NonlinearExpr[]
    crd = PMD.var(pm, nw, :crd)[id]
    cid = PMD.var(pm, nw, :cid)[id]

    for (idx,c) in enumerate(phases)
        push!(ccmd, JuMP.@expression(pm.model,  crd[idx]^2 + cid[idx]^2))
        JuMP.@constraint(pm.model, crd[idx]^2 + cid[idx]^2 == load_current[idx]^2)
    end
    PMD.sol(pm, nw, :load, id)[:ccmd] = JuMP.Containers.DenseAxisArray(ccmd, phases)
end

"""
	function constraint_mc_load_current_wye(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
		a::Vector{<:Real},
		alpha::Vector{<:Real},
		b::Vector{<:Real},
		beta::Vector{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
create non-linear expressions for the terminal current flows `:crd_bus` and `:cid_bus`
of wye-connected loads
"""
function constraint_mc_load_current_wye(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, load_current, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = PMD.var(pm, nw, :vr, bus_id)
    vi = PMD.var(pm, nw, :vi, bus_id)

    crd = JuMP.NonlinearExpr[]
    cid = JuMP.NonlinearExpr[]
    ccmd = JuMP.NonlinearExpr[]

    phases = connections[1:end-1]
    n      = connections[end]

    for (idx, p) in enumerate(phases)
        push!(crd, JuMP.@expression(pm.model,
             a[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1)
            +b[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1)
        ))
        push!(cid, JuMP.@expression(pm.model,
             a[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1)
            -b[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1)
        ))
        push!(ccmd, JuMP.@expression(pm.model,
            (a[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1)
            +b[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1))^2 
            +
            (a[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1)
            -b[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1))^2
        ))
        JuMP.@constraint(pm.model, load_current[idx]^2 ==
            (a[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1)
            +b[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1))^2 
            +
            (a[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1)
            -b[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1))^2
            )
    end  

    PMD.var(pm, nw, :crd)[id] = crd
    PMD.var(pm, nw, :cid)[id] = cid

    crd_bus_n = JuMP.@expression(pm.model, -sum(crd[i] for i in 1:length(phases)))
    cid_bus_n = JuMP.@expression(pm.model, -sum(cid[i] for i in 1:length(phases)))

    PMD.var(pm, nw, :crd_bus)[id] = crd_bus = PMD._merge_bus_flows(pm, [crd..., crd_bus_n], connections)
    PMD.var(pm, nw, :cid_bus)[id] = cid_bus = PMD._merge_bus_flows(pm, [cid..., cid_bus_n], connections)

    if report
        pd_bus = JuMP.NonlinearExpr[]
        qd_bus = JuMP.NonlinearExpr[]
        for (idx,c) in enumerate(connections)
            push!(pd_bus, JuMP.@expression(pm.model,  vr[c]*crd_bus[c]+vi[c]*cid_bus[c]))
            push!(qd_bus, JuMP.@expression(pm.model, -vr[c]*cid_bus[c]+vi[c]*crd_bus[c]))
        end

        PMD.sol(pm, nw, :load, id)[:pd_bus] = JuMP.Containers.DenseAxisArray(pd_bus, connections)
        PMD.sol(pm, nw, :load, id)[:qd_bus] = JuMP.Containers.DenseAxisArray(qd_bus, connections)

        PMD.sol(pm, nw, :load, id)[:crd] = JuMP.Containers.DenseAxisArray(crd, connections)
        PMD.sol(pm, nw, :load, id)[:cid] = JuMP.Containers.DenseAxisArray(cid, connections)
        PMD.sol(pm, nw, :load, id)[:ccmd] = JuMP.Containers.DenseAxisArray(ccmd, connections)

        PMD.sol(pm, nw, :load, id)[:crd_bus] = crd_bus
        PMD.sol(pm, nw, :load, id)[:cid_bus] = cid_bus

        pd = JuMP.NonlinearExpr[]
        qd = JuMP.NonlinearExpr[]
        for (idx, p) in enumerate(phases)
            push!(pd, JuMP.@expression(pm.model, a[idx]*(vr[p]^2+vi[p]^2)^(alpha[idx]/2) ))
            push!(qd, JuMP.@expression(pm.model, b[idx]*(vr[p]^2+vi[p]^2)^(beta[idx]/2)  ))
        end
        PMD.sol(pm, nw, :load, id)[:pd] = JuMP.Containers.DenseAxisArray(pd, connections)
        PMD.sol(pm, nw, :load, id)[:qd] = JuMP.Containers.DenseAxisArray(qd, connections)
    end
end


"""
	function constraint_mc_load_current_delta(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
		a::Vector{<:Real},
		alpha::Vector{<:Real},
		b::Vector{<:Real},
		beta::Vector{<:Real};
		report::Bool=true
	)

For IVR models with explicit neutrals,
create non-linear expressions for the terminal current flows `:crd_bus` and `:cid_bus`
of delta-connected loads
"""
function constraint_mc_load_current_delta(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = PMD.var(pm, nw, :vr, bus_id)
    vi = PMD.var(pm, nw, :vi, bus_id)


    ph = connections
    ph_next = [connections[2:end]..., connections[1]]
    P = length(ph)
    idxs = 1:P
    idxs_prev = [idxs[end], idxs[1:end-1]...]

    vrd = [vr[c]-vr[d] for (c,d) in zip(ph,ph_next)]
    vid = [vi[c]-vi[d] for (c,d) in zip(ph,ph_next)]

    crd = JuMP.@expression(pm.model, [i in 1:P],
        a[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1)
       +b[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1)
    )
    cid = JuMP.@expression(pm.model, [i in 1:P],
        a[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1)
       -b[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1)
    )

    crd_bus = JuMP.@expression(pm.model, [i in 1:P], crd[i]-crd[idxs_prev[i]])
    cid_bus = JuMP.@expression(pm.model, [i in 1:P], cid[i]-cid[idxs_prev[i]])

    PMD.var(pm, nw, :crd_bus)[id] = PMD._merge_bus_flows(pm, crd_bus, connections)
    PMD.var(pm, nw, :cid_bus)[id] = PMD._merge_bus_flows(pm, cid_bus, connections)

    if report
        pd_bus = JuMP.@expression(pm.model, [i in 1:P],  vr[i]*crd_bus[i]+vi[i]*cid_bus[i])
        qd_bus = JuMP.@expression(pm.model, [i in 1:P], -vr[i]*cid_bus[i]+vi[i]*crd_bus[i])

        PMD.sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        PMD.sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        pd = JuMP.@expression(pm.model, [i in 1:P], a[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2) )
        qd = JuMP.@expression(pm.model, [i in 1:P], b[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2)  )
        PMD.sol(pm, nw, :load, id)[:pd] = pd
        PMD.sol(pm, nw, :load, id)[:qd] = qd
    end
end


"""
	function constraint_mc_inverter_dc_link_ripple_power(
		pm::AbstractNLExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
        pdcmin::Vector{<:Real},
		pdcmax::Vector{<:Real},
		report::Bool=true
	)

For IVR models with explicit neutrals,
creates non-linear expressions for the inverter dc link power `:pdc_link_sqr`
of wye-connected generators as a function of voltage and current
"""
function constraint_mc_inverter_dc_link_ripple_power(pm::PMD.AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pdcmin::Real, pdcmax::Real; report::Bool=true)
    # bus_id = 1
    vr = PMD.var(pm, nw, :vr, bus_id)
    vi = PMD.var(pm, nw, :vi, bus_id)
    crg = PMD.var(pm, nw, :crg, id)
    cig = PMD.var(pm, nw, :cig, id)
    crg_bus = PMD.var(pm, nw, :crg_bus)[id]
    cig_bus = PMD.var(pm, nw, :cig_bus)[id]
    
    phases = connections[1:end-1]
    n      = connections[end]

    
    if pdcmax > 0
        if pdcmax < Inf
            JuMP.@constraint(pm.model, pdcmax^2 >= sum( vr[p]*crg_bus[idx] - vi[p]*cig_bus[idx] for (idx, p) in enumerate(connections) )^2
                                                    + 
                                                   sum( vr[p]*cig_bus[idx] + vi[p]*crg_bus[idx] for (idx, p) in enumerate(connections) )^2
                            )
        end
    elseif pdcmax == 0
        JuMP.@constraint(pm.model, sum( vr[p]*crg_bus[idx] - vi[p]*cig_bus[idx] for (idx, p) in enumerate(connections) ) == 0)
        JuMP.@constraint(pm.model, sum( vr[p]*cig_bus[idx] + vi[p]*crg_bus[idx] for (idx, p) in enumerate(connections) ) == 0)
    end
    
    pdc_link_sqr = JuMP.@expression(pm.model,  
        sum( vr[p]*crg_bus[idx] - vi[p]*cig_bus[idx] for (idx, p) in enumerate(connections) )^2
        + 
        sum( vr[p]*cig_bus[idx] + vi[p]*crg_bus[idx] for (idx, p) in enumerate(connections) )^2
    )
    
    PMD.var(pm, nw, :pdc_link_sqr)[id] = pdc_link_sqr
    
    if report
        PMD.sol(pm, nw, :gen, id)[:pdc_link_sqr] = pdc_link_sqr
    end

end




"""
	function constraint_mc_inverter_branch_dc_link_ripple_power(
		pm::AbstractNLExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
        pdcmin::Vector{<:Real},
		pdcmax::Vector{<:Real},
		report::Bool=true
	)

For IVR models with explicit neutrals,
creates non-linear expressions for the inverter dc link power `:pdc_link_sqr`
of wye-connected generators as a function of voltage and current
"""
function constraint_mc_inverter_branch_dc_link_ripple_power(pm::PMD.AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, pdcmin::Real, pdcmax::Real; report::Bool=true)
    ### use from side voltage and currents
    vr = [PMD.var(pm, nw, :vr, f_idx[1])[t] for t in f_connections]
    vi = [PMD.var(pm, nw, :vi, f_idx[1])[t] for t in f_connections]
    cr = PMD.var(pm, nw, :cr, f_idx)
    ci = PMD.var(pm, nw, :ci, f_idx)
    
    phases = f_connections[1:end-1]
    n      = f_connections[end]
    

    pdc_link_sqr = JuMP.@expression(pm.model,  
        sum( vr[p]*cr[idx] - vi[p]*ci[idx] for (idx, p) in enumerate(f_connections) )^2
        + 
        sum( vr[p]*ci[idx] + vi[p]*cr[idx] for (idx, p) in enumerate(f_connections) )^2
    )

    # if pdcmin !== pdcmax
        if pdcmin > -Inf
            JuMP.@constraint(pm.model, pdcmin^2 <= sum( vr[p]*cr[idx] - vi[p]*ci[idx] for (idx, p) in enumerate(f_connections) )^2
                                                    + 
                                                    sum( vr[p]*ci[idx] + vi[p]*cr[idx] for (idx, p) in enumerate(f_connections) )^2
                            )
        end
        if pdcmax < Inf
            JuMP.@constraint(pm.model, pdcmax^2 >= sum( vr[p]*cr[idx] - vi[p]*ci[idx] for (idx, p) in enumerate(f_connections) )^2
                                                    + 
                                                    sum( vr[p]*ci[idx] + vi[p]*cr[idx] for (idx, p) in enumerate(f_connections) )^2
                            )
        end

    PMD.var(pm, nw, :pdc_link_sqr)[id] = pdc_link_sqr

    if report
        PMD.sol(pm, nw, :branch, id)[:pdc_link_sqr] = pdc_link_sqr
    end
end


###### voltage sequence components
"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_negative_sequence(pm::PMD.AbstractUnbalancedIVRModel, nw::Int, bus_id::Int, vmnegmax::Real)
    # if !haskey(PMD.var(pm, nw_id_default), :vmpossqr)
    #     PMD.var(pm, nw_id_default)[:vmpossqr] = Dict{Int, Any}()
    #     PMD.var(pm, nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
    # end
    (vr_a, vr_b, vr_c) = [PMD.var(pm, nw, :vr, bus_id)[i] for i in 1:3]
    (vi_a, vi_b, vi_c) = [PMD.var(pm, nw, :vi, bus_id)[i] for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)
    
    # real and imaginary components of U-
    vreneg = JuMP.@expression(pm.model,
        (vr_a + a2re*vr_b - a2im*vi_b + are*vr_c - aim*vi_c)/3
    )
    vimneg = JuMP.@expression(pm.model,
        (vi_a + a2re*vi_b + a2im*vr_b + are*vi_c + aim*vr_c)/3
    )
    # square of magnitude of U-, |U-|^2
    vmnegsqr = JuMP.@expression(pm.model, vreneg^2 + vimneg^2)


    PMD.var(pm, nw, :vmnegsqr)[bus_id] = vmnegsqr
    PMD.sol(pm, nw, :bus, bus_id)[:vmnegsqr] = vmnegsqr

    vmneg = PMD.var(pm, nw, :vmneg)[bus_id] = JuMP.@variable(pm.model, base_name="$(nw)_vmneg_$bus_id", start = 0, lower_bound=0)
    PMD.sol(pm, nw, :bus, bus_id)[:vmneg] = vmneg
    
    JuMP.@constraint(pm.model, vmneg * vmneg == vmnegsqr)

    # # finally, apply constraint
    # JuMP.@constraint(pm.model, vmnegsqr <= vmnegmax^2)
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_positive_sequence(pm::PMD.AbstractUnbalancedIVRModel, nw::Int, bus_id::Int, vmposmax::Real)
    if !haskey(PMD.var(pm, nw_id_default), :vmpossqr)
        PMD.var(pm, nw_id_default)[:vmpossqr] = Dict{Int, Any}()
        PMD.var(pm, nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vr_a, vr_b, vr_c) = [PMD.var(pm, nw, :vr, bus_id)[i] for i in 1:3]
    (vi_a, vi_b, vi_c) = [PMD.var(pm, nw, :vi, bus_id)[i] for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)
    # real and imaginary components of U+
    vrepos = JuMP.@expression(pm.model,
        (vr_a + are*vr_b - aim*vi_b + a2re*vr_c - a2im*vi_c)/3
    )
    vimpos = JuMP.@expression(pm.model,
        (vi_a + are*vi_b + aim*vr_b + a2re*vi_c + a2im*vr_c)/3
    )
    # square of magnitude of U+, |U+|^2
    vmpossqr = JuMP.@expression(pm.model, vrepos^2+vimpos^2)
    # finally, apply constraint
    JuMP.@constraint(pm.model, vmpossqr <= vmposmax^2)
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_zero_sequence(pm::PMD.AbstractUnbalancedIVRModel, nw::Int, bus_id::Int, vmzeromax::Real)
    if !haskey(PMD.var(pm, nw_id_default), :vmpossqr)
        PMD.var(pm, nw_id_default)[:vmpossqr] = Dict{Int, Any}()
        PMD.var(pm, nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vr_a, vr_b, vr_c) = [PMD.var(pm, nw, :vr, bus_id)[i] for i in 1:3]
    (vi_a, vi_b, vi_c) = [PMD.var(pm, nw, :vi, bus_id)[i] for i in 1:3]
    # real and imaginary components of U+
    vrezero = JuMP.@expression(pm.model,
        (vr_a + vr_b + vr_c)/3
    )
    vimzero = JuMP.@expression(pm.model,
        (vi_a + vi_b + vi_c)/3
    )
    # square of magnitude of U+, |U+|^2
    vmzerosqr = JuMP.@expression(pm.model, vrezero^2+vimzero^2)
    # finally, apply constraint
    JuMP.@constraint(pm.model, vmzerosqr <= vmzeromax^2)
end



"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_vuf(pm::PMD.AbstractUnbalancedIVRModel, nw::Int, bus_id::Int, vufmax::Real)
    if !haskey(PMD.var(pm, nw_id_default), :vmpossqr)
        PMD.var(pm, nw_id_default)[:vmpossqr] = Dict{Int, Any}()
        PMD.var(pm, nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vr_a, vr_b, vr_c) = [PMD.var(pm, nw, :vr, bus_id)[i] for i in 1:3]
    (vi_a, vi_b, vi_c) = [PMD.var(pm, nw, :vi, bus_id)[i] for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)
    # real and imaginary components of U+
    vrepos = JuMP.@expression(pm.model,
        (vr_a + are*vr_b - aim*vi_b + a2re*vr_c - a2im*vi_c)/3
    )
    vimpos = JuMP.@expression(pm.model,
        (vi_a + are*vi_b + aim*vr_b + a2re*vi_c + a2im*vr_c)/3
    )
    # square of magnitude of U+, |U+|^2
    vmpossqr = JuMP.@expression(pm.model, vrepos^2+vimpos^2)
    # real and imaginary components of U-
    vreneg = JuMP.@expression(pm.model,
        (vr_a + a2re*vr_b - a2im*vi_b + are*vr_c - aim*vi_c)/3
    )
    vimneg = JuMP.@expression(pm.model,
        (vi_a + a2re*vi_b + a2im*vr_b + are*vi_c + aim*vr_c)/3
    )
    # square of magnitude of U-, |U-|^2
    vmnegsqr = JuMP.@expression(pm.model, vreneg^2+vimneg^2)
    # finally, apply constraint
    JuMP.@constraint(pm.model, vmnegsqr <= vufmax^2*vmpossqr)
    # DEBUGGING: save references for post check
    #PMD.var(pm, nw_id_default, :vmpossqr)[bus_id] = vmpossqr
    #PMD.var(pm, nw_id_default, :vmnegsqr)[bus_id] = vmnegsqr
end



"""
    constraint_mc_bus_voltage_balance(pm::AbstractUnbalancedACRModel, bus_id::Int; nw=nw_id_default)::Nothing

Template function for bus voltage balance constraints.
"""
function constraint_mc_bus_voltage_balance(pm::PMD.AbstractUnbalancedACRModel, bus_id::Int; nw=PMD.nw_id_default)::Nothing
    # @assert(length(PMD.ref(pm, nw, :conductor_ids))==3)

    bus = PMD.ref(pm, nw, :bus, bus_id)
    constraint_mc_bus_voltage_magnitude_negative_sequence(pm, nw, bus_id, 0)

    # if haskey(bus, "vm_vuf_max")
    #     constraint_mc_bus_voltage_magnitude_vuf(pm, nw, bus_id, bus["vm_vuf_max"])
    # end

    # if haskey(bus, "vm_seq_neg_max")
    #     constraint_mc_bus_voltage_magnitude_negative_sequence(pm, nw, bus_id, bus["vm_seq_neg_max"])
    # end

    # if haskey(bus, "vm_seq_pos_max")
    #     constraint_mc_bus_voltage_magnitude_positive_sequence(pm, nw, bus_id, bus["vm_seq_pos_max"])
    # end

    # if haskey(bus, "vm_seq_zero_max")
    #     constraint_mc_bus_voltage_magnitude_zero_sequence(pm, nw, bus_id, bus["vm_seq_zero_max"])
    # end

    # if haskey(bus, "vm_ll_min")|| haskey(bus, "vm_ll_max")
    #     vm_ll_min = haskey(bus, "vm_ll_min") ? bus["vm_ll_min"] : fill(0, 3)
    #     vm_ll_max = haskey(bus, "vm_ll_max") ? bus["vm_ll_max"] : fill(Inf, 3)
    #     PMD.constraint_mc_bus_voltage_magnitude_ll(pm, nw, bus_id, vm_ll_min, vm_ll_max)
    # end
    nothing
end

