# BRANCH - Constraints

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
function constraint_mc_branch_current_limit(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector, t_connections::Vector, c_rating::Vector{<:Real}; report::Bool=true)
    cr_fr = _PMD.var(pm, nw, :cr, f_idx)
    ci_fr = _PMD.var(pm, nw, :ci, f_idx)
    cr_to = _PMD.var(pm, nw, :cr, t_idx)
    ci_to = _PMD.var(pm, nw, :ci, t_idx)

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
function constraint_mc_thermal_limit_from(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})
    vr_fr = [_PMD.var(pm, nw, :vr, f_idx[1])[t] for t in f_connections]
    vi_fr = [_PMD.var(pm, nw, :vi, f_idx[1])[t] for t in f_connections]
    cr_fr = _PMD.var(pm, nw, :cr, f_idx)
    ci_fr = _PMD.var(pm, nw, :ci, f_idx)

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
function constraint_mc_thermal_limit_to(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})
    vr_to = [_PMD.var(pm, nw, :vr, t_idx[1])[t] for t in t_connections]
    vi_to = [_PMD.var(pm, nw, :vi, t_idx[1])[t] for t in t_connections]
    cr_to = _PMD.var(pm, nw, :cr, t_idx)
    ci_to = _PMD.var(pm, nw, :ci, t_idx)

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
function constraint_mc_generator_power_wye(pm::_PMD.AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true)
    vr = _PMD.var(pm, nw, :vr, bus_id)
    vi = _PMD.var(pm, nw, :vi, bus_id)
    crg = _PMD.var(pm, nw, :crg, id)
    cig = _PMD.var(pm, nw, :cig, id)

    phases = connections[1:end-1]
    n      = connections[end]

    pg = JuMP.NonlinearExpr[]
    qg = JuMP.NonlinearExpr[]

    for (idx, p) in enumerate(phases)
        push!(pg, JuMP.@expression(pm.model,  (vr[p]-vr[n])*crg[idx]+(vi[p]-vi[n])*cig[idx]))
        push!(qg, JuMP.@expression(pm.model, -(vr[p]-vr[n])*cig[idx]+(vi[p]-vi[n])*crg[idx]))
    end

    for (idx, p) in enumerate(phases)
        if pmin[idx]>-Inf
            JuMP.@constraint(pm.model, pmin[idx] .<= (vr[p]-vr[n])*crg[idx]  + (vi[p]-vi[n])*cig[idx])
        end
        if pmax[idx]< Inf
            JuMP.@constraint(pm.model, pmax[idx] .>= (vr[p]-vr[n])*crg[idx]  + (vi[p]-vi[n])*cig[idx])
        end
        if qmin[idx]>-Inf
            JuMP.@constraint(pm.model, qmin[idx] .<= (vi[p]-vi[n])*crg[idx]  - (vr[p]-vr[n])*cig[idx])
        end
        if qmax[idx]< Inf
            JuMP.@constraint(pm.model, qmax[idx] .>= (vi[p]-vi[n])*crg[idx]  - (vr[p]-vr[n])*cig[idx])
        end
    end

    _PMD.var(pm, nw, :pg)[id] = pg
    _PMD.var(pm, nw, :qg)[id] = qg

    if report
        _PMD.sol(pm, nw, :gen, id)[:pg] = pg
        _PMD.sol(pm, nw, :gen, id)[:qg] = qg
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
function constraint_mc_generator_power_delta(pm::_PMD.AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true)
    vr = _PMD.var(pm, nw, :vr, bus_id)
    vi = _PMD.var(pm, nw, :vi, bus_id)
    crg = _PMD.var(pm, nw, :crg, id)
    cig = _PMD.var(pm, nw, :cig, id)

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

    _PMD.var(pm, nw, :pg)[id] = JuMP.Containers.DenseAxisArray(pg, connections)
    _PMD.var(pm, nw, :qg)[id] = JuMP.Containers.DenseAxisArray(qg, connections)

    if report
        _PMD.sol(pm, nw, :gen, id)[:pg] = pg
        _PMD.sol(pm, nw, :gen, id)[:qg] = qg
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
creates non-linear expressions for the inverter dc link power `:pdc_link`
of wye-connected generators as a function of voltage and current
"""
function constraint_mc_inverter_dc_link_ripple_power(pm::_PMD.AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pdcmin::Real, pdcmax::Real; report::Bool=true)
    vr = _PMD.var(pm, nw, :vr, bus_id)
    vi = _PMD.var(pm, nw, :vi, bus_id)
    crg = _PMD.var(pm, nw, :crg, id)
    cig = _PMD.var(pm, nw, :cig, id)

    phases = connections[1:end-1]
    n      = connections[end]
    
    pdc_link = JuMP.@expression(pm.model,  
        sqrt(sum( ((vr[p]-vr[n])*crg[idx]-(vi[p]-vi[n])*cig[idx])^2 + 
                ((vr[p]-vr[n])*cig[idx]+(vi[p]-vi[n])*crg[idx])^2 
            for (idx, p) in enumerate(phases))
            )
        )
    
    if pdcmin > -Inf
        JuMP.@constraint(pm.model, pdcmin^2 <= sum( ((vr[p]-vr[n])*crg[idx]-(vi[p]-vi[n])*cig[idx])^2 + 
                                                        ((vr[p]-vr[n])*cig[idx]+(vi[p]-vi[n])*crg[idx])^2 
                                                    for (idx, p) in enumerate(phases))
                        )
    end
    if pdcmax < Inf
        JuMP.@constraint(pm.model, pdcmax^2 >= sum( ((vr[p]-vr[n])*crg[idx]-(vi[p]-vi[n])*cig[idx])^2 + 
                                                        ((vr[p]-vr[n])*cig[idx]+(vi[p]-vi[n])*crg[idx])^2 
                                                    for (idx, p) in enumerate(phases))
                        )
    end
    
    _PMD.var(pm, nw, :pdc_link)[id] = pdc_link

    if report
        _PMD.sol(pm, nw, :gen, id)[:pdc_link] = pdc_link
    end
end