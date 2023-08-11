"""
	function constraint_mc_voltage_absolute(
		pm::_PMD.RectangularVoltageExplicitNeutralModels,
		nw::Int,
		i::Int,
		terminals::Vector{Int},
		grounded::Vector{Bool},
		vmin::Vector{<:Real},
		vmax::Vector{<:Real};
		report::Bool=true
	)

Imposes absolute voltage magnitude bounds for models with explicit neutrals
"""
function constraint_mc_voltage_absolute(pm::_PMD.RectangularVoltageExplicitNeutralModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, vmin::Vector{<:Real}, vmax::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    ungrounded_terminals = terminals[(!).(grounded)]
    for (idx,t) in enumerate(terminals)
        if !grounded[idx]
            if vmax[idx] < Inf
                JuMP.@constraint(pm.model, vr[t]^2+vi[t]^2 <= vmax[idx]^2)
            end
            if vmin[idx] > 0.0
                JuMP.@constraint(pm.model, vr[t]^2+vi[t]^2 >= vmin[idx]^2)
            end
        end
    end
end


"""
	function constraint_mc_voltage_pairwise(
		pm::_PMD.RectangularVoltageExplicitNeutralModels,
		nw::Int,
		i::Int,
		terminals::Vector{Int},
		grounded::Vector{Bool},
		vm_pair_lb::Vector,
		vm_pair_ub::Vector;
		report::Bool=true
	)

Imposes pairwise voltage magnitude bounds, i.e. magnitude bounds on the voltage between to terminals, for models with explicit neutrals
"""
function constraint_mc_voltage_pairwise(pm::_PMD.RectangularVoltageExplicitNeutralModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, vm_pair_lb::Vector{<:Tuple{Any,Any,Real}}, vm_pair_ub::Vector{<:Tuple{Any,Any,Real}}; report::Bool=true)
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    for (a,b,lb) in vm_pair_lb
        if lb > 0.0
            JuMP.@constraint(pm.model, (vr[a]-vr[b])^2 + (vi[a]-vi[b])^2 >= lb^2)
        end
    end

    for (a,b,ub) in vm_pair_ub
        if ub < Inf
            JuMP.@constraint(pm.model, (vr[a]-vr[b])^2 + (vi[a]-vi[b])^2 <= ub^2)
        end
    end
end


"""
	function constraint_mc_theta_ref(
		pm::_PMD.RectangularVoltageExplicitNeutralModels,
		nw::Int,
		i::Int,
		va::Vector{<:Real},
		terminals::Vector{Int},
		grounded::Vector{Bool}
	)

Creates phase angle constraints at bus `i`
"""
function constraint_mc_theta_ref(pm::_PMD.RectangularVoltageExplicitNeutralModels, nw::Int, i::Int, va::Vector{<:Real}, terminals::Vector{Int}, grounded::Vector{Bool})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    # deal with cases first where tan(theta)==Inf or tan(theta)==0
    for idx in findall[(!).(grounded)]
        t = terminals[idx]
        if va[t] == pi/2
            JuMP.@constraint(pm.model, vr[t] == 0)
            JuMP.@constraint(pm.model, vi[t] >= 0)
        elseif va[t] == -pi/2
            JuMP.@constraint(pm.model, vr[t] == 0)
            JuMP.@constraint(pm.model, vi[t] <= 0)
        elseif va[t] == 0
            JuMP.@constraint(pm.model, vr[t] >= 0)
            JuMP.@constraint(pm.model, vi[t] == 0)
        elseif va[t] == pi
            JuMP.@constraint(pm.model, vr[t] >= 0)
            JuMP.@constraint(pm.model, vi[t] == 0)
        else
            JuMP.@constraint(pm.model, vi[t] == tan(va[t])*vr[t])
            # va_ref also implies a sign for vr, vi
            if 0<=va[t] && va[t] <= pi
                JuMP.@constraint(pm.model, vi[t] >= 0)
            else
                JuMP.@constraint(pm.model, vi[t] <= 0)
            end
        end
    end
end



"""
	function constraint_mc_voltage_fixed(
		pm::_PMD.RectangularVoltageExplicitNeutralModels,
		nw::Int,
		i::Int,
		vm::Vector{<:Real},
		va::Vector{<:Real},
		terminals::Vector{Int},
		grounded::Vector{Bool}
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
"""
function constraint_mc_voltage_fixed(pm::_PMD.RectangularVoltageExplicitNeutralModels, nw::Int, i::Int, vm::Vector{<:Real}, va::Vector{<:Real}, terminals::Vector{Int}, grounded::Vector{Bool})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    idxs = findall((!).(grounded))

    JuMP.@constraint(pm.model, [i in idxs], vr[terminals[i]]==vm[i]*cos(va[i]))
    JuMP.@constraint(pm.model, [i in idxs], vi[terminals[i]]==vm[i]*sin(va[i]))
end


"""
	function constraint_mc_voltage_magnitude_fixed(
		pm::_PMD.RectangularVoltageExplicitNeutralModels,
		nw::Int,
		i::Int,
		vm::Vector{<:Real},
		va::Vector{<:Real},
		terminals::Vector{Int},
		grounded::Vector{Bool}
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
"""
function constraint_mc_voltage_magnitude_fixed(pm::_PMD.RectangularVoltageExplicitNeutralModels, nw::Int, i::Int, vm::Vector{<:Real}, terminals::Vector{Int}, grounded::Vector{Bool})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    @assert iszero(vm[grounded]) "Infeasible model; the voltage magnitude of a grounded terminal is fixed to a non-zero value."
    idxs = findall((!).(grounded))

    JuMP.@constraint(pm.model, [i in idxs], vr[terminals[i]].^2 + vi[terminals[i]].^2 == vm[i].^2)
end



# GENERATOR - Constraints - Non-linear

"""
	function constraint_mc_generator_power_wye(
		pm::_PMD.AbstractNLExplicitNeutralIVRModel,
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
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    phases = connections[1:end-1]
    n      = connections[end]

    pg = Vector{JuMP.NonlinearExpression}([])
    qg = Vector{JuMP.NonlinearExpression}([])

    for (idx, p) in enumerate(phases)
        push!(pg, JuMP.@NLexpression(pm.model,  (vr[p]-vr[n])*crg[idx]+(vi[p]-vi[n])*cig[idx]))
        push!(qg, JuMP.@NLexpression(pm.model, -(vr[p]-vr[n])*cig[idx]+(vi[p]-vi[n])*crg[idx]))
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

    var(pm, nw, :pg)[id] = pg
    var(pm, nw, :qg)[id] = qg

    if report
        _PMD.sol(pm, nw, :gen, id)[:pg] = pg
        _PMD.sol(pm, nw, :gen, id)[:qg] = qg
    end
end


"""
	function constraint_mc_generator_power_delta(
		pm::_PMD.AbstractNLExplicitNeutralIVRModel,
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
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)

    nph = length(pmin)

    vrg = Dict()
    vig = Dict()
    for (idx,c,d) in zip(1:nph, connections, [connections[2:end]..., connections[1]])
        vrg[idx] = JuMP.@NLexpression(pm.model, vr[c]-vr[d])
        vig[idx] = JuMP.@NLexpression(pm.model, vi[c]-vi[d])
    end

    pg = Vector{JuMP.NonlinearExpression}([])
    qg = Vector{JuMP.NonlinearExpression}([])
    for idx in 1:nph
        push!(pg, JuMP.@NLexpression(pm.model,  vrg[idx]*crg[idx]+vig[idx]*cig[idx]))
        push!(qg, JuMP.@NLexpression(pm.model, -vrg[idx]*cig[idx]+vig[idx]*crg[idx]))
    end

    JuMP.@NLconstraint(pm.model, [i in 1:nph], pmin[i] <= pg[i])
    JuMP.@NLconstraint(pm.model, [i in 1:nph], pmax[i] >= pg[i])
    JuMP.@NLconstraint(pm.model, [i in 1:nph], qmin[i] <= qg[i])
    JuMP.@NLconstraint(pm.model, [i in 1:nph], qmax[i] >= qg[i])

    var(pm, nw, :pg)[id] = JuMP.Containers.DenseAxisArray(pg, connections)
    var(pm, nw, :qg)[id] = JuMP.Containers.DenseAxisArray(qg, connections)

    if report
        _PMD.sol(pm, nw, :gen, id)[:pg] = pg
        _PMD.sol(pm, nw, :gen, id)[:qg] = qg
    end
end


"""
	function constraint_mc_generator_current_wye(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		connections::Vector{Int};
		report::Bool=true,
		bounded::Bool=true
	)

For IVR models with explicit neutrals,
creates expressions for the terminal current flows `:crg_bus` and `:cig_bus` of wye-connected generators
"""
function constraint_mc_generator_current_wye(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)
    var(pm, nw, :crg_bus)[id] = _PMD._merge_bus_flows(pm, [crg..., -sum(crg)], connections)
    var(pm, nw, :cig_bus)[id] = _PMD._merge_bus_flows(pm, [cig..., -sum(cig)], connections)
end


"""
	function constraint_mc_generator_current_delta(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		connections::Vector{Int};
		report::Bool=true,
		bounded::Bool=true
	)

For IVR models with explicit neutrals,
creates expressions for the terminal current flows `:crg_bus` and `:cig_bus` of delta-connected generators
"""
function constraint_mc_generator_current_delta(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crg = var(pm, nw, :crg, id)
    cig = var(pm, nw, :cig, id)
    Md = _get_delta_transformation_matrix(length(connections))
    var(pm, nw, :crg_bus)[id] = _PMD._merge_bus_flows(pm, Md'*crg, connections)
    var(pm, nw, :cig_bus)[id] = _PMD._merge_bus_flows(pm, Md'*cig, connections)
end


# TRANSFORMER - Constraints

"""
	function constraint_mc_transformer_voltage_yy(
		pm::_PMD.RectangularVoltageExplicitNeutralModels,
		nw::Int,
		trans_id::Int,
		f_bus::Int,
		t_bus::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int},
		pol::Int,
		tm_set::Vector{<:Real},
		tm_fixed::Vector{Bool},
		tm_scale::Real
	)

For rectangular voltage models with explicit neutrals,
links the voltage of the from-side and to-side transformer windings
for wye-wye connected transformers

```
(vr_fr_P-vr_fr_n) == scale * (vr_to_P.-vr_to_n)
(vi_fr_P-vi_fr_n) == scale * (vi_to_P.-vi_to_n)
```
"""
function constraint_mc_transformer_voltage_yy(pm::_PMD.RectangularVoltageExplicitNeutralModels, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vr_fr_P = [var(pm, nw, :vr, f_bus)[c] for c in f_connections[1:end-1]]
    vi_fr_P = [var(pm, nw, :vi, f_bus)[c] for c in f_connections[1:end-1]]
    vr_fr_n = var(pm, nw, :vr, f_bus)[f_connections[end]]
    vi_fr_n = var(pm, nw, :vi, f_bus)[f_connections[end]]
    vr_to_P = [var(pm, nw, :vr, t_bus)[c] for c in t_connections[1:end-1]]
    vi_to_P = [var(pm, nw, :vi, t_bus)[c] for c in t_connections[1:end-1]]
    vr_to_n = var(pm, nw, :vr, t_bus)[t_connections[end]]
    vi_to_n = var(pm, nw, :vi, t_bus)[t_connections[end]]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
    scale = (tm_scale*pol).*tm_set

    JuMP.@constraint(pm.model, (vr_fr_P.-vr_fr_n) .== scale.*(vr_to_P.-vr_to_n))
    JuMP.@constraint(pm.model, (vi_fr_P.-vi_fr_n) .== scale.*(vi_to_P.-vi_to_n))
end


"""
	function constraint_mc_transformer_voltage_dy(
		pm::_PMD.RectangularVoltageExplicitNeutralModels,
		nw::Int,
		trans_id::Int,
		f_bus::Int,
		t_bus::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int},
		pol::Int,
		tm_set::Vector{<:Real},
		tm_fixed::Vector{Bool},
		tm_scale::Real
	)

For rectangular voltage models with explicit neutrals,
links the voltage of the from-side and to-side transformer windings
for delta-wye connected transformers

```
Md*vr_fr_P == scale * (vr_to_P - vr_to_n)
Md*vi_fr_P == scale * (vi_to_P - vi_to_n)
```
"""
function constraint_mc_transformer_voltage_dy(pm::_PMD.RectangularVoltageExplicitNeutralModels, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    vr_fr_P = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr_P = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]
    vr_to_P = [var(pm, nw, :vr, t_bus)[c] for c in t_connections[1:end-1]]
    vi_to_P = [var(pm, nw, :vi, t_bus)[c] for c in t_connections[1:end-1]]
    vr_to_n = var(pm, nw, :vr, t_bus)[t_connections[end]]
    vi_to_n = var(pm, nw, :vi, t_bus)[t_connections[end]]

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
    scale = (tm_scale*pol).*tm_set

    n_phases = length(tm)
    Md = _get_delta_transformation_matrix(n_phases)

    JuMP.@constraint(pm.model, Md*vr_fr_P .== scale.*(vr_to_P .- vr_to_n))
    JuMP.@constraint(pm.model, Md*vi_fr_P .== scale.*(vi_to_P .- vi_to_n))
end


"""
	function constraint_mc_transformer_current_yy(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
		nw::Int,
		trans_id::Int,
		f_bus::Int,
		t_bus::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int},
		pol::Int,
		tm_set::Vector{<:Real},
		tm_fixed::Vector{Bool},
		tm_scale::Real
	)

For IVR models with explicit neutrals,
links the current variables of the from-side and to-side transformer windings,
and creates expressions for the terminal current flows
for wye-wye connected transformers

```
scale*cr_fr_P + cr_to_P == 0
scale*ci_fr_P + ci_to_P == 0
```
"""
function constraint_mc_transformer_current_yy(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    cr_fr_P = var(pm, nw, :crt, f_idx)
    ci_fr_P = var(pm, nw, :cit, f_idx)
    cr_to_P = var(pm, nw, :crt, t_idx)
    ci_to_P = var(pm, nw, :cit, t_idx)

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
    scale = (tm_scale*pol).*tm_set

    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ ci_to_P .== 0)

    var(pm, nw, :crt_bus)[f_idx] = _PMD._merge_bus_flows(pm, [cr_fr_P..., -sum(cr_fr_P)], f_connections)
    var(pm, nw, :cit_bus)[f_idx] = _PMD._merge_bus_flows(pm, [ci_fr_P..., -sum(ci_fr_P)], f_connections)
    var(pm, nw, :crt_bus)[t_idx] = _PMD._merge_bus_flows(pm, [cr_to_P..., -sum(cr_to_P)], t_connections)
    var(pm, nw, :cit_bus)[t_idx] = _PMD._merge_bus_flows(pm, [ci_to_P..., -sum(ci_to_P)], t_connections)
end


"""
	function constraint_mc_transformer_current_dy(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
		nw::Int,
		trans_id::Int,
		f_bus::Int,
		t_bus::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int},
		pol::Int,
		tm_set::Vector{<:Real},
		tm_fixed::Vector{Bool},
		tm_scale::Real
	)

For IVR models with explicit neutrals,
links the current variables of the from-side and to-side transformer windings,
and creates expressions for the terminal current flows
for delta-wye connected transformers

```
scale*cr_fr_P + cr_to_P == 0
scale*ci_fr_P + ci_to_P == 0
```
"""
function constraint_mc_transformer_current_dy(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, trans_id::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, pol::Int, tm_set::Vector{<:Real}, tm_fixed::Vector{Bool}, tm_scale::Real)
    cr_fr_P = var(pm, nw, :crt, f_idx)
    ci_fr_P = var(pm, nw, :cit, f_idx)
    cr_to_P = var(pm, nw, :crt, t_idx)
    ci_to_P = var(pm, nw, :cit, t_idx)

    # construct tm as a parameter or scaled variable depending on whether it is fixed or not
    tm = [tm_fixed[idx] ? tm_set[idx] : var(pm, nw, :tap, trans_id)[idx] for idx in 1:length(tm_fixed)]
    scale = (tm_scale*pol).*tm_set

    n_phases = length(tm)
    Md = _get_delta_transformation_matrix(n_phases)

    JuMP.@constraint(pm.model, scale.*cr_fr_P .+ cr_to_P .== 0)
    JuMP.@constraint(pm.model, scale.*ci_fr_P .+ ci_to_P .== 0)

    var(pm, nw, :crt_bus)[f_idx] = _PMD._merge_bus_flows(pm, Md'*cr_fr_P, f_connections)
    var(pm, nw, :cit_bus)[f_idx] = _PMD._merge_bus_flows(pm, Md'*ci_fr_P, f_connections)
    var(pm, nw, :crt_bus)[t_idx] = _PMD._merge_bus_flows(pm, [cr_to_P..., -sum(cr_to_P)], t_connections)
    var(pm, nw, :cit_bus)[t_idx] = _PMD._merge_bus_flows(pm, [ci_to_P..., -sum(ci_to_P)], t_connections)
end


"""
	function constraint_mc_transformer_thermal_limit(
		pm::_PMD.AbstractNLExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		f_idx::Tuple,
		t_idx::Tuple,
		f_bus::Int,
		t_bus::Int,
		f_connections::Vector,
		t_connections::Vector,
		config::ConnConfig,
		sm_ub::Real;
		report::Bool=true
	)

For non-linear IVR models with explicit neutrals,
imposes a bound on the magnitude of the total apparent power at both windings.
Expressions are created for the transformer power variables.

```
sum(pt_fr)^2 + sum(qt_fr)^2 <= sm_ub^2
sum(pt_to)^2 + sum(qt_to)^2 <= sm_ub^2
```
"""
function constraint_mc_transformer_thermal_limit(pm::_PMD.AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int, f_idx::Tuple, t_idx::Tuple, f_bus::Int, t_bus::Int, f_connections::Vector, t_connections::Vector, config, sm_ub::Real; report::Bool=true)
    vr_fr = var(pm, nw, :vr, f_bus)
    vi_fr = var(pm, nw, :vi, f_bus)
    vr_to = var(pm, nw, :vr, t_bus)
    vi_to = var(pm, nw, :vi, t_bus)

    crt_fr = var(pm, nw, :crt, f_idx)
    cit_fr = var(pm, nw, :cit, f_idx)
    crt_to = var(pm, nw, :crt, t_idx)
    cit_to = var(pm, nw, :cit, t_idx)

    if config==_PMD.WYE || length(crt_fr)==1
        P_fr = f_connections[1:end-1]
        n_fr = f_connections[end]
        vrt_fr = [vr_fr[p]-vr_fr[n_fr] for p in P_fr]
        vit_fr = [vi_fr[p]-vi_fr[n_fr] for p in P_fr]
    elseif config==_PMD.DELTA && length(crt_fr)==3
        M = _get_delta_transformation_matrix(3)
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
    pt_fr = JuMP.@NLexpression(pm.model, [i in idxs],  vrt_fr[i]*crt_fr[i] + vit_fr[i]*cit_fr[i])
    qt_fr = JuMP.@NLexpression(pm.model, [i in idxs], -vrt_fr[i]*cit_fr[i] + vit_fr[i]*crt_fr[i])
    pt_to = JuMP.@NLexpression(pm.model, [i in idxs],  vrt_to[i]*crt_to[i] + vit_to[i]*cit_to[i])
    qt_to = JuMP.@NLexpression(pm.model, [i in idxs], -vrt_to[i]*cit_to[i] + vit_to[i]*crt_to[i])

    if sm_ub<Inf
        JuMP.@NLconstraint(pm.model, sum(pt_fr[i] for i in idxs)^2 + sum(qt_fr[i] for i in idxs)^2 <= sm_ub^2)
        JuMP.@NLconstraint(pm.model, sum(pt_to[i] for i in idxs)^2 + sum(qt_to[i] for i in idxs)^2 <= sm_ub^2)
    end

    if report
        _PMD.sol(pm, nw, :transformer, id)[:pf] = pt_fr
        _PMD.sol(pm, nw, :transformer, id)[:qf] = qt_fr
        _PMD.sol(pm, nw, :transformer, id)[:pt] = pt_to
        _PMD.sol(pm, nw, :transformer, id)[:qt] = qt_to
        _PMD.sol(pm, nw, :transformer, id)[:smtot_fr] = JuMP.@NLexpression(pm.model, sqrt(sum(pt_fr[i] for i in idxs)^2 + sum(qt_fr[i] for i in idxs)^2))
        _PMD.sol(pm, nw, :transformer, id)[:smtot_to] = JuMP.@NLexpression(pm.model, sqrt(sum(pt_to[i] for i in idxs)^2 + sum(qt_to[i] for i in idxs)^2))
    end
end


"""
	function constraint_mc_switch_current(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		f_idx::Tuple{Int,Int,Int},
		t_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int};
		report::Bool=true
	)

For IVR models with explicit neutrals,
create expressions for the terminal current flows `:crsw_bus` and `cisw_bus`,
and link the from-side to the to-side switch current
"""
function constraint_mc_switch_current(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, id::Int, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}; report::Bool=true)
    crsw_fr = var(pm, nw, :crsw, f_idx)
    cisw_fr = var(pm, nw, :cisw, f_idx)
    crsw_to = var(pm, nw, :crsw, t_idx)
    cisw_to = var(pm, nw, :cisw, t_idx)

    JuMP.@constraint(pm.model, crsw_fr .+ crsw_to .== 0)
    JuMP.@constraint(pm.model, cisw_fr .+ cisw_to .== 0)

    var(pm, nw, :crsw_bus)[f_idx] = _PMD._merge_bus_flows(pm, crsw_fr, f_connections)
    var(pm, nw, :cisw_bus)[f_idx] = _PMD._merge_bus_flows(pm, cisw_fr, f_connections)
    var(pm, nw, :crsw_bus)[t_idx] = _PMD._merge_bus_flows(pm, crsw_to, t_connections)
    var(pm, nw, :cisw_bus)[t_idx] = _PMD._merge_bus_flows(pm, cisw_to, t_connections)
end


"""
	function constraint_mc_switch_current_limit(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		connections::Vector{Int},
		rating::Vector{<:Real}
	)

For IVR models with explicit neutrals,
imposes a bound on the switch current magnitude per conductor.
Note that a bound on the from-side implies the same bound on the to-side current,
so it suffices to apply this only explicitly at the from-side.
"""
function constraint_mc_switch_current_limit(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, connections::Vector{Int}, rating::Vector{<:Real})::Nothing
    crsw = var(pm, nw, :crsw, f_idx)
    cisw = var(pm, nw, :cisw, f_idx)

    mu_cm_fr = JuMP.ConstraintRef[]
    for idx in 1:length(rating)
        if rating[idx] < Inf
            push!(mu_cm_fr, JuMP.@constraint(pm.model, crsw[idx]^2 + cisw[idx]^2 <= rating[idx]^2))
        end
    end

    _PMD.con(pm, nw, :mu_cm_switch)[f_idx] = mu_cm_fr

    nothing
end


"""
	function constraint_mc_switch_thermal_limit(
		pm::_PMD.AbstractNLExplicitNeutralIVRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		rating::Vector{<:Real}
	)

For IVR models with explicit neutrals,
imposes a bound on the switch power magnitude per conductor.
Note that a bound on the from-side implies the same bound on the to-side power
when the switch is closed (equal voltages), and also when it is open since the
power then equals zero on both ends.
"""
function constraint_mc_switch_thermal_limit(pm::_PMD.AbstractNLExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})::Nothing
    vr_fr = [var(pm, nw, :vr, f_idx[2])[t] for t in f_connections]
    vi_fr = [var(pm, nw, :vi, f_idx[2])[t] for t in f_connections]
    crsw_fr = var(pm, nw, :crsw, f_idx)
    cisw_fr = var(pm, nw, :cisw, f_idx)

    mu_sm_fr = JuMP.ConstraintRef[]
    for idx in 1:length(rating)
        if rating[idx] < Inf
            psw_fr_idx = JuMP.@NLexpression(pm.model,  vr_fr[idx]*crsw_fr[idx] + vi_fr[idx]*cisw_fr[idx])
            qsw_fr_idx = JuMP.@NLexpression(pm.model, -vr_fr[idx]*cisw_fr[idx] + vi_fr[idx]*crsw_fr[idx])

            push!(mu_sm_fr, JuMP.@NLconstraint(pm.model, psw_fr_idx^2 + qsw_fr_idx^2 <= rating[idx]^2))
        end
    end

    _PMD.con(pm, nw, :mu_sm_switch)[f_idx] = mu_sm_fr

    nothing
end


"""
	function constraint_mc_load_current_wye(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
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
function constraint_mc_load_current_wye(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)

    crd = Vector{JuMP.NonlinearExpression}([])
    cid = Vector{JuMP.NonlinearExpression}([])

    phases = connections[1:end-1]
    n      = connections[end]

    for (idx, p) in enumerate(phases)
        push!(crd, JuMP.@NLexpression(pm.model,
             a[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1)
            +b[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1)
        ))
        push!(cid, JuMP.@NLexpression(pm.model,
             a[idx]*(vi[p]-vi[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(alpha[idx]/2-1)
            -b[idx]*(vr[p]-vr[n])*((vr[p]-vr[n])^2+(vi[p]-vi[n])^2)^(beta[idx]/2 -1)
        ))
    end

    var(pm, nw, :crd)[id] = crd
    var(pm, nw, :cid)[id] = cid

    crd_bus_n = JuMP.@NLexpression(pm.model, -sum(crd[i] for i in 1:length(phases)))
    cid_bus_n = JuMP.@NLexpression(pm.model, -sum(cid[i] for i in 1:length(phases)))

    var(pm, nw, :crd_bus)[id] = crd_bus = _PMD._merge_bus_flows(pm, [crd..., crd_bus_n], connections)
    var(pm, nw, :cid_bus)[id] = cid_bus = _PMD._merge_bus_flows(pm, [cid..., cid_bus_n], connections)

    if report
        pd_bus = Vector{JuMP.NonlinearExpression}([])
        qd_bus = Vector{JuMP.NonlinearExpression}([])
        for (idx,c) in enumerate(connections)
            push!(pd_bus, JuMP.@NLexpression(pm.model,  vr[c]*crd_bus[c]+vi[c]*cid_bus[c]))
            push!(qd_bus, JuMP.@NLexpression(pm.model, -vr[c]*cid_bus[c]+vi[c]*crd_bus[c]))
        end

        _PMD.sol(pm, nw, :load, id)[:pd_bus] = JuMP.Containers.DenseAxisArray(pd_bus, connections)
        _PMD.sol(pm, nw, :load, id)[:qd_bus] = JuMP.Containers.DenseAxisArray(qd_bus, connections)

        _PMD.sol(pm, nw, :load, id)[:crd] = JuMP.Containers.DenseAxisArray(crd, connections)
        _PMD.sol(pm, nw, :load, id)[:cid] = JuMP.Containers.DenseAxisArray(cid, connections)

        _PMD.sol(pm, nw, :load, id)[:crd_bus] = crd_bus
        _PMD.sol(pm, nw, :load, id)[:cid_bus] = cid_bus

        pd = Vector{JuMP.NonlinearExpression}([])
        qd = Vector{JuMP.NonlinearExpression}([])
        for (idx, p) in enumerate(phases)
            push!(pd, JuMP.@NLexpression(pm.model, a[idx]*(vr[p]^2+vi[p]^2)^(alpha[idx]/2) ))
            push!(qd, JuMP.@NLexpression(pm.model, b[idx]*(vr[p]^2+vi[p]^2)^(beta[idx]/2)  ))
        end
        _PMD.sol(pm, nw, :load, id)[:pd] = JuMP.Containers.DenseAxisArray(pd, connections)
        _PMD.sol(pm, nw, :load, id)[:qd] = JuMP.Containers.DenseAxisArray(qd, connections)
    end
end


"""
	function constraint_mc_load_current_delta(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
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
function constraint_mc_load_current_delta(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}; report::Bool=true)
    vr = var(pm, nw, :vr, bus_id)
    vi = var(pm, nw, :vi, bus_id)


    ph = connections
    ph_next = [connections[2:end]..., connections[1]]
    P = length(ph)
    idxs = 1:P
    idxs_prev = [idxs[end], idxs[1:end-1]...]

    vrd = [vr[c]-vr[d] for (c,d) in zip(ph,ph_next)]
    vid = [vi[c]-vi[d] for (c,d) in zip(ph,ph_next)]

    crd = JuMP.@NLexpression(pm.model, [i in 1:P],
        a[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1)
       +b[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1)
    )
    cid = JuMP.@NLexpression(pm.model, [i in 1:P],
        a[i]*vid[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2-1)
       -b[i]*vrd[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2 -1)
    )

    crd_bus = JuMP.@NLexpression(pm.model, [i in 1:P], crd[i]-crd[idxs_prev[i]])
    cid_bus = JuMP.@NLexpression(pm.model, [i in 1:P], cid[i]-cid[idxs_prev[i]])

    var(pm, nw, :crd_bus)[id] = _PMD._merge_bus_flows(pm, crd_bus, connections)
    var(pm, nw, :cid_bus)[id] = _PMD._merge_bus_flows(pm, cid_bus, connections)

    if report
        pd_bus = JuMP.@NLexpression(pm.model, [i in 1:P],  vr[i]*crd_bus[i]+vi[i]*cid_bus[i])
        qd_bus = JuMP.@NLexpression(pm.model, [i in 1:P], -vr[i]*cid_bus[i]+vi[i]*crd_bus[i])

        _PMD.sol(pm, nw, :load, id)[:pd_bus] = pd_bus
        _PMD.sol(pm, nw, :load, id)[:qd_bus] = qd_bus

        pd = JuMP.@NLexpression(pm.model, [i in 1:P], a[i]*(vrd[i]^2+vid[i]^2)^(alpha[i]/2) )
        qd = JuMP.@NLexpression(pm.model, [i in 1:P], b[i]*(vrd[i]^2+vid[i]^2)^(beta[i]/2)  )
        _PMD.sol(pm, nw, :load, id)[:pd] = pd
        _PMD.sol(pm, nw, :load, id)[:qd] = qd
    end
end


# BRANCH - Constraints

"""
	function constraint_mc_current_from(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
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
function constraint_mc_current_from(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, f_bus::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, g_sh_fr::Matrix{<:Real}, b_sh_fr::Matrix{<:Real}; report::Bool=true)
    vr_fr = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]

    cr_fr =  var(pm, nw, :cr, f_idx)
    ci_fr =  var(pm, nw, :ci, f_idx)

    csr_fr =  var(pm, nw, :csr, f_idx[1])
    csi_fr =  var(pm, nw, :csi, f_idx[1])

    JuMP.@constraint(pm.model, cr_fr .== csr_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr)
    JuMP.@constraint(pm.model, ci_fr .== csi_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr)

    var(pm, nw, :cr_bus)[f_idx] = cr_bus_fr = _PMD._merge_bus_flows(pm, cr_fr, f_connections)
    var(pm, nw, :ci_bus)[f_idx] = ci_bus_fr = _PMD._merge_bus_flows(pm, ci_fr, f_connections)

    if report
        _PMD.sol(pm, nw, :branch, f_idx[1])[:pf] =  cr_fr.*vr_fr .+ ci_fr.*vi_fr
        _PMD.sol(pm, nw, :branch, f_idx[1])[:qf] = -cr_fr.*vi_fr .+ ci_fr.*vr_fr
    end

end


"""
	function constraint_mc_current_to(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
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
function constraint_mc_current_to(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, t_bus, f_idx::Tuple{Int,Int,Int}, t_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, g_sh_to::Matrix{<:Real}, b_sh_to::Matrix{<:Real}; report::Bool=true)
    vr_to = [var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[c] for c in t_connections]

    cr_to = var(pm, nw, :cr, t_idx)
    ci_to = var(pm, nw, :ci, t_idx)

    csr_to = -var(pm, nw, :csr, f_idx[1])
    csi_to = -var(pm, nw, :csi, f_idx[1])

    JuMP.@constraint(pm.model, cr_to .== csr_to + g_sh_to*vr_to - b_sh_to*vi_to)
    JuMP.@constraint(pm.model, ci_to .== csi_to + g_sh_to*vi_to + b_sh_to*vr_to)

    var(pm, nw, :cr_bus)[t_idx] = cr_bus_to = _PMD._merge_bus_flows(pm, cr_to, t_connections)
    var(pm, nw, :ci_bus)[t_idx] = ci_bus_to = _PMD._merge_bus_flows(pm, ci_to, t_connections)

    if report
        _PMD.sol(pm, nw, :branch, t_idx[1])[:pt] =  cr_to.*vr_to .+ ci_to.*vi_to
        _PMD.sol(pm, nw, :branch, t_idx[1])[:qt] = -cr_to.*vi_to .+ ci_to.*vr_to
    end
end


"""
	function constraint_mc_bus_voltage_drop(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
		nw::Int,
		i::Int,
		f_bus::Int,
		t_bus::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		t_connections::Vector{Int},
		r::Matrix{<:Real},
		x::Matrix{<:Real}
	)

For IVR models with explicit neutrals,
defines voltage drop over a branch, linking from and to side complex voltage.

```
vr_to == vr_fr - r*csr_fr + x*csi_fr
vi_to == vi_fr - r*csi_fr - x*csr_fr
```
"""
function constraint_mc_bus_voltage_drop(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, r::Matrix{<:Real}, x::Matrix{<:Real})
    vr_fr = [var(pm, nw, :vr, f_bus)[c] for c in f_connections]
    vi_fr = [var(pm, nw, :vi, f_bus)[c] for c in f_connections]

    vr_to = [var(pm, nw, :vr, t_bus)[c] for c in t_connections]
    vi_to = [var(pm, nw, :vi, t_bus)[c] for c in t_connections]

    csr_fr = var(pm, nw, :csr, f_idx[1])
    csi_fr = var(pm, nw, :csi, f_idx[1])

    JuMP.@constraint(pm.model, vr_to .== vr_fr - r*csr_fr + x*csi_fr)
    JuMP.@constraint(pm.model, vi_to .== vi_fr - r*csi_fr - x*csr_fr)
end


"""
	function constraint_mc_branch_current_limit(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
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
    cr_fr = var(pm, nw, :cr, f_idx)
    ci_fr = var(pm, nw, :ci, f_idx)
    cr_to = var(pm, nw, :cr, t_idx)
    ci_to = var(pm, nw, :ci, t_idx)

    cnds_finite_rating = [c for (c,r) in enumerate(c_rating) if r<Inf]
    JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_fr[c]^2+ci_fr[c]^2 <= c_rating[c]^2)
    JuMP.@constraint(pm.model, [c in cnds_finite_rating], cr_to[c]^2+ci_to[c]^2 <= c_rating[c]^2)
end


"""
	function constraint_mc_thermal_limit_from(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
		nw::Int,
		f_idx::Tuple{Int,Int,Int},
		f_connections::Vector{Int},
		rate_a::Vector{<:Real}
	)

For IVR models with explicit neutrals,
imposes a bound on the from-side line power magnitude.
"""
function constraint_mc_thermal_limit_from(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rate_a::Vector{<:Real})
    vr_fr = [var(pm, nw, :vr, f_idx[1])[t] for t in f_connections]
    vi_fr = [var(pm, nw, :vi, f_idx[1])[t] for t in f_connections]
    cr_fr = var(pm, nw, :cr, f_idx)
    ci_fr = var(pm, nw, :ci, f_idx)

    for idx in 1:length(rate_a)
        if rate_a[idx]<Inf
            pf_idx = JuMP.@NLexpression(pm.model,  vr_fr[idx]*cr_fr[idx] + vi_fr[idx]*ci_fr[idx])
            qf_idx = JuMP.@NLexpression(pm.model, -vr_fr[idx]*ci_fr[idx] + vi_fr[idx]*cr_fr[idx])

            JuMP.@NLconstraint(pm.model, pf_idx^2 + qf_idx^2 <= rate_a[idx]^2)
        end
    end
end


"""
	function constraint_mc_thermal_limit_to(
		pm::_PMD.AbstractExplicitNeutralIVRModel,
		nw::Int,
		t_idx::Tuple{Int,Int,Int},
		t_connections::Vector{Int},
		rate_a::Vector{<:Real}
	)

For IVR models with explicit neutrals,
imposes a bound on the to-side line power magnitude.
"""
function constraint_mc_thermal_limit_to(pm::_PMD.AbstractExplicitNeutralIVRModel, nw::Int, t_idx::Tuple{Int,Int,Int}, t_connections::Vector{Int}, rate_a::Vector{<:Real})
    vr_to = [var(pm, nw, :vr, t_idx[1])[t] for t in t_connections]
    vi_to = [var(pm, nw, :vi, t_idx[1])[t] for t in t_connections]
    cr_to = var(pm, nw, :cr, t_idx)
    ci_to = var(pm, nw, :ci, t_idx)

    for idx in 1:length(rate_a)
        if rate_a[idx]<Inf
            pt_idx = JuMP.@NLexpression(pm.model,  vr_to[idx]*cr_to[idx] + vi_to[idx]*ci_to[idx])
            qt_idx = JuMP.@NLexpression(pm.model, -vr_to[idx]*ci_to[idx] + vi_to[idx]*cr_to[idx])

            JuMP.@NLconstraint(pm.model, pt_idx^2 + qt_idx^2 <= rate_a[idx]^2)
        end
    end
end


"""
	function constraint_mc_current_balance(
		pm::_PMD.RectangularVoltageExplicitNeutralModels,
		nw::Int,
		i::Int,
		terminals::Vector{Int},
		grounded::Vector{Bool},
		bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
		bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
		bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}},
		bus_gens::Vector{Tuple{Int,Vector{Int}}},
		bus_storage::Vector{Tuple{Int,Vector{Int}}},
		bus_loads::Vector{Tuple{Int,Vector{Int}}},
		bus_shunts::Vector{Tuple{Int,Vector{Int}}}
	)

Kirchhoff's current law applied to buses
`sum(cr + im*ci) = 0`
"""
function constraint_mc_current_balance(pm::_PMD.RectangularVoltageExplicitNeutralModels, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = var(pm, nw, :vr, i)
    vi = var(pm, nw, :vi, i)

    cr    = get(var(pm, nw), :cr_bus,   Dict()); _PMD._check_var_keys(cr, bus_arcs, "real current", "branch")
    ci    = get(var(pm, nw), :ci_bus,   Dict()); _PMD._check_var_keys(ci, bus_arcs, "imaginary current", "branch")
    crd   = get(var(pm, nw), :crd_bus,  Dict()); _PMD._check_var_keys(crd, bus_loads, "real current", "load")
    cid   = get(var(pm, nw), :cid_bus,  Dict()); _PMD._check_var_keys(cid, bus_loads, "imaginary current", "load")
    crg   = get(var(pm, nw), :crg_bus,  Dict()); _PMD._check_var_keys(crg, bus_gens, "real current", "generator")
    cig   = get(var(pm, nw), :cig_bus,  Dict()); _PMD._check_var_keys(cig, bus_gens, "imaginary current", "generator")
    crs   = get(var(pm, nw), :crs_bus,  Dict()); _PMD._check_var_keys(crs, bus_storage, "real currentr", "storage")
    cis   = get(var(pm, nw), :cis_bus,  Dict()); _PMD._check_var_keys(cis, bus_storage, "imaginary current", "storage")
    crsw  = get(var(pm, nw), :crsw_bus, Dict()); _PMD._check_var_keys(crsw, bus_arcs_sw, "real current", "switch")
    cisw  = get(var(pm, nw), :cisw_bus, Dict()); _PMD._check_var_keys(cisw, bus_arcs_sw, "imaginary current", "switch")
    crt   = get(var(pm, nw), :crt_bus,  Dict()); _PMD._check_var_keys(crt, bus_arcs_trans, "real current", "transformer")
    cit   = get(var(pm, nw), :cit_bus,  Dict()); _PMD._check_var_keys(cit, bus_arcs_trans, "imaginary current", "transformer")

    Gt, Bt = _PMD._build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx, t) in ungrounded_terminals
        @smart_constraint(pm.model,  [cr, crd, crg, crs, crsw, crt, vr],
                                      sum(cr[a][t] for (a, conns) in bus_arcs if t in conns)
                                    + sum(crsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                    + sum(crt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                    ==
                                      sum(crg[g][t]         for (g, conns) in bus_gens if t in conns)
                                    - sum(crs[s][t]         for (s, conns) in bus_storage if t in conns)
                                    - sum(crd[d][t]         for (d, conns) in bus_loads if t in conns)
                                    - sum( Gt[idx,jdx]*vr[u] -Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals) # shunts
                                    )
        @smart_constraint(pm.model, [ci, cid, cig, cis, cisw, cit, vi],
                                      sum(ci[a][t] for (a, conns) in bus_arcs if t in conns)
                                    + sum(cisw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                    + sum(cit[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                    ==
                                      sum(cig[g][t]         for (g, conns) in bus_gens if t in conns)
                                    - sum(cis[s][t]         for (s, conns) in bus_storage if t in conns)
                                    - sum(cid[d][t]         for (d, conns) in bus_loads if t in conns)
                                    - sum( Gt[idx,jdx]*vi[u] +Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals) # shunts
                                    )
    end
end