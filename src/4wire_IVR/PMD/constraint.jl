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
    Md = _PMD._get_delta_transformation_matrix(n_phases)

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
    Md = _PMD._get_delta_transformation_matrix(n_phases)

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