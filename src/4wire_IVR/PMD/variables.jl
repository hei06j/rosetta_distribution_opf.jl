# TRANSFORMER

# TRANSFORMER - Variables

"""
	function variable_mc_transformer_current(
		pm::_PMD.AbstractExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

For IVR models with explicit neutrals,
create transformer current variables `:crt` and `:cit`,
and placeholder dictionaries for the terminal current flows `:crt_bus` and `:cit_bus`
"""
function variable_mc_transformer_current(pm::_PMD.AbstractExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_transformer_current_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_transformer_current_imaginary(pm; nw=nw, bounded=bounded, report=report)

    var(pm, nw)[:crt_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:cit_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end


"""
	function variable_mc_transformer_current_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates transformer real current variables `:crt` for models with explicit neutrals
"""
function variable_mc_transformer_current_real(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in _PMD.ref(pm, nw, :transformer))
    crt = var(pm, nw)[:crt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:int_dim[l]], base_name="$(nw)_crt_$((l,i,j))",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :transformer, l), "crt_start", c, 0.0)
        ) for (l,i,j) in _PMD.ref(pm, nw, :arcs_transformer)
    )

    report && _IM.sol_component_value_edge(pm, _PMD.pmd_it_sym, nw, :transformer, :cr_fr, :cr_to, _PMD.ref(pm, nw, :arcs_transformer_from), _PMD.ref(pm, nw, :arcs_transformer_to), crt)
end


"""
	function variable_mc_transformer_current_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates transformer imaginary current variables `:cit` for models with explicit neutrals
"""
function variable_mc_transformer_current_imaginary(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in _PMD.ref(pm, nw, :transformer))
    cit = var(pm, nw)[:cit] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:int_dim[l]], base_name="$(nw)_cit_$((l,i,j))",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :transformer, l), "cit_start", c, 0.0)
        ) for (l,i,j) in _PMD.ref(pm, nw, :arcs_transformer)
    )

    report && _IM.sol_component_value_edge(pm, _PMD.pmd_it_sym, nw, :transformer, :ci_fr, :ci_to, _PMD.ref(pm, nw, :arcs_transformer_from), _PMD.ref(pm, nw, :arcs_transformer_to), cit)
end


"""
	function variable_mc_transformer_power_real(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates transformer active power variables `:pt` for models with explicit neutrals
"""
function variable_mc_transformer_power_real(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in _PMD.ref(pm, nw, :transformer))
    pt = var(pm, nw)[:pt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:int_dim[l]], base_name="$(nw)_pt_$((l,i,j))",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :transformer, l), "pt_start", c, 0.0)
        ) for (l,i,j) in _PMD.ref(pm, nw, :arcs_transformer)
    )

    if bounded
        for (l,i,j) in _PMD.ref(pm, nw, :arcs_transformer_from)
            trans = _PMD.ref(pm, nw, :transformer, l)
            f_bus = _PMD.ref(pm, nw, :bus, i)
            t_bus = _PMD.ref(pm, nw, :bus, j)
            sm_ub = trans["sm_ub"]
            set_lower_bound(pt[(l,i,j)], -sm_ub)
            set_upper_bound(pt[(l,i,j)],  sm_ub)
            set_lower_bound(pt[(l,j,i)], -sm_ub)
            set_upper_bound(pt[(l,j,i)],  sm_ub)
        end
    end

    report && _IM.sol_component_value_edge(pm, _PMD.pmd_it_sym, nw, :transformer, :pf, :pt, _PMD.ref(pm, nw, :arcs_transformer_from), _PMD.ref(pm, nw, :arcs_transformer_to), pt)
end


"""
	function variable_mc_transformer_power_imaginary(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates transformer reactive power variables `:qt` for models with explicit neutrals
"""
function variable_mc_transformer_power_imaginary(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(l => _infer_int_dim_transformer(trans, false) for (l,trans) in _PMD.ref(pm, nw, :transformer))
    qt = var(pm, nw)[:qt] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:int_dim[l]], base_name="$(nw)_qt_$((l,i,j))",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :transformer, l), "qt_start", c, 0.0)
        ) for (l,i,j) in _PMD.ref(pm, nw, :arcs_transformer)
    )

    if bounded
        for (l,i,j) in _PMD.ref(pm, nw, :arcs_transformer_from)
            trans = _PMD.ref(pm, nw, :transformer, l)
            f_bus = _PMD.ref(pm, nw, :bus, i)
            t_bus = _PMD.ref(pm, nw, :bus, j)
            sm_ub = trans["sm_ub"]
            set_lower_bound(qt[(l,i,j)], -sm_ub)
            set_upper_bound(qt[(l,i,j)],  sm_ub)
            set_lower_bound(qt[(l,j,i)], -sm_ub)
            set_upper_bound(qt[(l,j,i)],  sm_ub)
        end
    end

    report && _IM.sol_component_value_edge(pm, _PMD.pmd_it_sym, nw, :transformer, :qf, :qt, _PMD.ref(pm, nw, :arcs_transformer_from), _PMD.ref(pm, nw, :arcs_transformer_to), qt)
end


# TRANSFORMER -  Variable - Non-linear

"""
	function variable_mc_transformer_power(
		pm::_PMD.AbstractNLExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

For non-linear IVR models with explicit neutrals,
no power variables are required.
"""
function variable_mc_transformer_power(pm::_PMD.AbstractNLExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # do nothing
end


# SWITCH

# SWITCH - Variables

"""
	function variable_mc_switch_current(
		pm::_PMD.AbstractExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

For IVR models with explicit neutrals,
creates switch current variables `:crs` and `:cis`,
and placeholder dictionaries for the terminal current flows `:crsw_bus` and `:cisw_bus`
"""
function variable_mc_switch_current(pm::_PMD.AbstractExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_switch_current_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_switch_current_imaginary(pm; nw=nw, bounded=bounded, report=report)

    var(pm, nw)[:crsw_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:cisw_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end


"""
	function variable_mc_switch_current_real(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For models with explicit neutrals,
creates switch real current variables `:crsw` for models with explicit neutrals.
"""
function variable_mc_switch_current_real(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(switch["f_connections"]) for (l,switch) in _PMD.ref(pm, nw, :switch))
    crsw = var(pm, nw)[:crsw] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_crsw_$((l,i,j))",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :switch, l), "crsw_start", c, 0.0)
        ) for (l,i,j) in _PMD.ref(pm, nw, :arcs_switch)
    )

    if bounded
        for (l,i,j) in _PMD.ref(pm, nw, :arcs_switch)
            cmax = _PMD.ref(pm, nw, :switch, l)["current_rating"]
            set_upper_bound.(crsw[(l,i,j)],  cmax)
            set_lower_bound.(crsw[(l,i,j)], -cmax)
        end
    end

    report && _IM.sol_component_value_edge(pm, _PMD.pmd_it_sym, nw, :switch, :cr_fr, :cr_to, _PMD.ref(pm, nw, :arcs_switch_from), _PMD.ref(pm, nw, :arcs_switch_to), crsw)
end


"""
	function variable_mc_switch_current_imaginary(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For models with explicit neutrals,
creates switch imaginary current variables `:cisw` for models with explicit neutrals.
"""
function variable_mc_switch_current_imaginary(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(switch["f_connections"]) for (l,switch) in _PMD.ref(pm, nw, :switch))
    cisw = var(pm, nw)[:cisw] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_cisw_$((l,i,j))",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :switch, l), "cisw_start", c, 0.0)
        ) for (l,i,j) in _PMD.ref(pm, nw, :arcs_switch)
    )

    if bounded
        for (l,i,j) in _PMD.ref(pm, nw, :arcs_switch)
            cmax = _PMD.ref(pm, nw, :switch, l)["current_rating"]
            set_upper_bound.(cisw[(l,i,j)],  cmax)
            set_lower_bound.(cisw[(l,i,j)], -cmax)
        end
    end

    report && _IM.sol_component_value_edge(pm, _PMD.pmd_it_sym, nw, :switch, :ci_fr, :ci_to, _PMD.ref(pm, nw, :arcs_switch_from), _PMD.ref(pm, nw, :arcs_switch_to), cisw)
end