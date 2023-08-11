# BUS

# BUS - Variables

"""
	function variable_mc_bus_voltage(
		pm::_PMD.RectangularVoltageExplicitNeutralModels;
		nw=nw_id_default,
		bounded::Bool=true,
	)

Creates rectangular voltage variables `:vr` and `:vi` for models with explicit neutrals
"""
function variable_mc_bus_voltage(pm::_PMD.RectangularVoltageExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_bus_voltage_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_bus_voltage_imaginary(pm; nw=nw, bounded=bounded, report=report)
end


"""
	function variable_mc_bus_voltage_real(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,å
		bounded::Bool=true,
		report::Bool=true
	)

Creates real voltage variables `:vr` for models with explicit neutrals
"""
function variable_mc_bus_voltage_real(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in _PMD.ref(pm, nw, :bus))

    vr = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vr_$(i)",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :bus, i), "vr_start", t, 0.0)
        ) for i in _PMD.ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in _PMD.ref(pm, nw, :bus)
            for (idx,t) in enumerate(terminals[i])
                set_lower_bound.(vr[i][t], -bus["vmax"][idx])
                set_upper_bound.(vr[i][t],  bus["vmax"][idx])
            end
        end
    end

    # perfectly grounded terminals get a constant 0 instead of a variable
    vr = var(pm, nw)[:vr] = Dict(i=>JuMP.Containers.DenseAxisArray(
        Vector{JuMP.AffExpr}([t in vr[i].axes[1] ? vr[i][t] : 0.0 for t in bus["terminals"]]), bus["terminals"]
        ) for (i, bus) in _PMD.ref(pm, nw, :bus))

    report && _IM.sol_component_value(pm, _PMD.pmd_it_sym, nw, :bus, :vr, _PMD.ids(pm, nw, :bus), vr)
end


"""
	function variable_mc_bus_voltage_imaginary(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates imaginary voltage variables `:vr` for models with explicit neutrals
"""
function variable_mc_bus_voltage_imaginary(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i,bus) in _PMD.ref(pm, nw, :bus))
    vi = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vi_$(i)",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :bus, i), "vi_start", t, 0.0)
        ) for i in _PMD.ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in _PMD.ref(pm, nw, :bus)
            for (idx,t) in enumerate(terminals[i])
                set_lower_bound.(vi[i][t], -bus["vmax"][idx])
                set_upper_bound.(vi[i][t],  bus["vmax"][idx])
            end
        end
    end

    # perfectly grounded terminals get a constant 0 instead of a variable
    vi = var(pm, nw)[:vi] = Dict(i=>JuMP.Containers.DenseAxisArray(
        Vector{JuMP.AffExpr}([t in vi[i].axes[1] ? vi[i][t] : 0.0 for t in bus["terminals"]]), bus["terminals"]
        ) for (i, bus) in _PMD.ref(pm, nw, :bus))

    report && _IM.sol_component_value(pm, _PMD.pmd_it_sym, nw, :bus, :vi, _PMD.ids(pm, nw, :bus), vi)
end



# BRANCH

# BRANCH - Variables

"""
	function variable_mc_branch_current_real(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates branch real current variables `:cr` for models with explicit neutrals
"""
function variable_mc_branch_current_real(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in _PMD.ref(pm, nw, :branch))
    cr = var(pm, nw)[:cr] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_cr_$((l,i,j))",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :branch, l), "cr_start", c, 0.0)
        ) for (l,i,j) in _PMD.ref(pm, nw, :arcs_branch)
    )

    if bounded
        for (l,i,j) in _PMD.ref(pm, nw, :arcs_branch)
            cmax = _PMD.ref(pm, nw, :branch, l)["c_rating_a"]
            set_upper_bound.(cr[(l,i,j)],  cmax)
            set_lower_bound.(cr[(l,i,j)], -cmax)
        end
    end

    report && _IM.sol_component_value_edge(pm, _PMD.pmd_it_sym, nw, :branch, :cr_fr, :cr_to, _PMD.ref(pm, nw, :arcs_branch_from), _PMD.ref(pm, nw, :arcs_branch_to), cr)
end


"""
	function variable_mc_branch_current_imaginary(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates branch imaginary current variables `:ci` for models with explicit neutrals
"""
function variable_mc_branch_current_imaginary(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in _PMD.ref(pm, nw, :branch))
    ci = var(pm, nw)[:ci] = Dict((l,i,j) => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_ci_$((l,i,j))",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :branch, l), "ci_start", c, 0.0)
        ) for (l,i,j) in _PMD.ref(pm, nw, :arcs_branch)
    )

    if bounded
        for (l,i,j) in _PMD.ref(pm, nw, :arcs_branch)
            cmax = _PMD.ref(pm, nw, :branch, l)["c_rating_a"]
            set_upper_bound.(ci[(l,i,j)],  cmax)
            set_lower_bound.(ci[(l,i,j)], -cmax)
        end
    end

    report && _IM.sol_component_value_edge(pm, _PMD.pmd_it_sym, nw, :branch, :ci_fr, :ci_to, _PMD.ref(pm, nw, :arcs_branch_from), _PMD.ref(pm, nw, :arcs_branch_to), ci)
end


"""
	function variable_mc_branch_current_series_real(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates branch real series current variables `:csr` for models with explicit neutrals
"""
function variable_mc_branch_current_series_real(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in _PMD.ref(pm, nw, :branch))
    csr = var(pm, nw)[:csr] = Dict(l => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_csr_$(l)",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :branch, l), "csr_start", c, 0.0)
        ) for (l,i,j) in _PMD.ref(pm, nw, :arcs_branch)
    )

    report && _IM.sol_component_value(pm, _PMD.pmd_it_sym, nw, :branch, :csr_fr, _PMD.ids(pm, nw, :branch), csr)
end


"""
	function variable_mc_branch_current_series_imaginary(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates branch imaginary series current variables `:csi` for models with explicit neutrals
"""
function variable_mc_branch_current_series_imaginary(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in _PMD.ref(pm, nw, :branch))
    csi = var(pm, nw)[:csi] = Dict(l => JuMP.@variable(pm.model,
            [c in 1:nconds[l]], base_name="$(nw)_csi_$(l)",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :branch, l), "csi_start", c, 0.0)
        ) for (l,i,j) in _PMD.ref(pm, nw, :arcs_branch)
    )

    report && _IM.sol_component_value(pm, _PMD.pmd_it_sym, nw, :branch, :csi_fr, _PMD.ids(pm, nw, :branch), csi)
end


"""
	function variable_mc_branch_current(
		pm::_PMD.AbstractExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For IVR models with explicit neutrals,
creates total current variables `:cr` and `:ci`,
series current variables `:csr` and `:csi`,
and placeholder dictionaries for the terminal current flows `:cr_bus` and `:ci_bus`
"""
function variable_mc_branch_current(pm::_PMD.AbstractExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_branch_current_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_branch_current_imaginary(pm; nw=nw, bounded=bounded, report=report)

    variable_mc_branch_current_series_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_branch_current_series_imaginary(pm; nw=nw, bounded=bounded, report=report)

    var(pm, nw)[:cr_bus] = Dict{Tuple{Int,Int,Int}, Any}()
    var(pm, nw)[:ci_bus] = Dict{Tuple{Int,Int,Int}, Any}()
end




# LOAD - Variables

"""
	function variable_mc_load_current(
		pm::_PMD.AbstractExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For IVR models with explicit neutrals,
creates placeholder dictionaries for the load current `:crd` and `:cid`,
and for the terminal current flows `:crd_bus` and `:cid_bus`
"""
function variable_mc_load_current(pm::_PMD.AbstractExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    var(pm, nw)[:crd] = Dict{Int, Any}()
    var(pm, nw)[:cid] = Dict{Int, Any}()
    var(pm, nw)[:crd_bus] = Dict{Int, Any}()
    var(pm, nw)[:cid_bus] = Dict{Int, Any}()
end


"""
	function variable_mc_load_power(
		pm::_PMD.AbstractNLExplicitNeutralIVRModel;
		nw=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

For non-linear IVR models with explicit neutrals,
creates placeholder dictionaries for the load power `:pd` and `:qd`,
and for the terminal power flows `:pd_bus` and `:qd_bus`
"""
function variable_mc_load_power(pm::_PMD.AbstractNLExplicitNeutralIVRModel; nw=nw_id_default, bounded::Bool=true, report::Bool=true)
    var(pm, nw)[:pd] = Dict{Int, Any}()
    var(pm, nw)[:qd] = Dict{Int, Any}()
    var(pm, nw)[:pd_bus] = Dict{Int, Any}()
    var(pm, nw)[:qd_bus] = Dict{Int, Any}()
end




# GENERATOR

# GENERATOR - Variables

"""
	function variable_mc_generator_current_real(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates generator real current variables `:crg` for models with explicit neutrals
"""
function variable_mc_generator_current_real(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _PMD._infer_int_dim_unit(gen, false) for (i,gen) in _PMD.ref(pm, nw, :gen))
    crg = var(pm, nw)[:crg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_crg_$(i)",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :gen, i), "crg_start", c, 0.0)
        ) for i in _PMD.ids(pm, nw, :gen)
    )

    report && _IM.sol_component_value(pm, _PMD.pmd_it_sym, nw, :gen, :crg, _PMD.ids(pm, nw, :gen), crg)
end


"""
	function variable_mc_generator_current_imaginary(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates generator imaginary current variables `:cig` for models with explicit neutrals
"""
function variable_mc_generator_current_imaginary(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _PMD._infer_int_dim_unit(gen, false) for (i,gen) in _PMD.ref(pm, nw, :gen))
    cig = var(pm, nw)[:cig] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_cig_$(i)",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :gen, i), "cig_start", c, 0.0)
        ) for i in _PMD.ids(pm, nw, :gen)
    )

    report && _IM.sol_component_value(pm, _PMD.pmd_it_sym, nw, :gen, :cig, _PMD.ids(pm, nw, :gen), cig)
end


"""
	function variable_mc_generator_power_real(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates generator active power variables `:pg` for models with explicit neutrals
"""
function variable_mc_generator_power_real(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _PMD._infer_int_dim_unit(gen, false) for (i,gen) in _PMD.ref(pm, nw, :gen))
    pg = var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_pg_$(i)",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :gen, i), "pg_start", c, 0.0)
        ) for i in _PMD.ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in _PMD.ref(pm, nw, :gen)
            set_lower_bound.(pg[i], gen["pmin"])
            set_upper_bound.(pg[i], gen["pmax"])
        end
    end

    report && _IM.sol_component_value(pm, _PMD.pmd_it_sym, nw, :gen, :pg, _PMD.ids(pm, nw, :gen), pg)
end


"""
	function variable_mc_generator_power_imaginary(
		pm::_PMD.ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates generator reactive power variables `:qg` for models with explicit neutrals
"""
function variable_mc_generator_power_imaginary(pm::_PMD.ExplicitNeutralModels; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => _PMD._infer_int_dim_unit(gen, false) for (i,gen) in _PMD.ref(pm, nw, :gen))
    qg = var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_qg_$(i)",
            start = _PMD.comp_start_value(_PMD.ref(pm, nw, :gen, i), "qg_start", c, 0.0)
        ) for i in _PMD.ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in _PMD.ref(pm, nw, :gen)
            set_lower_bound.(qg[i], gen["qmin"])
            set_upper_bound.(qg[i], gen["qmax"])
        end
    end

    report && _IM.sol_component_value(pm, _PMD.pmd_it_sym, nw, :gen, :qg, _PMD.ids(pm, nw, :gen), qg)
end


""
function variable_mc_generator_current(pm::_PMD.AbstractUnbalancedIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_generator_current_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_generator_current_imaginary(pm; nw=nw, bounded=bounded, report=report)

    var(pm, nw)[:crg_bus] = Dict{Int, Any}()
    var(pm, nw)[:cig_bus] = Dict{Int, Any}()

    # store active and reactive power expressions for use in objective + post processing
    var(pm, nw)[:pg] = Dict{Int, Any}()
    var(pm, nw)[:qg] = Dict{Int, Any}()
end


# GENERATOR - Variables - Non-linear

"""
	function variable_mc_generator_power(
		pm::_PMD.AbstractNLExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true
	)

For IVR models with explicit neutrals,
no power variables are required
"""
function variable_mc_generator_power(pm::_PMD.AbstractNLExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # do nothing
end


# GENERATOR - Variables - Quadratic

"""
	function variable_mc_generator_power(
		pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true
	)

For quadratic IVR models with explicit neutrals,
creates generator power variables `:pg` and `:qg`
"""
function variable_mc_generator_power(pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_generator_power_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_generator_power_imaginary(pm; nw=nw, bounded=bounded, report=report)
end


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


# TRANSFORMER - Variables

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
            start = comp_start_value(_PMD.ref(pm, nw, :transformer, l), "cit_start", c, 0.0)
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


# TRANSFORMER - Variable - Quadratic

"""
	function variable_mc_transformer_power(
		pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
	)

For quadratic IVR models with explicit neutrals,
creates transformer power variables `:pt` and `:qt`
"""
function variable_mc_transformer_power(pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_transformer_power_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_transformer_power_imaginary(pm; nw=nw, bounded=bounded, report=report)
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






