# Helper functions for AbstractPowerModel `var` access.
var(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int=nw_id_default) = _IM.var(pm, _PMD.pmd_it_sym, nw)
var(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, key::Symbol) = _IM.var(pm, _PMD.pmd_it_sym, nw, key)
var(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, key::Symbol, idx::Any) = _IM.var(pm, _PMD.pmd_it_sym, nw, key, idx)
var(pm::_PMD.AbstractUnbalancedPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.var(pm, _PMD.pmd_it_sym, key; nw = nw)
var(pm::_PMD.AbstractUnbalancedPowerModel, key::Symbol, idx::Any; nw::Int = nw_id_default) = _IM.var(pm, _PMD.pmd_it_sym, key, idx; nw = nw)

# ref = _IM.build_ref(data_math, _PMD.ref_add_core!, _PMD._pmd_global_keys, _PMD.pmd_it_name)[:it][:pmd][:nw][0]

"""
	function set_lower_bound(
		x::JuMP.VariableRef,
		bound::Real
	)

Local wrapper method for JuMP.set_lower_bound, which skips NaN and infinite (-Inf only)
"""
function set_lower_bound(x::JuMP.VariableRef, bound::Real)
    if !(isnan(bound) || bound==-Inf)
        JuMP.set_lower_bound(x, bound)
    end
end


"""
	function set_lower_bound(
		xs::Vector{JuMP.VariableRef},
		bound::Real
	)

Local wrapper method for JuMP.set_lower_bound, which skips NaN and infinite (-Inf only).
Note that with this signature, the bound is applied to every variable in the vector.
"""
function set_lower_bound(xs::Vector{JuMP.VariableRef}, bound::Real)
    for x in xs
        set_lower_bound(x, bound)
    end
end


"""
	function set_upper_bound(
		x::JuMP.VariableRef,
		bound
	)

Local wrapper method for JuMP.set_upper_bound, which skips NaN and infinite (+Inf only)
"""
function set_upper_bound(x::JuMP.VariableRef, bound::Real)
    if !(isnan(bound) || bound==Inf)
        JuMP.set_upper_bound(x, bound)
    end
end


"""
	function set_upper_bound(
		xs::Vector{JuMP.VariableRef},
		bound::Real
	)

Local wrapper method for JuMP.set_upper_bound, which skips NaN and infinite (+Inf only).
Note that with this signature, the bound is applied to every variable in the vector.
"""
function set_upper_bound(xs::Vector{JuMP.VariableRef}, bound::Real)
    for x in xs
        set_upper_bound(x, bound)
    end
end


"""
    @smart_constraint model::JuMP.Model vars::Vector expr::JuMP.Expression

Detection of whether a constraint should be NL or not"
"""
macro smart_constraint(model, vars, expr)
    esc(quote
        if _PMD._has_nl_expression($vars)
            JuMP.@NLconstraint($model, $expr)
        else
            JuMP.@constraint($model, $expr)
        end
    end)
end



"infer the internal dimension of a winding, load or generator based on the connections and the configuration"
function _infer_int_dim(connections::Vector, configuration, kron_reduced)
    if configuration==_PMD.WYE
        if kron_reduced
            return length(connections)
        else
            return length(connections)-1
        end
    else # DELTA
        if length(connections)==2
            return 1
        elseif length(connections)==3
            return 3
        else
            error("Only 1 and 3 phase delta-connections are supported.")
        end
    end
end

"infer the internal dimension for a unit, i.e. any one-port component with `connections` and `configuration` properties"
function _infer_int_dim_unit(unit::Dict{String,<:Any}, kron_reduced)
    return _infer_int_dim(unit["connections"], unit["configuration"], kron_reduced)
end


"infer the internal dimension for a transformer (only in the MATHEMATICAL data model format)"
function _infer_int_dim_transformer(trans::Dict{String,<:Any}, kron_reduced)
    return _infer_int_dim(trans["f_connections"], trans["configuration"], kron_reduced)
end


"Merges flow variables that enter the same terminals, i.e. multiple neutrals of an underground cable connected to same neutral terminal"
function _merge_bus_flows(model, flows::Vector, connections::Vector)::JuMP.Containers.DenseAxisArray
    flows_merged = []
    conns_unique = unique(connections)
    for t in conns_unique
        idxs = findall(connections.==t)
        flows_t = flows[idxs]
        if length(flows_t)==1
            flows_merged_t = flows_t[1]
        elseif any(isa(a, JuMP.NonlinearExpression) for a in flows_t)
            flows_merged_t = JuMP.@NLexpression(model, sum(flows_t[i] for i in 1:length(flows_t)))
        else
            flows_merged_t = sum(flows_t)
        end
        push!(flows_merged, flows_merged_t)
    end
    JuMP.Containers.DenseAxisArray(flows_merged, conns_unique)
end


"checks if a sufficient number of variables exist for the given keys collection"
function _check_var_keys(vars, keys, var_name, comp_name)
    if length(vars) < length(keys)
        error("$(var_name) decision variables appear to be missing for $(comp_name) components")
    end
end


"helper function to build bus shunt matrices for power balance constraints"
function _build_bus_shunt_matrices(ref, terminals::Vector{Int}, bus_shunts::Vector{<:Tuple{Int,Vector{Int}}})::Tuple{Matrix{<:Real},Matrix{<:Real}}
    ncnds = length(terminals)
    Gs = fill(0.0, ncnds, ncnds)
    Bs = fill(0.0, ncnds, ncnds)
    for (i, connections) in bus_shunts
        shunt = ref[:shunt][i]
        for (idx,c) in enumerate(connections)
            for (jdx,d) in enumerate(connections)
                Gs[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] += shunt["gs"][idx,jdx]
                Bs[findfirst(isequal(c), terminals),findfirst(isequal(d), terminals)] += shunt["bs"][idx,jdx]
            end
        end
    end

    return (Gs, Bs)
end