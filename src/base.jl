# Helper functions for AbstractPowerModel `var` access.
var(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int=nw_id_default) = _IM.var(pm, _PMD.pmd_it_sym, nw)
var(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, key::Symbol) = _IM.var(pm, _PMD.pmd_it_sym, nw, key)
var(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, key::Symbol, idx::Any) = _IM.var(pm, _PMD.pmd_it_sym, nw, key, idx)
var(pm::_PMD.AbstractUnbalancedPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.var(pm, _PMD.pmd_it_sym, key; nw = nw)
var(pm::_PMD.AbstractUnbalancedPowerModel, key::Symbol, idx::Any; nw::Int = nw_id_default) = _IM.var(pm, _PMD.pmd_it_sym, key, idx; nw = nw)

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