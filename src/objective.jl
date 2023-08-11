
"""
    objective_mc_min_fuel_cost(pm::AbstractUnbalancedPowerModel)

Standard fuel cost minimization objective
"""
function objective_mc_min_fuel_cost(pm::_PMD.AbstractUnbalancedPowerModel; report::Bool=true)
    model = _PMD.check_gen_cost_models(pm)

    if model == 1
        return _PMD.objective_mc_min_fuel_cost_pwl(pm; report=report)
    elseif model == 2
        return _PMD.objective_mc_min_fuel_cost_polynomial(pm; report=report)
    else
        error("Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end



# """
#     objective_mc_min_fuel_cost_pwl(pm::AbstractUnbalancedPowerModel)

# Fuel cost minimization objective with piecewise linear terms
# """
# function objective_mc_min_fuel_cost_pwl(pm::_PMD.AbstractUnbalancedPowerModel; report::Bool=true)
#     objective_mc_variable_pg_cost(pm; report=report)

#     return JuMP.@objective(pm.model, Min,
#         sum(
#             sum( var(pm, n, :pg_cost, i) for (i,gen) in nw_ref[:gen])
#         for (n, nw_ref) in nws(pm))
#     )
# end


# """
#     objective_mc_min_fuel_cost_polynomial(pm::AbstractUnbalancedPowerModel)

# Fuel cost minimization objective for polynomial terms
# """
# function objective_mc_min_fuel_cost_polynomial(pm::AbstractUnbalancedPowerModel; report::Bool=true)
#     order = calc_max_cost_index(pm.data)-1

#     if order <= 2
#         return _objective_mc_min_fuel_cost_polynomial_linquad(pm; report=report)
#     else
#         return _objective_mc_min_fuel_cost_polynomial_nl(pm; report=report)
#     end
# end


