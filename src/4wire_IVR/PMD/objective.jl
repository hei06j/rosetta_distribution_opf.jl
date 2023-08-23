"gen connections adaptation of min fuel cost polynomial linquad objective"
function objective_mc_min_fuel_cost_polynomial_linquad(pm)
    pg_contains_nl_exp = any(x<:JuMP.NonlinearExpression for x in vcat([typeof.(isa(pg, JuMP.Containers.DenseAxisArray) ? pg.data : pg) for nw in _PMD.nw_ids(pm) for (id,pg) in _PMD.var(pm, nw, :pg)]...))
    gen_cost = Dict()

    if !pg_contains_nl_exp
        for (n, nw_ref) in _PMD.nws(pm)
            for (i,gen) in nw_ref[:gen]
                pg = sum(_PMD.var(pm, n, :pg, i))

                if length(gen["cost"]) == 1
                    gen_cost[(n,i)] = gen["cost"][1]
                elseif length(gen["cost"]) == 2
                    gen_cost[(n,i)] = gen["cost"][1]*pg + gen["cost"][2]
                elseif length(gen["cost"]) == 3
                    gen_cost[(n,i)] = gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3]
                else
                    gen_cost[(n,i)] = 0.0
                end
            end
        end

        return JuMP.@objective(pm.model, Min,
            sum(
                sum( gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] )
            for (n, nw_ref) in _PMD.nws(pm))
        )
    else
        for (n, nw_ref) in _PMD.nws(pm)
            for (i,gen) in nw_ref[:gen]
                bus = gen["gen_bus"]

                #to avoid function calls inside of @NLconstraint:
                pg = _PMD.var(pm, n, :pg, i)
                pg = isa(pg, JuMP.Containers.DenseAxisArray) ? pg.data : pg

                int_dim = length(pg)
                if length(gen["cost"]) == 1
                    gen_cost[(n,i)] = gen["cost"][1]
                elseif length(gen["cost"]) == 2
                    gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*sum(pg[i] for i in 1:int_dim) + gen["cost"][2])
                elseif length(gen["cost"]) == 3
                    gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*sum(pg[i] for i in 1:int_dim)^2 + gen["cost"][2]*sum(pg[i] for i in 1:int_dim) + gen["cost"][3])
                else
                    gen_cost[(n,i)] = 0.0
                end
            end
        end

        return JuMP.@NLobjective(pm.model, Min,
            sum(
                sum(    gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] )
            for (n, nw_ref) in _PMD.nws(pm))
        )
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


