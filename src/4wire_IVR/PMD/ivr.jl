"""
	function build_mc_opf(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for OPF in current-voltage variable space with explicit neutrals
"""
function build_mc_opf(pm::PMD.AbstractExplicitNeutralIVRModel)
    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_load_current(pm)
    PMD.variable_mc_load_power(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_generator_power(pm)
    variable_mc_transformer_current(pm)
    variable_mc_transformer_power(pm)
    variable_mc_switch_current(pm)

    # Constraints
    for i in PMD.ids(pm, :bus)

        if i in PMD.ids(pm, :ref_buses)
            PMD.constraint_mc_voltage_reference(pm, i)
        end

        PMD.constraint_mc_voltage_absolute(pm, i)
        PMD.constraint_mc_voltage_pairwise(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
        PMD.constraint_mc_generator_current(pm, id)
    end

    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
        PMD.constraint_mc_load_current(pm, id)
    end

    for i in PMD.ids(pm, :transformer)
        constraint_mc_transformer_voltage(pm, i)
        constraint_mc_transformer_current(pm, i)

        constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_current_from(pm, i)
        PMD.constraint_mc_current_to(pm, i)
        PMD.constraint_mc_bus_voltage_drop(pm, i)

        PMD.constraint_mc_branch_current_limit(pm, i)
        PMD.constraint_mc_thermal_limit_from(pm, i)
        PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    # for i in PMD.ids(pm, :switch)
    #     constraint_mc_switch_current(pm, i)
    #     constraint_mc_switch_state(pm, i)

    #     constraint_mc_switch_current_limit(pm, i)
    #     constraint_mc_switch_thermal_limit(pm, i)
    # end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
    end

    # Objective
    PMD.objective_mc_min_fuel_cost_polynomial_linquad(pm)
end