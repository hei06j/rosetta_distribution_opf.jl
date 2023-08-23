"""
	function build_mc_opf(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for OPF in current-voltage variable space with explicit neutrals
"""
function build_mc_opf(pm::_PMD.AbstractExplicitNeutralIVRModel)
    # Variables
    variable_mc_bus_voltage(pm)
    variable_mc_branch_current(pm)
    variable_mc_load_current(pm)
    variable_mc_load_power(pm)
    variable_mc_generator_current(pm)
    variable_mc_generator_power(pm)
    variable_mc_transformer_current(pm)
    variable_mc_transformer_power(pm)
    variable_mc_switch_current(pm)

    # Constraints
    for i in _PMD.ids(pm, :bus)

        if i in _PMD.ids(pm, :ref_buses)
            constraint_mc_voltage_reference(pm, i)
        end

        constraint_mc_voltage_absolute(pm, i)
        constraint_mc_voltage_pairwise(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in _PMD.ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
        constraint_mc_generator_current(pm, id)
    end

    for id in _PMD.ids(pm, :load)
        constraint_mc_load_power(pm, id)
        constraint_mc_load_current(pm, id)
    end

    for i in _PMD.ids(pm, :transformer)
        constraint_mc_transformer_voltage(pm, i)
        constraint_mc_transformer_current(pm, i)

        constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in _PMD.ids(pm, :branch)
        constraint_mc_current_from(pm, i)
        constraint_mc_current_to(pm, i)
        constraint_mc_bus_voltage_drop(pm, i)

        constraint_mc_branch_current_limit(pm, i)
        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    # for i in _PMD.ids(pm, :switch)
    #     constraint_mc_switch_current(pm, i)
    #     constraint_mc_switch_state(pm, i)

    #     constraint_mc_switch_current_limit(pm, i)
    #     constraint_mc_switch_thermal_limit(pm, i)
    # end

    for i in _PMD.ids(pm, :bus)
        constraint_mc_current_balance(pm, i)
    end

    # Objective
    objective_mc_min_fuel_cost_polynomial_linquad(pm)
end