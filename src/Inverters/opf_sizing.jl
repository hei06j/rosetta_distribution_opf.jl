"""
	function variable_mc_converter_pdclink(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates generator active power variables `:pg` for models with explicit neutrals
"""
function variable_mc_converter_pdclink(pm::PMD.ExplicitNeutralModels; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    converter_ids = [1]

    pdclink_sqr = PMD.var(pm, nw)[:pdclink_sqr] = Dict(i => JuMP.@variable(pm.model,
            base_name="$(nw)_pdclink_sqr",
        ) for i in converter_ids
    )

    if bounded
        for (i,gen) in PMD.ref(pm, nw, :gen)
            if i in converter_ids
                PMD.set_lower_bound.(pdclink_sqr[i], gen["pdcmin"])
                PMD.set_upper_bound.(pdclink_sqr[i], gen["pdcmax"])
            end
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :gen, :pdclink_sqr, converter_ids, pdclink_sqr)
end



"""
    function constraint_mc_inverter_dc_link_ripple_power(
        pm::ExplicitNeutralModels,
        id::Int;
        nw::Int=PMD.nw_id_default,
        report::Bool=true
    )

Constrains generator power variables for models with explicit neutrals.
"""
function constraint_mc_converter_pdclink(pm::PMD.ExplicitNeutralModels, id::Int; nw::Int=PMD.nw_id_default, report::Bool=true)
    generator = PMD.ref(pm, nw, :gen, id)
    bus = PMD.ref(pm, nw,:bus, generator["gen_bus"])

    pdcmin = get(generator, "pdcmin", 0) #-Inf)
    pdcmax = get(generator, "pdcmax", Inf)
    
    constraint_mc_converter_pdclink(pm, nw, id, bus["index"], generator["connections"], pdcmin, pdcmax; report=report)
end


"""
	function constraint_mc_inverter_dc_link_ripple_power(
		pm::AbstractNLExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
        pdcmin::Vector{<:Real},
		pdcmax::Vector{<:Real},
		report::Bool=true
	)

For IVR models with explicit neutrals,
creates non-linear expressions for the inverter dc link power `:pdc_link_sqr`
of wye-connected generators as a function of voltage and current
"""
function constraint_mc_converter_pdclink(pm::PMD.AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pdcmin::Real, pdcmax::Real; report::Bool=true)
    # bus_id = 1
    vr = PMD.var(pm, nw, :vr, bus_id)
    vi = PMD.var(pm, nw, :vi, bus_id)
    crg = PMD.var(pm, nw, :crg, id)
    cig = PMD.var(pm, nw, :cig, id)
    crg_bus = PMD.var(pm, nw, :crg_bus)[id]
    cig_bus = PMD.var(pm, nw, :cig_bus)[id]
    pdclink_sqr = PMD.var(pm, nw, :pdclink_sqr, id)
    
    @show pdclink_sqr
    @show connections
    @show [vr[p]*crg_bus[idx] for (idx, p) in enumerate(connections)]

    JuMP.@constraint(pm.model,  pdclink_sqr ==
        sum( vr[p]*crg_bus[idx] .- vi[p]*cig_bus[idx] for (idx, p) in enumerate(connections) )^2
        + 
        sum( vr[p]*cig_bus[idx] .+ vi[p]*crg_bus[idx] for (idx, p) in enumerate(connections) )^2
    )
end



"""
function build_mc_opf_mx(
    pm::AbstractExplicitNeutralIVRModel
)

constructor for OPF in current-voltage variable space with explicit neutrals including reconfigurable inverters
"""
function build_mc_opf_sizing(pm::PMD.AbstractExplicitNeutralIVRModel)
    converter_branches = [branch["index"] for (i, branch) in pm.data["branch"] if occursin("inverter_branch", branch["name"])]
    converter_ids, converter_bus_ids, converter_branch_ids = get_pv_bus_branch(pm.ref[:it][:pmd][:nw][0])  # maybe output branch_ids and arcs seperately?
    branch_converters = Dict(branch_id[1] => converter_ids[i] for (i, branch_id) in enumerate(converter_branch_ids))


    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_load_current(pm)
    PMD.variable_mc_load_power(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_generator_power(pm)
    # variable_mc_generator_current(pm)
    # variable_mc_generator_power(pm)
    PMD.variable_mc_transformer_current(pm)
    PMD.variable_mc_transformer_power(pm)
    PMD.variable_mc_switch_current(pm)

    if pm.setting["dc_link"]
        variable_mc_converter_pdclink(pm)
        # PMD.var(pm, 0)[:pdc_link_sqr] = Dict{Int, Any}()
    end

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
        if id ∈ converter_ids  # Generators connected with inverter
            constraint_mc_generator_power(pm, id)
            PMD.constraint_mc_generator_current(pm, id)
            # constraint_mc_generator_current(pm, id)
            constraint_mc_generator_current_limit(pm, id)
            
            if pm.setting["dc_link"]
                # constraint_mc_inverter_dc_link_ripple_power(pm, id)
                constraint_mc_converter_pdclink(pm, id)
            end

        else  # Other generators
            PMD.constraint_mc_generator_power(pm, id)
            PMD.constraint_mc_generator_current(pm, id)
        end
    end

    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
        PMD.constraint_mc_load_current(pm, id)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_voltage(pm, i)
        PMD.constraint_mc_transformer_current(pm, i)

        PMD.constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :branch)

        if i ∈ converter_branches
            PMD.constraint_mc_current_from(pm, i)
            PMD.constraint_mc_current_to(pm, i)
            PMD.constraint_mc_bus_voltage_drop(pm, i)
            PMD.constraint_mc_branch_current_limit(pm, i)

        else  # normal branch
            PMD.constraint_mc_current_from(pm, i)
            PMD.constraint_mc_current_to(pm, i)
            PMD.constraint_mc_bus_voltage_drop(pm, i)
            PMD.constraint_mc_branch_current_limit(pm, i)
            # PMD.constraint_mc_thermal_limit_from(pm, i)
            # PMD.constraint_mc_thermal_limit_to(pm, i)
        end

    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_current(pm, i)
        PMD.constraint_mc_switch_state(pm, i)

        PMD.constraint_mc_switch_current_limit(pm, i)
        PMD.constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
    end

    # Objective
    # PMD.objective_mc_min_fuel_cost(pm)
    # objective_mc_min_IUF(pm)
    objective_mc_min_max_phase_current(pm)
    # objective_mc_min_ref_branch_loss(pm)
end
