
"""
function variable_mc_generator_power_rating(
    pm::ExplicitNeutralModels;
    nw::Int=nw_id_default,
    bounded::Bool=true,
    report::Bool=true
)

Creates generator active power variables `:pg` for models with explicit neutrals
"""
function variable_mc_generator_power_rating(pm::PMD.ExplicitNeutralModels; nw::Int=PMD.nw_id_default, report::Bool=true)
    converter_ids = [1]  # TODO fix this
    
    # int_dim = Dict(i => PMD._infer_int_dim_unit(gen, false) for (i,gen) in PMD.ref(pm, nw, :gen))
    srating = Dict(i => JuMP.@variable(pm.model,
            base_name="srating_$(i)",
        ) for i in converter_ids
    )

    for (n, nw) in PMD.nws(pm)
        PMD.var(pm, n)[:srating] = srating
    end

    for (i,gen) in PMD.ref(pm, nw, :gen)
        if i in converter_ids
            PMD.set_lower_bound.(PMD.var(pm, nw)[:srating][i], gen["srating_min"])
            PMD.set_upper_bound.(PMD.var(pm, nw)[:srating][i], gen["srating_max"])
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :gen, :srating, converter_ids, srating)
end


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
    converter_ids = [1] # TODO fix this

    pdclink_sqr = PMD.var(pm, nw)[:pdclink_sqr] = Dict(i => JuMP.@variable(pm.model,
            base_name="$(nw)_pdclink_sqr",
        ) for i in converter_ids
    )

    if bounded
        for (i,gen) in PMD.ref(pm, nw, :gen)
            if i in converter_ids
                # PMD.set_lower_bound.(pdclink_sqr[i], gen["pdcmin"]^2)
                # PMD.set_upper_bound.(pdclink_sqr[i], gen["pdcmax"]^2)
                PMD.set_lower_bound.(pdclink_sqr[i], 0)
                PMD.set_upper_bound.(pdclink_sqr[i], gen["pdcrating_max"]^2)
            end
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :gen, :pdclink_sqr, converter_ids, pdclink_sqr)
end



"""
	function variable_mc_converter_pdclink(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		report::Bool=true
	)

Creates generator active power variables `:pg` for models with explicit neutrals
"""
function variable_mc_converter_pdcrating(pm::PMD.ExplicitNeutralModels; nw::Int=PMD.nw_id_default, report::Bool=true)
    converter_ids = [1] # TODO fix this

    pdcrating = Dict(i => JuMP.@variable(pm.model,
            base_name="pdcrating_$(i)",
        ) for i in converter_ids
    )

    for (n, nw) in PMD.nws(pm)
        PMD.var(pm, n)[:pdcrating] = pdcrating
    end

    for (i,gen) in PMD.ref(pm, nw, :gen)
        if i in converter_ids
            PMD.set_lower_bound.(PMD.var(pm, nw)[:pdcrating][i], gen["pdcrating_min"])
            PMD.set_upper_bound.(PMD.var(pm, nw)[:pdcrating][i], gen["pdcrating_max"])
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :gen, :pdcrating, converter_ids, pdcrating)
end


function constraint_mc_converter_pdclink_rating(pm::PMD.ExplicitNeutralModels, id::Int; nw::Int=PMD.nw_id_default)
    constraint_mc_converter_pdclink_rating(pm, nw, id)
end

function constraint_mc_converter_pdclink_rating(pm::PMD.AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int)
    pdclink_sqr = PMD.var(pm, nw, :pdclink_sqr, id)
    pdcrating = PMD.var(pm, nw, :pdcrating, id)

    JuMP.@constraint(pm.model, pdclink_sqr <= pdcrating^2)
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


# GENERATOR - Constraints

"""
    function constraint_mc_generator_power_rating(
        pm::ExplicitNeutralModels,
        id::Int;
        nw::Int=PMD.nw_id_default,
        report::Bool=true
    )

Constrains generator power variables for models with explicit neutrals.
"""
function constraint_mc_generator_power_rating(pm::PMD.ExplicitNeutralModels, id::Int; nw::Int=PMD.nw_id_default, report::Bool=true)
    generator = PMD.ref(pm, nw, :gen, id)
    bus = PMD.ref(pm, nw,:bus, generator["gen_bus"])

    configuration = generator["configuration"]

    N = length(generator["connections"])

    if configuration==PMD.WYE || length(pmin)==1
        constraint_mc_generator_power_rating_wye(pm, nw, id, bus["index"], generator["connections"]; report=report)
    else
        # constraint_mc_generator_power_rating_delta(pm, nw, id, bus["index"], generator["connections"], pmin, pmax, qmin, qmax; report=report)
    end
end


function constraint_mc_generator_current_limit_rating(pm::PMD.AbstractExplicitNeutralIVRModel, id::Int; nw::Int=PMD.nw_id_default, report::Bool=true, bounded::Bool=true)
    generator = PMD.ref(pm, nw, :gen, id)
    bus = PMD.ref(pm, nw,:bus, generator["gen_bus"])

    # if haskey(generator, "c_rating") && any(generator["c_rating"] .< Inf)
    constraint_mc_generator_current_limit_rating(pm, nw, id, bus["index"], generator["connections"]; report=report, bounded=bounded)
    # end
end


"""
	function constraint_mc_generator_power_rating_wye(
		pm::AbstractNLExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		bus_id::Int,
		connections::Vector{Int},
		report::Bool=true
	)

For IVR models with explicit neutrals,
creates non-linear expressions for the generator power `:pd` and `:qd`
of wye-connected generators as a function of voltage and current
"""
function constraint_mc_generator_power_rating_wye(pm::PMD.AbstractNLExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}; report::Bool=true)
    vr = PMD.var(pm, nw, :vr, bus_id)
    vi = PMD.var(pm, nw, :vi, bus_id)
    crg = PMD.var(pm, nw, :crg, id)
    cig = PMD.var(pm, nw, :cig, id)

    phases = connections[1:end-1]
    n      = connections[end]

    pg = JuMP.NonlinearExpr[]
    qg = JuMP.NonlinearExpr[]

    srating = PMD.var(pm, nw, :srating, id)
    pmin = zero(3)              # TODO when considering batteries, this should be -srating/3 * ones(3)
    pmax = srating/3 * ones(3)
    qmin = -srating/3 * ones(3)
    qmax = srating/3 * ones(3)

    for (idx, p) in enumerate(phases)
        push!(pg, JuMP.@expression(pm.model,   (vr[p]-vr[n])*crg[idx]  + (vi[p]-vi[n])*cig[idx]))
        push!(qg, JuMP.@expression(pm.model, - (vr[p]-vr[n])*cig[idx]  + (vi[p]-vi[n])*crg[idx]))
    end

    JuMP.@constraint(pm.model, sum(pmin) <= sum(pg))
    JuMP.@constraint(pm.model, sum(pmax) >= sum(pg))
    JuMP.@constraint(pm.model, sum(qmin) <= sum(qg))
    JuMP.@constraint(pm.model, sum(qmax) >= sum(qg))

    PMD.var(pm, nw, :pg)[id] = pg
    PMD.var(pm, nw, :qg)[id] = qg

    if report
        PMD.sol(pm, nw, :gen, id)[:pg] = pg
        PMD.sol(pm, nw, :gen, id)[:qg] = qg
    end
end


"""
	function constraint_mc_generator_current_limit_rating(
		pm::AbstractExplicitNeutralIVRModel,
		nw::Int,
		id::Int,
		connections::Vector{Int};
		report::Bool=true,
		bounded::Bool=true
	)

For IVR models with explicit neutrals,
creates expressions for the terminal current flows `:crg_bus` and `:cig_bus` of wye-connected generators
"""
function constraint_mc_generator_current_limit_rating(pm::PMD.AbstractExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
    crg_bus = PMD.var(pm, nw, :crg_bus)[id]
    cig_bus = PMD.var(pm, nw, :cig_bus)[id]
    
    vr = PMD.var(pm, nw, :vr, bus_id)
    vi = PMD.var(pm, nw, :vi, bus_id)

    phases = connections[1:end-1]
    n      = connections[end]

    srating = PMD.var(pm, nw, :srating, id)

    JuMP.@constraint(pm.model, ((vr[phases].-vr[n]).^2 .+ (vi[phases].-vi[n]).^2) .* (crg_bus[phases].^2 .+ cig_bus[phases].^2) .<= (srating*ones(3)/3).^2)

    # @assert length(c_rating) == length(crg_bus)
    # cnds_finite_nonzero_rating = [c for (c,r) in enumerate(c_rating) if (r<Inf && r!==0)]
    # cnds_zero_rating = [c for (c,r) in enumerate(c_rating) if r==0]

    # JuMP.@constraint(pm.model, [c in cnds_finite_nonzero_rating], crg_bus[c]^2+cig_bus[c]^2 <= c_rating[c]^2)
    # JuMP.@constraint(pm.model, [c in cnds_zero_rating], crg_bus[c] == 0)
    # JuMP.@constraint(pm.model, [c in cnds_zero_rating], cig_bus[c] == 0)
end


function objective_mc_min_sizing(pm::PMD.AbstractUnbalancedPowerModel)
    id = 1
    n = 1
    obj = PMD.var(pm, n, :srating, id) + PMD.var(pm, n, :pdcrating, id)
    # obj = sum(
    #        PMD.var(pm, n, :srating, id) + PMD.var(pm, n, :pdcrating, id)
    #     for (n, nw_ref) in PMD.nws(pm))

    return JuMP.@objective(pm.model, Min, obj)
end

"""
function build_mc_opf_sizing(
    pm::AbstractExplicitNeutralIVRModel
)

constructor for OPF sizing in current-voltage variable space with explicit neutrals including reconfigurable inverters
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
    variable_mc_generator_power_rating(pm)

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

            constraint_mc_generator_power_rating(pm, id)
            PMD.constraint_mc_generator_current(pm, id)
            # constraint_mc_generator_current(pm, id)
            constraint_mc_generator_current_limit_rating(pm, id)
            
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





"""
function build_mn_mc_opf_sizing(
    pm::AbstractExplicitNeutralIVRModel
)

constructor for multi-network OPF sizing in current-voltage variable space with explicit neutrals including reconfigurable inverters
"""
function build_mn_mc_opf_sizing(pm::PMD.AbstractExplicitNeutralIVRModel)
    converter_branches = [branch["index"] for (i, branch) in pm.data["nw"]["1"]["branch"] if occursin("inverter_branch", branch["name"])]
    converter_ids, converter_bus_ids, converter_branch_ids = get_pv_bus_branch(pm.ref[:it][:pmd][:nw][1])  # maybe output branch_ids and arcs seperately?
    branch_converters = Dict(branch_id[1] => converter_ids[i] for (i, branch_id) in enumerate(converter_branch_ids))


    variable_mc_generator_power_rating(pm; nw=1)
    variable_mc_converter_pdcrating(pm; nw=1)

    # Variables
    for (n, network) in PMD.nws(pm)
        PMD.variable_mc_bus_voltage(pm; nw=n)
        PMD.variable_mc_branch_current(pm; nw=n)
        PMD.variable_mc_load_current(pm; nw=n)
        PMD.variable_mc_load_power(pm; nw=n)
        PMD.variable_mc_generator_current(pm; nw=n)
        PMD.variable_mc_generator_power(pm; nw=n)
        # variable_mc_generator_current(pm; nw=n)
        # variable_mc_generator_power(pm; nw=n)
        PMD.variable_mc_transformer_current(pm; nw=n)
        PMD.variable_mc_transformer_power(pm; nw=n)
        PMD.variable_mc_switch_current(pm; nw=n)

        if pm.setting["dc_link"]
            variable_mc_converter_pdclink(pm; nw=n)
            # PMD.var(pm, 0)[:pdc_link_sqr] = Dict{Int, Any}()
        end

        # Constraints
        for i in PMD.ids(pm, n, :bus)

            if i in PMD.ids(pm, n, :ref_buses)
                PMD.constraint_mc_voltage_reference(pm, i; nw=n)
            end

            PMD.constraint_mc_voltage_absolute(pm, i; nw=n)
            PMD.constraint_mc_voltage_pairwise(pm, i; nw=n)
        end

        # components should be constrained before KCL, or the bus current variables might be undefined

        for id in PMD.ids(pm, n, :gen)
            if id ∈ converter_ids  # Generators connected with inverter

                constraint_mc_generator_power_rating(pm, id; nw=n)
                PMD.constraint_mc_generator_current(pm, id; nw=n)
                # constraint_mc_generator_current(pm, id; nw=n)
                constraint_mc_generator_current_limit_rating(pm, id; nw=n)
                
                if pm.setting["dc_link"]
                    # constraint_mc_inverter_dc_link_ripple_power(pm, id; nw=n)
                    constraint_mc_converter_pdclink(pm, id; nw=n)
                    constraint_mc_converter_pdclink_rating(pm, id; nw=n)
                end

            else  # Other generators
                PMD.constraint_mc_generator_power(pm, id; nw=n)
                PMD.constraint_mc_generator_current(pm, id; nw=n)
            end
        end

        for id in PMD.ids(pm, n, :load)
            PMD.constraint_mc_load_power(pm, id; nw=n)
            PMD.constraint_mc_load_current(pm, id; nw=n)
        end

        for i in PMD.ids(pm, n, :transformer)
            PMD.constraint_mc_transformer_voltage(pm, i; nw=n)
            PMD.constraint_mc_transformer_current(pm, i; nw=n)

            PMD.constraint_mc_transformer_thermal_limit(pm, i; nw=n)
        end

        for i in PMD.ids(pm, n, :branch)

            if i ∈ converter_branches
                PMD.constraint_mc_current_from(pm, i; nw=n)
                PMD.constraint_mc_current_to(pm, i; nw=n)
                PMD.constraint_mc_bus_voltage_drop(pm, i; nw=n)
                PMD.constraint_mc_branch_current_limit(pm, i; nw=n)

            else  # normal branch
                PMD.constraint_mc_current_from(pm, i; nw=n)
                PMD.constraint_mc_current_to(pm, i; nw=n)
                PMD.constraint_mc_bus_voltage_drop(pm, i; nw=n)
                PMD.constraint_mc_branch_current_limit(pm, i; nw=n)
                # PMD.constraint_mc_thermal_limit_from(pm, i; nw=n)
                # PMD.constraint_mc_thermal_limit_to(pm, i; nw=n)
            end

        end

        for i in PMD.ids(pm, n, :switch)
            PMD.constraint_mc_switch_current(pm, i; nw=n)
            PMD.constraint_mc_switch_state(pm, i; nw=n)

            PMD.constraint_mc_switch_current_limit(pm, i; nw=n)
            PMD.constraint_mc_switch_thermal_limit(pm, i; nw=n)
        end

        for i in PMD.ids(pm, n, :bus)
            PMD.constraint_mc_current_balance(pm, i; nw=n)
        end

    end

    # Objective
    objective_mc_min_sizing(pm)
    # PMD.objective_mc_min_fuel_cost(pm)
    # objective_mc_min_IUF(pm)
    # objective_mc_min_max_phase_current(pm)
    # objective_mc_min_ref_branch_loss(pm)
end
