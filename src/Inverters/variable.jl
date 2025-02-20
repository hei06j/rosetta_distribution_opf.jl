# Inverter - Variables
"given a variable that is indexed by component ids, builds the standard solution structure"
function sol_component_value_comp(aim::IM.AbstractInfrastructureModel, it::Symbol, n::Int, comp_name::Symbol, field_name::Symbol, comp_id, variable)
    @assert !haskey(IM.sol(aim, it, n, comp_name, comp_id), field_name)
    IM.sol(aim, it, n, comp_name, comp_id)[field_name] = variable
end

"given a variable that is indexed by component ids, builds the standard solution structure"
function sol_component_value(aim::IM.AbstractInfrastructureModel, it::Symbol, n::Int, comp_name::Symbol, field_name::Symbol, comp_ids, variables)
    for i in comp_ids
        @assert !haskey(IM.sol(aim, it, n, comp_name, i), field_name)
        IM.sol(aim, it, n, comp_name, i)[field_name] = variables[i]
    end
    return 1
end


function variable_reconfigurable_inverter(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    pv_gen_ids = [i for (i, gen) in pm.data["gen"] if !occursin("source", gen["name"])]
    # connections = Dict(i => pm.data["gen"]["$i"]["connections"] for i in pv_gen_ids)
    connections = Dict(i => length(pm.data["gen"]["$i"]["connections"]) for i in pv_gen_ids)
    m_legs = Dict(i => pm.data["gen"]["$i"]["m_legs"] for i in pv_gen_ids)

    bg = PMD.var(pm, nw)[:bg] = Dict(parse(Int, i) => JuMP.@variable(pm.model, [1:connections[i], 1:m_legs[i]], base_name="bg_$i", Bin) for i in pv_gen_ids)
    for i in pv_gen_ids
        report && sol_component_value_comp(pm, PMD.pmd_it_sym, nw, :gen, :bg, parse.(Int, i), bg[parse.(Int, i)])
    end
end


# GENERATOR - Variables
""
function variable_mc_generator_current(pm::PMD.AbstractUnbalancedIVRModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_generator_current_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_generator_current_imaginary(pm; nw=nw, bounded=bounded, report=report)

    PMD.var(pm, nw)[:crg_bus] = Dict{Int, Any}()
    PMD.var(pm, nw)[:cig_bus] = Dict{Int, Any}()

    # store active and reactive power expressions for use in objective + post processing
    PMD.var(pm, nw)[:pg] = Dict{Int, Any}()
    PMD.var(pm, nw)[:qg] = Dict{Int, Any}()
end


"""
	function variable_mc_generator_current_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates generator real current variables `:crg` for models with explicit neutrals
"""
function variable_mc_generator_current_real(pm::PMD.ExplicitNeutralModels; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict()
    for (i,gen) in PMD.ref(pm, nw, :gen)
        if occursin("pv", gen["name"])
            int_dim[i] = PMD._infer_int_dim_unit(gen, true)
        else
            int_dim[i] = PMD._infer_int_dim_unit(gen, false)
        end
    end
    # int_dim = Dict(i => PMD._infer_int_dim_unit(gen, true) for (i,gen) in PMD.ref(pm, nw, :gen))
    crg = PMD.var(pm, nw)[:crg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_crg_$(i)",
            start = PMD.comp_start_value(PMD.ref(pm, nw, :gen, i), "crg_start", c, 0.0)
        ) for i in PMD.ids(pm, nw, :gen)
    )

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :gen, :crg, PMD.ids(pm, nw, :gen), crg)
end


"""
	function variable_mc_generator_current_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates generator imaginary current variables `:cig` for models with explicit neutrals
"""
function variable_mc_generator_current_imaginary(pm::PMD.ExplicitNeutralModels; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict()
    for (i,gen) in PMD.ref(pm, nw, :gen)
        if occursin("pv", gen["name"])
            int_dim[i] = PMD._infer_int_dim_unit(gen, true)
        else
            int_dim[i] = PMD._infer_int_dim_unit(gen, false)
        end
    end
    # int_dim = Dict(i => PMD._infer_int_dim_unit(gen, true) for (i,gen) in PMD.ref(pm, nw, :gen))
    cig = PMD.var(pm, nw)[:cig] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_cig_$(i)",
            start = PMD.comp_start_value(PMD.ref(pm, nw, :gen, i), "cig_start", c, 0.0)
        ) for i in PMD.ids(pm, nw, :gen)
    )

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :gen, :cig, PMD.ids(pm, nw, :gen), cig)
end




"create variables for generators, delegate to PowerModels"
function variable_mc_generator_power(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_generator_power_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_generator_power_imaginary(pm; nw=nw, bounded=bounded, report=report)
end

"""
	function variable_mc_generator_power_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates generator active power variables `:pg` for models with explicit neutrals
"""
function variable_mc_generator_power_real(pm::PMD.ExplicitNeutralModels; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => PMD._infer_int_dim_unit(gen, false) for (i,gen) in PMD.ref(pm, nw, :gen))
    pg = PMD.var(pm, nw)[:pg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_pg_$(i)",
            start = PMD.comp_start_value(PMD.ref(pm, nw, :gen, i), "pg_start", c, 0.0)
        ) for i in PMD.ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in PMD.ref(pm, nw, :gen)
            PMD.set_lower_bound.(pg[i], gen["pmin"])
            PMD.set_upper_bound.(pg[i], gen["pmax"])
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :gen, :pg, PMD.ids(pm, nw, :gen), pg)
end


"""
	function variable_mc_generator_power_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates generator reactive power variables `:qg` for models with explicit neutrals
"""
function variable_mc_generator_power_imaginary(pm::PMD.ExplicitNeutralModels; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => PMD._infer_int_dim_unit(gen, false) for (i,gen) in PMD.ref(pm, nw, :gen))
    qg = PMD.var(pm, nw)[:qg] = Dict(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_qg_$(i)",
            start = PMD.comp_start_value(PMD.ref(pm, nw, :gen, i), "qg_start", c, 0.0)
        ) for i in PMD.ids(pm, nw, :gen)
    )

    if bounded
        for (i,gen) in PMD.ref(pm, nw, :gen)
            PMD.set_lower_bound.(qg[i], gen["qmin"])
            PMD.set_upper_bound.(qg[i], gen["qmax"])
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :gen, :qg, PMD.ids(pm, nw, :gen), qg)
end



"""
	function variable_mc_load_current_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates load real current variables `:cmd` for models with explicit neutrals
"""
function variable_mc_load_current_magnitude(pm::PMD.ExplicitNeutralModels; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    int_dim = Dict(i => PMD._infer_int_dim_unit(load, false) for (i,load) in PMD.ref(pm, nw, :load))

    cmd = PMD.var(pm, nw)[:cmd] = Dict{Int,Any}(i => JuMP.@variable(pm.model,
            [c in 1:int_dim[i]], base_name="$(nw)_cmd_$(i)",
            start = PMD.comp_start_value(ref(pm, nw, :load, i), "cmd_start", c, 0.0)
        ) for i in PMD.ids(pm, nw, :load)
    )

    report && IM.sol_component_value(pm, pmd_it_sym, nw, :load, :cmd, PMD.ids(pm, nw, :load), cmd)
end