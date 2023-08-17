
# LOAD - Constraints - Quadratic

# """
# 	function constraint_mc_load_power(
# 		pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel,
# 		id::Int;
# 		nw::Int=nw_id_default,
# 		report::Bool=true
# 	)

# For quadratic IVR models with explicit neutrals,
# link the load power variables `:pd` and `:qd` to the voltage,
# and link together the power, voltage and current variables
# """
# function constraint_mc_load_power(pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true)
#     load = ref(pm, nw, :load, id)
#     bus = ref(pm, nw,:bus, load["load_bus"])

#     configuration = load["configuration"]
#     int_dim = _PMD._infer_int_dim_unit(load, false)
#     a, alpha, b, beta = _PMD._load_expmodel_params(load, bus)

#     # Note that one-dimensional delta loads are handled as wye-connected loads.
#     # The distinction between one-dimensional wye and delta loads is purely semantic
#     # when neutrals are modeled explicitly.
#     if configuration==_PMD.WYE || int_dim==1
#         constraint_mc_load_power_wye(pm, nw, id, load["load_bus"], load["connections"], load["model"], a, b; report=report)
#     else
#         constraint_mc_load_power_delta(pm, nw, id, load["load_bus"], load["connections"], load["model"], a, b; report=report)
#     end
# end


# """
# 	function constraint_mc_load_power_wye(
# 		pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel,
# 		nw::Int,
# 		id::Int,
# 		bus_id::Int,
# 		connections::Vector{Int},
# 		model::LoadModel,
# 		a::Vector{<:Real},
# 		b::Vector{<:Real};
# 		report::Bool=true,
# 		bounded::Bool=true
# 	)

# For quadratic IVR models with explicit neutrals,
# link the load power variables `:pd` and `:qd` to the voltage,
# and link together the power, voltage and current variables
# for wye-connected loads
# """
# function constraint_mc_load_power_wye(pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, model::_PMD.LoadModel, a::Vector{<:Real}, b::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
#     vr = var(pm, nw, :vr, bus_id)
#     vi = var(pm, nw, :vi, bus_id)

#     crd = var(pm, nw, :crd, id)
#     cid = var(pm, nw, :cid, id)

#     phases = connections[1:end-1]
#     n      = connections[end]

#     vr_pn = [vr[p]-vr[n] for p in phases]
#     vi_pn = [vi[p]-vi[n] for p in phases]

#     if model==POWER
#         pd = a
#         qd = b
#     elseif model==IMPEDANCE
#         pd = a .* (vr_pn.^2 .+ vi_pn.^2)
#         qd = b .* (vr_pn.^2 .+ vi_pn.^2)
#     elseif model==CURRENT
#         pd = var(pm, nw, :pd, id)
#         qd = var(pm, nw, :qd, id)
#         JuMP.@constraint(pm.model, pd.^2 .== a.^2 .* (vr_pn.^2 .+ vi_pn.^2))
#         JuMP.@constraint(pm.model, sign.(a).*pd .>= 0)
#         JuMP.@constraint(pm.model, qd.^2 .== b.^2 .* (vr_pn.^2 .+ vi_pn.^2))
#         JuMP.@constraint(pm.model, sign.(b).*qd .>= 0)
#     else
#         error("Load model $model for load $id is not supported by this formulation.")
#     end

#     JuMP.@constraint(pm.model, pd .==  vr_pn.*crd .+ vi_pn.*cid)
#     JuMP.@constraint(pm.model, qd .== -vr_pn.*cid .+ vi_pn.*crd)

#     # constant current loads are already reported through variable function
#     if report && model!=CURRENT
#         _PMD.sol(pm, nw, :load, id)[:pd] = pd
#         _PMD.sol(pm, nw, :load, id)[:qd] = qd
#     end
# end


# """
# 	function constraint_mc_load_power_delta(
# 		pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel,
# 		nw::Int,
# 		id::Int,
# 		bus_id::Int,
# 		connections::Vector{Int},
# 		model::LoadModel,
# 		a::Vector{<:Real},
# 		b::Vector{<:Real};
# 		report::Bool=true,
# 		bounded::Bool=true
# 	)

# For quadratic IVR models with explicit neutrals,
# link the load power variables `:pd` and `:qd` to the voltage,
# and link together the power, voltage and current variables
# for delta-connected loads
# """
# function constraint_mc_load_power_delta(pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, model::_PMD.LoadModel, a::Vector{<:Real}, b::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
#     vr = var(pm, nw, :vr, bus_id)
#     vi = var(pm, nw, :vi, bus_id)

#     crd = var(pm, nw, :crd, id)
#     cid = var(pm, nw, :cid, id)

#     phases = connections

#     Md = _get_delta_transformation_matrix(length(connections))
#     vrd = Md*[vr[p] for p in phases]
#     vid = Md*[vi[p] for p in phases]

#     if model==POWER
#         pd = a
#         qd = b
#     elseif model==IMPEDANCE
#         pd = a .* (vrd.^2 .+ vid.^2)
#         qd = b .* (vrd.^2 .+ vid.^2)
#     elseif model==CURRENT
#         pd = var(pm, nw, :pd, id)
#         qd = var(pm, nw, :qd, id)
#         JuMP.@constraint(pm.model, pd.^2 .== a.^2 .* (vrd.^2 .+ vid.^2))
#         JuMP.@constraint(pm.model, sign.(a).*pd .>= 0)
#         JuMP.@constraint(pm.model, qd.^2 .== b.^2 .* (vrd.^2 .+ vid.^2))
#         JuMP.@constraint(pm.model, sign.(b).*qd .>= 0)
#     else
#         error("Load model $model for load $id is not supported by this formulation.")
#     end

#     JuMP.@constraint(pm.model, pd .==  vrd.*crd .+ vid.*cid)
#     JuMP.@constraint(pm.model, qd .== -vrd.*cid .+ vid.*crd)

#     # constant current loads are already reported through variable function
#     if report && model!=CURRENT
#         _PMD.sol(pm, nw, :load, id)[:pd] = pd
#         _PMD.sol(pm, nw, :load, id)[:qd] = qd
#     end
# end


# """
# 	function constraint_mc_load_current(
# 		pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel,
# 		id::Int;
# 		nw::Int=nw_id_default,
# 		report::Bool=true,
# 		bounded::Bool=true
# 	)

# For quadratic IVR models with explicit neutrals,
# create expressions for the terminal current flows `:crd_bus` and `:cid_bus`
# """
# function constraint_mc_load_current(pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel, id::Int; nw::Int=nw_id_default, report::Bool=true, bounded::Bool=true)
#     load = ref(pm, nw, :load, id)

#     int_dim = _PMD._infer_int_dim_unit(load, false)
#     # Note that one-dimensional delta loads are handled as wye-connected loads.
#     # The distinction between one-dimensional wye and delta loads is purely semantic
#     # when neutrals are modeled explicitly.
#     if get(load, "configuration", _PMD.WYE) == _PMD.WYE || int_dim==1
#         constraint_mc_load_current_wye(pm, nw, id, load["connections"]; report=report, bounded=bounded)
#     else
#         constraint_mc_load_current_delta(pm, nw, id, load["connections"]; report=report, bounded=bounded)
#     end
# end


# """
# 	function constraint_mc_load_current_wye(
# 		pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel,
# 		nw::Int,
# 		id::Int,
# 		connections::Vector{Int};
# 		report::Bool=true,
# 		bounded::Bool=true
# 	)

# For quadratic IVR models with explicit neutrals,
# create expressions for the terminal current flows `:crd_bus` and `:cid_bus`
# for wye-connected loads
# """
# function constraint_mc_load_current_wye(pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
#     crd = var(pm, nw, :crd, id)
#     cid = var(pm, nw, :cid, id)
#     var(pm, nw, :crd_bus)[id] = _PMD._merge_bus_flows(pm, [crd..., -sum(crd)], connections)
#     var(pm, nw, :cid_bus)[id] = _PMD._merge_bus_flows(pm, [cid..., -sum(cid)], connections)
# end


# """
# 	function constraint_mc_load_current_delta(
# 		pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel,
# 		nw::Int,
# 		id::Int,
# 		connections::Vector{Int};
# 		report::Bool=true,
# 		bounded::Bool=true
# 	)

# For quadratic IVR models with explicit neutrals,
# create expressions for the terminal current flows `:crd_bus` and `:cid_bus`
# for delta-connected loads
# """
# function constraint_mc_load_current_delta(pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, connections::Vector{Int}; report::Bool=true, bounded::Bool=true)
#     crd = var(pm, nw, :crd, id)
#     cid = var(pm, nw, :cid, id)
#     Md = _get_delta_transformation_matrix(length(connections))
#     var(pm, nw, :crd_bus)[id] = _PMD._merge_bus_flows(pm, Md'*crd, connections)
#     var(pm, nw, :cid_bus)[id] = _PMD._merge_bus_flows(pm, Md'*cid, connections)
# end



# # GENERATOR - Constraints - Quadratic

# """
# 	function constraint_mc_generator_power_wye(
# 		pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel,
# 		nw::Int,
# 		id::Int,
# 		bus_id::Int,
# 		connections::Vector{Int},
# 		pmin::Vector{<:Real},
# 		pmax::Vector{<:Real},
# 		qmin::Vector{<:Real},
# 		qmax::Vector{<:Real};
# 		report::Bool=true,
# 		bounded::Bool=true
# 	)

# For quadratic IVR models with explicit neutrals,
# links the generator power variables `:pd` and `:qd`
# of wye-connected generators to the voltage and current
# """
# function constraint_mc_generator_power_wye(pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
#     vr = var(pm, nw, :vr, bus_id)
#     vi = var(pm, nw, :vi, bus_id)
#     crg = var(pm, nw, :crg, id)
#     cig = var(pm, nw, :cig, id)

#     phases = connections[1:end-1]
#     n      = connections[end]

#     vr_pn = [vr[p]-vr[n] for p in phases]
#     vi_pn = [vi[p]-vi[n] for p in phases]

#     pg = var(pm, nw, :pg, id)
#     qg = var(pm, nw, :qg, id)

#     JuMP.@constraint(pm.model, pg .==  vr_pn.*crg .+ vi_pn.*cig)
#     JuMP.@constraint(pm.model, qg .== -vr_pn.*cig .+ vi_pn.*crg)
# end


# """
# 	function constraint_mc_generator_power_delta(
# 		pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel,
# 		nw::Int,
# 		id::Int,
# 		bus_id::Int,
# 		connections::Vector{Int},
# 		pmin::Vector{<:Real},
# 		pmax::Vector{<:Real},
# 		qmin::Vector{<:Real},
# 		qmax::Vector{<:Real};
# 		report::Bool=true,
# 		bounded::Bool=true
# 	)

# For quadratic IVR models with explicit neutrals,
# links the generator power variables `:pd` and `:qd`
# of delta-connected generators to the voltage and current
# """
# function constraint_mc_generator_power_delta(pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}; report::Bool=true, bounded::Bool=true)
#     vr = var(pm, nw, :vr, bus_id)
#     vi = var(pm, nw, :vi, bus_id)
#     crg = var(pm, nw, :crg, id)
#     cig = var(pm, nw, :cig, id)

#     ph = connections
#     ph_next = [connections[2:end]..., connections[1]]

#     vrg = [vr[c]-vr[d] for (c,d) in zip(ph,ph_next)]
#     vig = [vi[c]-vi[d] for (c,d) in zip(ph,ph_next)]

#     pg = var(pm, nw, :pg, id)
#     qg = var(pm, nw, :qg, id)

#     JuMP.@constraint(pm.model, pg .==  vrg.*crg .+ vig.*cig)
#     JuMP.@constraint(pm.model, qg .== -vrg.*cig .+ vig.*crg)
# end



# # TRANSFORMER - Constraint - Quadratic

# """
# 	function constraint_mc_transformer_thermal_limit(
# 		pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel,
# 		nw::Int,
# 		id::Int,
# 		f_idx::Tuple,
# 		t_idx::Tuple,
# 		f_bus::Int,
# 		t_bus::Int,
# 		f_connections::Vector,
# 		t_connections::Vector,
# 		config::ConnConfig,
# 		sm_ub::Real;
# 		report::Bool=true
# 	)

# For quadratic IVR models with explicit neutrals,
# imposes a bound on the magnitude of the total apparent power at both windings.

# ```
# sum(pt_fr)^2 + sum(qt_fr)^2 <= sm_ub^2
# sum(pt_to)^2 + sum(qt_to)^2 <= sm_ub^2
# ```
# """
# function constraint_mc_transformer_thermal_limit(pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel, nw::Int, id::Int, f_idx::Tuple, t_idx::Tuple, f_bus::Int, t_bus::Int, f_connections::Vector, t_connections::Vector, config, sm_ub::Real; report::Bool=true)
#     vr_fr = var(pm, nw, :vr, f_bus)
#     vi_fr = var(pm, nw, :vi, f_bus)
#     vr_to = var(pm, nw, :vr, t_bus)
#     vi_to = var(pm, nw, :vi, t_bus)

#     crt_fr = var(pm, nw, :crt, f_idx)
#     cit_fr = var(pm, nw, :cit, f_idx)
#     crt_to = var(pm, nw, :crt, t_idx)
#     cit_to = var(pm, nw, :cit, t_idx)

#     if config==_PMD.WYE || length(crt_fr)==1
#         P_fr = f_connections[1:end-1]
#         n_fr = f_connections[end]
#         vrt_fr = [vr_fr[p]-vr_fr[n_fr] for p in P_fr]
#         vit_fr = [vi_fr[p]-vi_fr[n_fr] for p in P_fr]
#     elseif config==_PMD.DELTA && length(crt_fr)==3
#         M = _get_delta_transformation_matrix(3)
#         vrt_fr = M*[vr_to[p] for p in f_connections]
#         vit_fr = M*[vi_to[p] for p in f_connections]
#     else
#         error("The configuration $config of dimension $(length(crt)) is not supported.")
#     end

#     P_to = t_connections[1:end-1]
#     n_to = t_connections[end]
#     vrt_to = [vr_to[p]-vr_to[n_to] for p in P_to]
#     vit_to = [vi_to[p]-vi_to[n_to] for p in P_to]

#     pt_fr = var(pm, nw, :pt, f_idx)
#     qt_fr = var(pm, nw, :qt, f_idx)
#     pt_to = var(pm, nw, :pt, t_idx)
#     qt_to = var(pm, nw, :qt, t_idx)

#     JuMP.@constraint(pm.model, pt_fr .==  vrt_fr.*crt_fr .+ vit_fr.*cit_fr)
#     JuMP.@constraint(pm.model, qt_fr .== -vrt_fr.*cit_fr .+ vit_fr.*crt_fr)
#     JuMP.@constraint(pm.model, pt_to .==  vrt_to.*crt_to .+ vit_to.*cit_to)
#     JuMP.@constraint(pm.model, qt_to .== -vrt_to.*cit_to .+ vit_to.*crt_to)

#     if sm_ub < Inf
#         JuMP.@constraint(pm.model, sum(pt_fr)^2+sum(qt_fr)^2 <= sm_ub^2)
#         JuMP.@constraint(pm.model, sum(pt_to)^2+sum(qt_to)^2 <= sm_ub^2)
#     end
# end


# """
#     function constraint_mc_switch_power(
#         pm::_PMD.ExplicitNeutralModels,
#         id::Int;
#         nw::Int=nw_id_default,
#         report::Bool=true
#     )

# For IVR models with explicit neutrals,
# link the switch power or create appropiate expressions for them
# """
# function constraint_mc_switch_power(pm::_PMD.ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, report::Bool=true)
#     switch = ref(pm, nw, :switch, id)
#     f_bus = switch["f_bus"]
#     t_bus = switch["t_bus"]
#     f_idx = (id, f_bus, t_bus)
#     t_idx = (id, t_bus, f_bus)

#     constraint_mc_switch_power(pm, nw, id, f_idx, t_idx, switch["f_connections"], switch["t_connections"])
# end


# """
# 	function constraint_mc_switch_thermal_limit(
# 		pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel,
# 		nw::Int,
# 		f_idx::Tuple{Int,Int,Int},
# 		f_connections::Vector{Int},
# 		rating::Vector{<:Real}
# 	)

# For quadratic IVR models with explicit neutrals,
# throw an error because this cannot be represented quadratically
# without introducing explicit power variables.
# """
# function constraint_mc_switch_thermal_limit(pm::_PMD.AbstractQuadraticExplicitNeutralIVRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, rating::Vector{<:Real})
#     @warn("""
#         A switch power bound cannot be represented quadratically in the default AbstractQuadraticExplicitNeutralIVRModel.
#         Either extend this quadratic formulation by including explicit switch power variables, or use AbstractNLExplicitNeutralIVRModel instead.
#         """)
# end

