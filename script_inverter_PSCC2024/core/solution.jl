function get_solutions(model, result)
    # JuMP.optimize!(model)
    # @assert(JuMP.termination_status(model) == LOCALLY_SOLVED)
    # cost = JuMP.objective_value(model)

    alpha = exp(im*2/3*pi)
    T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
    Tre = real.(T)
    Tim = imag.(T)
    
    pv_gen_id = [parse(Int,i) for (i,gen) in data_math["gen"] if occursin("pv", gen["name"])][1]
    gen_bus_id = data_math["gen"]["$pv_gen_id"]["gen_bus"]
    
    # vr_vals = [(i,bus["vr"]) for (i,bus) in result["solution"]["bus"]]
    # vi_vals = [(i,bus["vi"]) for (i,bus) in result["solution"]["bus"]]
    v_pv = result["solution"]["bus"]["$gen_bus_id"]["vr"] .+ im*result["solution"]["bus"]["$gen_bus_id"]["vi"]
    v_pv_012 = T * v_pv[1:3]
    vm_pv = abs.(v_pv)
    vm_pv_012 = abs.(v_pv_012)    # va1 = angle.(v_pv).*180/pi
    # va1_pn = va1[1:3] .- va1[4]
    # va1_pn_pp = Array(va1_pn[[1,2,3]]) .- Array(va1_pn[[2,3,1]])
    # va1_pp = Array(va1[[1,2,3]]) .- Array(va1[[2,3,1]])

    # crg_values = JuMP.value.(crg_bus)
    # cig_values = JuMP.value.(cig_bus)
    # crg_values[:,pv_gen_id] + im * cig_values[:,pv_gen_id]
    # cg_pv = [gen["crg"] + im * gen["cig"] for (i,gen) in result["solution"]["gen"] if parse(Int, i) ∈ pv_gen_id]
    cg_pv = result["solution"]["gen"]["$pv_gen_id"]["crg"] + im * result["solution"]["gen"]["$pv_gen_id"]["cig"]
    cg_pv_012 = T * cg_pv[1:3]
    cgm_pv = abs.(cg_pv)
    cgm_pv_012 = abs.(cg_pv_012)
    # cga1 = angle.(cg_pv).*180/pi

    cg_src = result["solution"]["gen"]["2"]["crg"] + im * result["solution"]["gen"]["2"]["cig"]
    cg_src_012 = T * cg_src[1:3]
    cgm_src = abs.(cg_src)
    cgm_src_012 = abs.(cg_src_012)

    cd = []#Array{Float64}(undef, 4, 0)
    cd_012 = []#Array{Float64}(undef, 3, 0)
    for load_id in keys(result["solution"]["load"])
        c = zeros(4) .+ im*zeros(4)
        connections = model.data["load"]["$load_id"]["connections"]
        c[connections] = result["solution"]["load"]["$load_id"]["crd_bus"] .+ im * result["solution"]["load"]["$load_id"]["cid_bus"]
        append!(cd, c)
        append!(cd_012, T*c[1:3])
    end

    # cd = [load["crd_bus"] .+ im * load["cid_bus"] for (i,load) in result["solution"]["load"]]
    # cd = sum(cd, dims=2)
    # cd_012 = T * cd[1:3]
    # cdm = abs.(cd)
    # cdm_012 = abs.(cd_012)
    
    # pv_branch = (1, 2, 1)
    _, _, _, pv_branch = RPMD.get_pv_bus_branch(model.ref[:it][:pmd][:nw][0])
    c_pv = []#Array{Float64}(undef, 4, 0)
    c_pv_012 = []#Array{Float64}(undef, 3, 0)
    for branch_id in pv_branch
        c = result["solution"]["branch"]["$branch_id"]["cr_fr"] .+ im * result["solution"]["branch"]["$branch_id"]["ci_fr"]
        append!(c_pv, c)
        append!(c_pv_012, T*c[1:3])
    end


    ref_gen, ref_bus, ref_arc, ref_branch = RPMD.get_ref_bus_branch(model.ref[:it][:pmd][:nw][0])
    c_ref = result["solution"]["branch"]["$ref_branch"]["cr_fr"] .+ im * result["solution"]["branch"]["$ref_branch"]["ci_fr"]
    c_ref_012 = T * c_ref[1:3]

    
    results = Dict()
    results["v_inv"] = round.(v_pv, digits=6)
    results["v_inv_012"] = round.(v_pv_012, digits=6)
    results["c_inv"] = round.(cg_pv, digits=6)
    results["c_inv_012"] = round.(cg_pv_012, digits=6)
    results["c_source"] = round.(cg_src, digits=6)
    results["c_source_012"] = round.(cg_src_012, digits=6)
    results["c_load"] = round.(cd, digits=6)
    results["c_load_012"] = round.(cd_012, digits=6)
    results["c_inverter_branch"] = round.(c_pv, digits=6)
    results["c_inverter_branch_012"] = round.(c_pv_012, digits=6)
    results["c_source_branch"] = round.(c_ref, digits=6)
    results["c_source_branch_012"] = round.(c_ref_012, digits=6)
    # results["sg_inv"] = Array(value.(pg)[:,1] .+ im * value.(qg)[:,1])
    # results["sg_source"] = Array(value.(pg)[:,2] .+ im * value.(qg)[:,2])
    
    return results
end
