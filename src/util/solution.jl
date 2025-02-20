function get_solutions(model, result)
    pv_gen_id = [i for (i,gen) in PMD.ref(model, 0, :gen) if occursin("pv", gen["name"])][1]
    gen_bus_id = PMD.ref(model, 0, :gen, pv_gen_id)["gen_bus"]
    
    v_pv = result["solution"]["bus"]["$gen_bus_id"]["vr"] .+ im*result["solution"]["bus"]["$gen_bus_id"]["vi"]
    v_pv_012 = sequence(v_pv[1:3])
    vm_pv = abs.(v_pv)
    vm_pv_012 = abs.(v_pv_012)    # va1 = angle.(v_pv).*180/pi
    # va1_pn = va1[1:3] .- va1[4]
    # va1_pn_pp = Array(va1_pn[[1,2,3]]) .- Array(va1_pn[[2,3,1]])
    # va1_pp = Array(va1[[1,2,3]]) .- Array(va1[[2,3,1]])

    cg_pv = result["solution"]["gen"]["$pv_gen_id"]["crg"] + im * result["solution"]["gen"]["$pv_gen_id"]["cig"]
    cg_pv_012 = sequence(cg_pv[1:3])
    cgm_pv = abs.(cg_pv)
    cgm_pv_012 = abs.(cg_pv_012)
    # cga1 = angle.(cg_pv).*180/pi

    cg_src = result["solution"]["gen"]["2"]["crg"] + im * result["solution"]["gen"]["2"]["cig"]
    cg_src_012 = sequence(cg_src[1:3])
    cgm_src = abs.(cg_src)
    cgm_src_012 = abs.(cg_src_012)

    cd = []#Array{Float64}(undef, 4, 0)
    cd_012 = []#Array{Float64}(undef, 3, 0)
    for load_id in keys(result["solution"]["load"])
        c = zeros(4) .+ im*zeros(4)
        connections = model.data["load"]["$load_id"]["connections"]
        c[connections] = result["solution"]["load"]["$load_id"]["crd_bus"] .+ im * result["solution"]["load"]["$load_id"]["cid_bus"]
        append!(cd, c)
        append!(cd_012, sequence(c[1:3]))
    end
    
    _, _, _, pv_branch = get_pv_bus_branch(PMD.ref(model, 0))
    c_pv = []#Array{Float64}(undef, 4, 0)
    c_pv_012 = []#Array{Float64}(undef, 3, 0)
    for branch_id in pv_branch
        c = result["solution"]["branch"]["$branch_id"]["cr_fr"] .+ im * result["solution"]["branch"]["$branch_id"]["ci_fr"]
        append!(c_pv, c)
        append!(c_pv_012, sequence(c[1:3]))
    end

    ref_gen, ref_bus, ref_arc, ref_branch = get_ref_bus_branch(PMD.ref(model, 0))
    c_ref = result["solution"]["branch"]["$ref_branch"]["cr_fr"] .+ im * result["solution"]["branch"]["$ref_branch"]["ci_fr"]
    c_ref_012 = sequence(c_ref[1:3])
    
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
