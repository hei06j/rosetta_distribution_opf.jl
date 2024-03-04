### objectives: cost, loss, VUF, VUF2, PVUR, LVUR, IUF, IUF2, PIUR, PPUR (x), PQUR (x)
if objective == "cost"
    obj_cost = JuMP.@expression(model, sum(gen["cost"][1]*sum(pg[:,i]) + gen["cost"][2] for (i,gen) in ref[:gen]))
    # JuMP.@constraint(model, obj_cost >= 0)
    JuMP.@objective(model, Min, obj_cost)

elseif objective == "loss"
    branch_id = [i for (i, branch) in ref[:branch] if occursin("internal", branch["name"])][1]
    br_r = ref[:branch][branch_id]["br_r"]
    br_x = ref[:branch][branch_id]["br_x"]
    JuMP.@objective(model, Min, sum(sqrt.(br_r.^2 .+ br_x.^2) * Array(csr[:,branch_id].^2 .+ csi[:,branch_id].^2) ))

elseif objective in ["VUF" "VUF2" "PVUR" "LVUR"]
    gen_buses = [i for (i, gen) in ref[:gen] if i!==1]
    # terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref[:bus])
    terminals = Dict(i => collect(1:3) for i in gen_buses)
    vm = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm_$i", lower_bound = 0 ) for i in gen_buses)
    vm = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm[i].axes[1] ? vm[i][t,i] : 0.0 for t in 1:n_ph, i in gen_buses]), 1:n_ph, gen_buses)
    vr_012 = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vr_012") for i in gen_buses)
    vr_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vr_012[i].axes[1] ? vr_012[i][t,i] : 0.0 for t in 1:n_ph, i in gen_buses]), 1:n_ph, gen_buses)
    vi_012 = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vi_012") for i in gen_buses)
    vi_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vi_012[i].axes[1] ? vi_012[i][t,i] : 0.0 for t in 1:n_ph, i in gen_buses]), 1:n_ph, gen_buses)
    vm_012 = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm_012") for i in gen_buses)
    vm_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm_012[i].axes[1] ? vm_012[i][t,i] : 0.0 for t in 1:n_ph, i in gen_buses]), 1:n_ph, gen_buses)

    # vm = JuMP.@variable(model, [t in 1:3,], base_name="vm", lower_bound = 0 )
    # vr_012 = JuMP.@variable(model, [t in 1:3], base_name="vr_012")
    # vi_012 = JuMP.@variable(model, [t in 1:3], base_name="vi_012")
    # vm_012 = JuMP.@variable(model, [t in 1:3], base_name="vm_012", lower_bound = 0)

    for i in gen_buses
        # i = 1
        bus = ref[:bus][i]
        # JuMP.@constraint(model, vm[:,i].^2 .== vr[:,i].^2 .+ vi[:,i].^2)
        terminals = ref[:bus][i]["terminals"][1:end-1]
        JuMP.@constraint(model, vm[terminals,i].^2 .== vr[terminals,i].^2 .+ vi[terminals,i].^2)
        JuMP.@constraint(model, vr_012[terminals,i] .== Tre * Array(vr[terminals,i]) .- Tim * Array(vi[terminals,i]))
        JuMP.@constraint(model, vi_012[terminals,i] .== Tre * Array(vi[terminals,i]) .+ Tim * Array(vr[terminals,i]))
        JuMP.@constraint(model, vm_012[terminals,i].^2 .== vr_012[terminals,i].^2 .+ vi_012[terminals,i].^2)
    end

    if objective == "VUF"
        vuf = JuMP.@expression(model, Array(vm_012[3,:] ./ vm_012[2,:]))
        JuMP.@objective(model, Min, vuf)

    elseif objective == "VUF2"
        JuMP.@objective(model, Min, sum(vm_012[3,:]))
        # vuf2 = JuMP.@expression(model, [vm_012[3,2], vm_012[3,3], vm_012[3,4], vm_012[3,5], vm_012[3,6], vm_012[3,7]])
        # JuMP.@objective(model, Min, vuf2)

    elseif objective == "PVUR"
        phase_voltage = JuMP.@variable(model, [i=1], base_name="phase_voltage_$i")
        # phase_voltage = JuMP.@variable(model, [i in keys(ref[:bus])], base_name="phase_voltage_$i")
        # for (i, bus) in ref[:bus]
            i = 1
            bus = ref[:bus][i]
            terminals = ref[:bus][i]["terminals"][1:end-1]
            JuMP.@constraint(model, [t in terminals], phase_voltage[i] >= vm[t,i] - sum(vm[terminals,i])/3)
        # end
        JuMP.@objective(model, Min, phase_voltage[1] / (sum(vm[terminals,1])/3) )

    elseif objective == "LVUR"
        # line_voltage = JuMP.@variable(model, [i in keys(ref[:bus])], base_name="line_voltage_$i")
        # terminals = Dict(i => bus["terminals"][(!).(bus["grounded"])] for (i, bus) in ref[:bus])
        # vm_ll = Dict(i => JuMP.@variable(model, [t in terminals[i], i], base_name="vm_ll_$i") for i in keys(ref[:bus]))
        # vm_ll = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([t in vm_ll[i].axes[1] ? vm_ll[i][t,i] : 0.0 for t in 1:n_ph, i in keys(ref[:bus])]), 1:n_ph, keys(ref[:bus]))
        line_voltage = JuMP.@variable(model, base_name="line_voltage")
        vm_ll = JuMP.@variable(model, [t in 1:3], base_name="vm_ll", lower_bound=0)
        vr_ll = JuMP.@variable(model, [t in 1:3], base_name="vr_ll")
        vi_ll = JuMP.@variable(model, [t in 1:3], base_name="vi_ll")
        # for (i, bus) in ref[:bus]
            # terminals = ref[:bus][i]["terminals"][1:end-1]
            # JuMP.@constraint(model, vm_ll[terminals,i] .== Array(vm[terminals,i]) .- Array(vm[terminals2,i]))
            # JuMP.@constraint(model, [t in terminals], line_voltage[i] >= vm_ll[t,i] - sum(vm_ll[terminals,i])/3)
            i = 1
            bus = ref[:bus][i]
            terminals = collect(1:3)
            terminals2 = [terminals[2:end]..., terminals[1]]
            # JuMP.@constraint(model, vm_ll[terminals] .== vm[terminals] .- vm[terminals2])
            JuMP.@constraint(model, vr_ll[terminals] .== Array(vr[terminals,i]) .- Array(vr[terminals2,i]))
            JuMP.@constraint(model, vi_ll[terminals] .== Array(vi[terminals,i]) .- Array(vi[terminals2,i]))
            JuMP.@constraint(model, vm_ll[terminals].^2 .== vr_ll[terminals].^2 .+ vi_ll[terminals].^2)

            JuMP.@constraint(model, line_voltage .>= vm_ll[terminals] .- sum(vm_ll[terminals])/3)
        # end
        JuMP.@objective(model, Min, line_voltage / (sum(vm_ll[terminals,1])/3) )
    end


elseif objective in ["IUF" "IUF2" "PIUR"]
    int_dim = Dict(i => RPMD._infer_int_dim_unit(gen, !(4 in gen["connections"])) for (i,gen) in ref[:gen])
    cmg = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cmg_$i", lower_bound=0) for i in keys(ref[:gen]))
    cmg = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? cmg[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
    crg_012 = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="crg_012_$i") for i in keys(ref[:gen]))
    crg_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? crg_012[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
    cig_012 = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cig_012_$i") for i in keys(ref[:gen]))
    cig_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? cig_012[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))
    cmg_012 = Dict(i => JuMP.@variable(model, [c in 1:int_dim[i]], base_name="cmg_012_$i", lower_bound=0) for i in keys(ref[:gen]))
    cmg_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in 1:int_dim[i] ? cmg_012[i][c] : 0.0 for c in 1:n_ph, i in keys(ref[:gen])]), 1:n_ph, keys(ref[:gen]))

    for (i, gen) in ref[:gen]
        if 4 in gen["connections"]
            phases = ref[:gen][i]["connections"][1:end-1]
        else
            phases = ref[:gen][i]["connections"]
        end
        JuMP.@constraint(model, cmg[phases,i].^2 .== crg_bus[phases,i].^2 .+ crg_bus[phases,i].^2)
        JuMP.@constraint(model, crg_012[phases,i] .== Tre * Array(crg_bus[phases,i]) .- Tim * Array(cig_bus[phases,i]))
        JuMP.@constraint(model, cig_012[phases,i] .== Tre * Array(cig_bus[phases,i]) .+ Tim * Array(crg_bus[phases,i]))
        JuMP.@constraint(model, cmg_012[phases,i].^2 .== crg_012[phases,i].^2 .+ cig_012[phases,i].^2)
    end

    if objective == "IUF"
        JuMP.@objective(model, Min, cmg_012[3,1] / cmg_012[2,1])

    elseif objective == "IUF2"
        JuMP.@objective(model, Min, cmg_012[3,1])
        
    elseif objective == "PIUR"
        phase_current = JuMP.@variable(model, [i in keys(ref[:gen])], base_name="phase_current_$i")
        # for (i, gen) in ref[:gen]
            i = 1
            gen = ref[:gen][i]
            connections = ref[:gen][i]["connections"][1:end-1]
            JuMP.@constraint(model, [t in connections], phase_current[i] >= cmg[t,i] - sum(cmg[connections,i])/3)
        # end
        JuMP.@objective(model, Min, phase_current[1] / (sum(cmg[connections,1])/3) )
    end


elseif objective in ["IUF_inv" "IUF2_inv"]
    branch_id = 1
    branch = ref[:branch][branch_id]
    arc =  (1, 2, 1)
    nconds = Dict(l => length(branch["f_connections"]) for (l,branch) in ref[:branch])
    conds = Dict(l => branch["f_connections"] for (l,branch) in ref[:branch])
    
    cm = Dict((l,i,j) => JuMP.@variable(model, [c in conds[l]], base_name="cm_$((l,i,j))", lower_bound=0) for (l,i,j) in [arc])
    cm = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in cm[(l,i,j)].axes[1] ? cm[(l,i,j)][c] : 0.0 for c in 1:n_ph, (l,i,j) in [arc]]), 1:n_ph, [arc])
    cr_012 = Dict((l,i,j) => JuMP.@variable(model, [c in conds[l]], base_name="cr_012_$((l,i,j))") for (l,i,j) in [arc])
    cr_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in cr_012[(l,i,j)].axes[1] ? cr_012[(l,i,j)][c] : 0.0 for c in 1:n_ph, (l,i,j) in [arc]]), 1:n_ph, [arc])
    ci_012 = Dict((l,i,j) => JuMP.@variable(model, [c in conds[l]], base_name="ci_012_$((l,i,j))") for (l,i,j) in [arc])
    ci_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in ci_012[(l,i,j)].axes[1] ? ci_012[(l,i,j)][c] : 0.0 for c in 1:n_ph, (l,i,j) in [arc]]), 1:n_ph, [arc])
    cm_012 = Dict((l,i,j) => JuMP.@variable(model, [c in conds[l]], base_name="cm_012_$((l,i,j))", lower_bound=0) for (l,i,j) in [arc])
    cm_012 = JuMP.Containers.DenseAxisArray(Matrix{JuMP.AffExpr}([c in cm_012[(l,i,j)].axes[1] ? cm_012[(l,i,j)][c] : 0.0 for c in 1:n_ph, (l,i,j) in [arc]]), 1:n_ph, [arc])

    phases = 1:3
    JuMP.@constraint(model, cm[phases,arc].^2 .== cr_bus[phases,arc].^2 .+ ci_bus[phases,arc].^2)
    JuMP.@constraint(model, cr_012[phases,arc] .== Tre * Array(cr_bus[phases,arc]) .- Tim * Array(ci_bus[phases,arc]))
    JuMP.@constraint(model, ci_012[phases,arc] .== Tre * Array(ci_bus[phases,arc]) .+ Tim * Array(cr_bus[phases,arc]))
    JuMP.@constraint(model, cm_012[phases,arc].^2 .== cr_012[phases,arc].^2 .+ ci_012[phases,arc].^2)

    if objective == "IUF"
        # JuMP.@objective(model, Min, cmg_012[3,1] / cmg_012[2,1])
        JuMP.@objective(model, Min, cm_012[3,arc] / cm_012[2,arc])

    elseif objective == "IUF2"
        # JuMP.@objective(model, Min, cmg_012[3,1])
        JuMP.@objective(model, Min, cm_012[3,arc])
    
    end



elseif objective == "PPUR"
    phase_p = JuMP.@variable(model, base_name="phase_p", lower_bound=0)
    # phase_p = JuMP.@variable(model, [i in keys(ref[:gen])], base_name="phase_p_$i")
    # for (i, gen) in ref[:gen]
        i = 1
        gen = ref[:gen][i]
        connections = gen["connections"][1:end-1]
        JuMP.@constraint(model, phase_p .>= pg[connections,i] .- sum(pg[connections,i])/3)
        # JuMP.@constraint(model, [t in connections], phase_p[i] >= pg[t,i] - sum(pg[connections,i])/3)
    # end
    # JuMP.@objective(model, Min, phase_p[1] / (sum(pg[connections,1])/3) )
    JuMP.@objective(model, Min, phase_p / (sum(pg[connections,1])/3) )


elseif objective == "PQUR"
    phase_q = JuMP.@variable(model, [i=1], base_name="phase_q_$i")
    # phase_q = JuMP.@variable(model, [i in keys(ref[:gen])], base_name="phase_q_$i")
    # for (i, gen) in ref[:gen]
        i = 1
        gen = ref[:gen][i]
        connections = ref[:gen][i]["connections"][1:end-1]
        JuMP.@constraint(model, [c in connections], phase_q[i] >= qg[connections,i] - sum(qg[connections,i])/3)
    # end
    JuMP.@objective(model, Min, phase_q[1] / (sum(qg[connections,1])/3) )

end