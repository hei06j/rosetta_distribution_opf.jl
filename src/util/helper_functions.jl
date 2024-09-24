function get_ref_bus_branch(ref)
    ref_bus = [i for (i,bus) in ref[:ref_buses]][1]
    ref_gen = ref[:bus_gens][ref_bus][1]
    ref_arc = [(branchid, fbus, tbus) for (branchid, fbus, tbus) in ref[:arcs_branch] if fbus == ref_bus][1]
    ref_branch = ref_arc[1]
    return ref_gen, ref_bus, ref_arc, ref_branch
end

function get_pv_bus_branch(ref)
    pv_buses = [busid for (busid,genid) in ref[:bus_gens] if (!isempty(genid) && busid ∉ keys(ref[:ref_buses]))]
    pv_genids = [genid[1] for (busid,genid) in ref[:bus_gens] if (!isempty(genid) && busid ∉ keys(ref[:ref_buses]))]
    pv_arcs = [ref[:bus_arcs_branch][i][1] for i in pv_buses]
    pv_branches = first.(pv_arcs)
    return pv_genids, pv_buses, pv_arcs, pv_branches
end
