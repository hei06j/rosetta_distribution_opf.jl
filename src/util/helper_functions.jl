function get_sequence_components(x)
    ### x must be complex
    alpha = exp(im*2/3*pi)
    T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
    Tre = real.(T)
    Tim = imag.(T)

    x_re = real.(x)
    x_im = imag.(x)
    
    x_seq_re = Tre * x_re .- Tim * x_im
    x_seq_im = Tre * x_im .+ Tim * x_re
    x_seq = x_seq_re + im*x_seq_im
    x_seq_m= abs.(x_seq)
    return x_seq_re, x_seq_im, x_seq_m
end

function sequence(x)
    @assert length(x) == 3

    alpha = exp(im*2/3*pi)
    T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
    Tre = real.(T)
    Tim = imag.(T)

    return T * x

end

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


function set_gen_max_powers!(data_math, gen_id, kva)
    data_math["gen"]["$gen_id"]["pmax"] = kva/3 * ones(3)
    data_math["gen"]["$gen_id"]["pmin"] = 0 * ones(3)
    data_math["gen"]["$gen_id"]["qmax"] = kva/3 * ones(3)
    data_math["gen"]["$gen_id"]["qmin"] = -kva/3 * ones(3)
end