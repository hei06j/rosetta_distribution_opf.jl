function add_inverter_losses!(data_math, gen_id; c_rating_a=0, reconfigurable=false, GFM=false, three_wire=false, dc_link=true)
    gen = data_math["gen"]["$gen_id"]
    old_gen_bus = copy(gen["gen_bus"])
    new_gen_bus = length(data_math["bus"]) + 1
    
    if three_wire
        gen["connections"] = gen["connections"][1:3]
        data_math["bus"]["$old_gen_bus"]["terminals"] = data_math["bus"]["$old_gen_bus"]["terminals"][1:4]
        data_math["bus"]["$old_gen_bus"]["grounded"] = data_math["bus"]["$old_gen_bus"]["grounded"][1:4]
    end

    data_math["bus"]["$new_gen_bus"] = deepcopy(data_math["bus"]["$old_gen_bus"])
    data_math["bus"]["$new_gen_bus"]["bus_i"] = new_gen_bus
    data_math["bus"]["$new_gen_bus"]["index"] = new_gen_bus
    data_math["bus"]["$new_gen_bus"]["name"] = "GFL_bus"
    data_math["bus"]["$new_gen_bus"]["type"] = "GFL"
    data_math["bus"]["$new_gen_bus"]["bus_type"] = 2
    data_math["bus"]["$old_gen_bus"]["bus_type"] = 1
    gen["gen_bus"] = new_gen_bus

    if GFM
        data_math["bus"]["$new_gen_bus"]["vm"] = [1. ; 1 ; 1 ; 0]
        data_math["bus"]["$new_gen_bus"]["va"] = [0. ; -2.094 ; 2.094 ; 0.]
        data_math["bus"]["$new_gen_bus"]["grounded"] = Bool[0, 0, 0, 1]
        data_math["bus"]["$new_gen_bus"]["name"] = "GFM_bus"
        data_math["bus"]["$new_gen_bus"]["type"] = "GFM"
        
        # if three_wire 
        #     data_math["bus"]["$new_gen_bus"]["vmin"] = data_math["bus"]["$new_gen_bus"]["vmin"][1:3]
        #     data_math["bus"]["$new_gen_bus"]["vmax"] = data_math["bus"]["$new_gen_bus"]["vmax"][1:3]
        #     data_math["bus"]["$new_gen_bus"]["vm"] = [1. ; 1 ; 1]
        #     data_math["bus"]["$new_gen_bus"]["va"] = [0. ; -2.094 ; 2.094]
        #     data_math["bus"]["$new_gen_bus"]["grounded"] = Bool[0, 0, 0]
        #     data_math["bus"]["$new_gen_bus"]["terminals"] = data_math["bus"]["$new_gen_bus"]["terminals"][1:3]
        # end
    end
    
    
    sbase = data_math["settings"]["sbase"]                          # p.u.
    sbace_factor = data_math["settings"]["power_scale_factor"]      # 
    vbase = [v for v in values(data_math["settings"]["vbases_default"])][1]
    vbase_factor = data_math["settings"]["voltage_scale_factor"]
    # vbase = 0.2309      # [kV]  data_math["settings"]["vbases_default"]["5"]
    Ibase = (sbase * sbace_factor) / (vbase * vbase_factor)  #[kA]
    zbase = (vbase * vbase_factor)^2 / (sbase * sbace_factor)
    vbase_max = vbase*1.1  # [V]

    Rf = 0#0.015
    Lf = 0#0.42E-3
    Cf = 0#0.33E-9
    new_branch_id = length(data_math["branch"]) + 1
    data_math["branch"]["$new_branch_id"] = Dict{String, Any}()
    data_math["branch"]["$new_branch_id"]["rate_a"] = [1, 1, 1, 1]*1000 #[Inf, Inf, Inf, Inf]
    data_math["branch"]["$new_branch_id"]["rate_a"] = [1, 1, 1, 1]*1000 #[Inf, Inf, Inf, Inf]
    data_math["branch"]["$new_branch_id"]["rate_b"] = [1, 1, 1, 1]*1000 #[Inf, Inf, Inf, Inf]
    data_math["branch"]["$new_branch_id"]["vbase"] = 0.23094
    data_math["branch"]["$new_branch_id"]["source_id"] = "pvsystem_$gen_id"
    data_math["branch"]["$new_branch_id"]["br_status"] = 1
    data_math["branch"]["$new_branch_id"]["angmin"] = [-1.0472, -1.0472, -1.0472, -1.0472]
    data_math["branch"]["$new_branch_id"]["angmax"] = [1.0472, 1.0472, 1.0472, 1.0472]
    data_math["branch"]["$new_branch_id"]["f_connections"] = [1, 2, 3, 4]
    data_math["branch"]["$new_branch_id"]["t_connections"] = [1, 2, 3, 4]
    data_math["branch"]["$new_branch_id"]["g_fr"] = zeros(4,4)
    data_math["branch"]["$new_branch_id"]["g_to"] = zeros(4,4)

    data_math["branch"]["$new_branch_id"]["index"] = new_branch_id
    data_math["branch"]["$new_branch_id"]["name"] = "inverter_branch_$gen_id"
    data_math["branch"]["$new_branch_id"]["br_r"] = diagm(Rf/zbase * ones(4))
    data_math["branch"]["$new_branch_id"]["br_x"] = diagm(Lf*2*pi*50/zbase * ones(4))
    data_math["branch"]["$new_branch_id"]["b_to"] = 0*diagm(Cf*2*pi*50*zbase * ones(4))  # TODO set this to zero for now, it should not be zero in reality.
    data_math["branch"]["$new_branch_id"]["b_fr"] = zeros(4,4)
    data_math["branch"]["$new_branch_id"]["f_bus"] = new_gen_bus
    data_math["branch"]["$new_branch_id"]["t_bus"] = old_gen_bus

    # data_math["branch"]["$new_branch_id"]["m_legs"] = 8  # add m_legs to branch data rather than gen data???

    if reconfigurable
        gen["m_legs"] = 8
    end

    if dc_link
        gen["pdcmin"] = -Inf
        gen["pdcmax"] = Inf
        # data_math["branch"]["$new_branch_id"]["pdcmin"] = -Inf
        # data_math["branch"]["$new_branch_id"]["pdcmax"] = Inf
    end

    # if three_wire 
    #     data_math["branch"]["$new_branch_id"]["f_connections"] = [1, 2, 3]
    #     data_math["branch"]["$new_branch_id"]["t_connections"] = [1, 2, 3]
    #     data_math["branch"]["$new_branch_id"]["br_r"] = diagm(Rf/zbase * ones(4))
    #     data_math["branch"]["$new_branch_id"]["br_x"] = diagm(Lf*2*pi*50/zbase * ones(4))
    #     data_math["branch"]["$new_branch_id"]["b_to"] = diagm(Cf*2*pi*50*zbase * ones(4))
    # end

    ### If the actual c_rating_a is not given, use the power rating
    if c_rating_a == 0
        pmax = gen["pmax"] * (sbase * sbace_factor)
        c_rating_a = pmax / vbase_max
        @show c_rating_a, pmax, gen["pmax"]
    end
    c_rating_a_pu = c_rating_a / Ibase
    # @show c_rating_a_pu, c_rating_a
    data_math["branch"]["$new_branch_id"]["c_rating_a"] = [c_rating_a_pu ; c_rating_a_pu[1]] * 5  # TODO the 5 multiplier is added to relax this limit, but need to be revisited
    data_math["branch"]["$new_branch_id"]["c_rating_b"] = data_math["branch"]["$new_branch_id"]["c_rating_a"]
    data_math["branch"]["$new_branch_id"]["c_rating_c"] = data_math["branch"]["$new_branch_id"]["c_rating_a"]

    return new_gen_bus, new_branch_id
end