
function add_inverter_losses(data_math, gen_id; GFM=false, three_wire=false)
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
    
    Rf = 0.015
    Lf = 0.42E-3
    Cf = 0.33E-9
    zbase = 230.94^2 / 1000
    new_branch_id = length(data_math["branch"]) + 1
    data_math["branch"]["$new_branch_id"] = deepcopy(data_math["branch"]["$(new_branch_id-1)"])
    data_math["branch"]["$new_branch_id"]["index"] = new_branch_id
    data_math["branch"]["$new_branch_id"]["name"] = "GFL_internal_z_$gen_id"
    data_math["branch"]["$new_branch_id"]["br_r"] = diagm(Rf/zbase * ones(4))
    data_math["branch"]["$new_branch_id"]["br_x"] = diagm(Lf*2*pi*50/zbase * ones(4))
    data_math["branch"]["$new_branch_id"]["b_to"] = diagm(Cf*2*pi*50*zbase * ones(4))
    data_math["branch"]["$new_branch_id"]["f_bus"] = new_gen_bus
    data_math["branch"]["$new_branch_id"]["t_bus"] = old_gen_bus

    # if three_wire 
    #     data_math["branch"]["$new_branch_id"]["f_connections"] = [1, 2, 3]
    #     data_math["branch"]["$new_branch_id"]["t_connections"] = [1, 2, 3]
    #     data_math["branch"]["$new_branch_id"]["br_r"] = diagm(Rf/zbase * ones(4))
    #     data_math["branch"]["$new_branch_id"]["br_x"] = diagm(Lf*2*pi*50/zbase * ones(4))
    #     data_math["branch"]["$new_branch_id"]["b_to"] = diagm(Cf*2*pi*50*zbase * ones(4))
    # end

    #### Irating = 52 A 
    #### Ibase = (Sbase * Sbace_Factor) / (Vbase * Vbase_Factor) = (1 * 1000) / (0.2309 * 1000) = 4.33
    #### Irating_pu = 52 / 4.33 = 12
    # a = 12  
    # Ibase = 4.33
    # vbase_max = 253
    # # c_rating_a = 3*gen["pmax"][1] * 1000 / (vbase_max*3) / Ibase
    # c_rating_a = 3*gen["smax"][1] * 1000 / (vbase_max*3) / Ibase
    # # c_rating_a = c_rating_a < 10. ? 10. : c_rating_a
    # data_math["branch"]["$new_branch_id"]["c_rating_a"] = c_rating_a * ones(4)
    # data_math["branch"]["$new_branch_id"]["c_rating_b"] = c_rating_a * ones(4)
    # data_math["branch"]["$new_branch_id"]["c_rating_c"] = c_rating_a * ones(4)

end