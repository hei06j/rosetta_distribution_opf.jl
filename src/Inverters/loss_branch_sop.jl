# function add_sop_inverter_losses!(data_math, gen_id1, gen_id2; c_rating_a=0, reconfigurable=false, dc_link=true)
#     ## add a branch for each inverter
#     f_bus = data_math["gen"]["$gen_id1"]["gen_bus"]
#     t_bus = data_math["gen"]["$gen_id2"]["gen_bus"]    

#     sbase = data_math["settings"]["sbase"]                          # p.u.
#     sbace_factor = data_math["settings"]["power_scale_factor"]      # 
#     vbase = [v for v in values(data_math["settings"]["vbases_default"])][1]
#     vbase_factor = data_math["settings"]["voltage_scale_factor"]
#     Ibase = (sbase * sbace_factor) / (vbase * vbase_factor)  #[kA]
#     zbase = (vbase * vbase_factor)^2 / (sbase * sbace_factor)
#     vbase_max = vbase*1.1  # [V]
        
#     ## add a branch connecting the two inverters
#     Rin = 0.015
#     Cin = 0.33E-9
#     Lf = 0.42E-3
#     Vdc = 0.6   # [kV]
#     Rin_Vdc = Rin / Vdc^2
#     Cin_Vdc = Cin / Vdc^2
#     new_branch_id = length(data_math["branch"]) + 1
#     branch_data = deepcopy(data_math["branch"]["$(new_branch_id-1)"])
#     data_math["branch"]["$new_branch_id"] = deepcopy(branch_data)
#     data_math["branch"]["$new_branch_id"]["index"] = new_branch_id
#     data_math["branch"]["$new_branch_id"]["name"] = "SOP_branch__$(gen_id1)_$(gen_id2)"
#     data_math["branch"]["$new_branch_id"]["source_id"] = "SOP_branch__$(gen_id1)_$(gen_id2)"
#     # data_math["branch"]["$new_branch_id"]["f_connections"] = branch_data["f_connections"][[1,4]]
#     # data_math["branch"]["$new_branch_id"]["t_connections"] = branch_data["t_connections"][[1,4]]
#     # data_math["branch"]["$new_branch_id"]["rate_a"] = [Inf, Inf] # branch_data["rate_a"][[1,4]]
#     # data_math["branch"]["$new_branch_id"]["rate_b"] = [Inf, Inf] # branch_data["rate_b"][[1,4]]
#     # data_math["branch"]["$new_branch_id"]["rate_c"] = [Inf, Inf] # branch_data["rate_c"][[1,4]]
#     # data_math["branch"]["$new_branch_id"]["c_rating_a"] = branch_data["c_rating_a"][[1,4]]
#     # data_math["branch"]["$new_branch_id"]["c_rating_b"] = branch_data["c_rating_a"][[1,4]] #branch_data["c_rating_b"][[1,4]]
#     # data_math["branch"]["$new_branch_id"]["c_rating_c"] = branch_data["c_rating_a"][[1,4]] #branch_data["c_rating_c"][[1,4]]
#     # data_math["branch"]["$new_branch_id"]["angmin"] = branch_data["angmin"][[1,4]]
#     # data_math["branch"]["$new_branch_id"]["angmax"] = branch_data["angmax"][[1,4]]
#     data_math["branch"]["$new_branch_id"]["br_r"] = diagm(Rin_Vdc/zbase * ones(4))
#     data_math["branch"]["$new_branch_id"]["br_x"] = zeros(4,4) #diagm(Lf*2*pi*50/zbase * ones(2))
#     data_math["branch"]["$new_branch_id"]["g_fr"] = zeros(4,4)
#     data_math["branch"]["$new_branch_id"]["g_to"] = zeros(4,4)
#     data_math["branch"]["$new_branch_id"]["b_fr"] = diagm(Cin_Vdc*2*pi*50*zbase * ones(4))
#     data_math["branch"]["$new_branch_id"]["b_to"] = diagm(Cin_Vdc*2*pi*50*zbase * ones(4))
#     data_math["branch"]["$new_branch_id"]["f_bus"] = f_bus
#     data_math["branch"]["$new_branch_id"]["t_bus"] = t_bus


#     ### If the actual c_rating_a is not given, select the higher pmax as the power rating of the SOP
#     if c_rating_a == 0
#         pmax1 = data_math["gen"]["$gen_id1"]["pmax"]
#         pmax2 = data_math["gen"]["$gen_id2"]["pmax"]
#         pmax = maximum([pmax1, pmax2]) * (sbase * sbace_factor)
#         c_rating_a = pmax / vbase_max #/ Ibase
#     end

#     c_rating_a_pu = c_rating_a / Ibase
#     data_math["branch"]["$new_branch_id"]["c_rating_a"] = [c_rating_a_pu ; c_rating_a_pu[1]]
    
#     if reconfigurable
#         data_math["branch"]["$new_branch_id"]["m_legs"] = 12
#     else
#         data_math["branch"]["$new_branch_id"]["m_legs"] = 8
#     end

#     if dc_link 
#         data_math["branch"]["$new_branch_id"]["pdcmin"] = -Inf
#         data_math["branch"]["$new_branch_id"]["pdcmax"] = Inf
#     end

#     delete!(data_math["gen"], "$gen_id1")
#     delete!(data_math["gen"], "$gen_id2")

#     if isempty([i for (i, gen) in data_math["gen"] if gen["gen_bus"] == f_bus])
#         data_math["bus"]["$f_bus"]["bus_type"] = 1
#     end
#     if isempty([i for (i, gen) in data_math["gen"] if gen["gen_bus"] == t_bus])
#         data_math["bus"]["$t_bus"]["bus_type"] = 1
#     end
    
#     ## set the generations output to zero, so that no exgenous input exists
#     # for id in [gen_id1, gen_id2]
#     #     # data_math["gen"]["$id"]["pmin"] *= 0
#     #     # data_math["gen"]["$id"]["pmax"] *= 0
#     #     # data_math["gen"]["$id"]["qmin"] *= 0
#     #     # data_math["gen"]["$id"]["qmax"] *= 0
#     #     if reconfigurable
#     #         data_math["gen"]["$id"]["m_legs"] = 8
#     #     end
#     # end
#     # data_math["gen"]["$gen_id1"]["inverter_branch"] = new_branch_id1
#     # data_math["gen"]["$gen_id2"]["inverter_branch"] = new_branch_id2
#     # delete!(data_math["gen"], "$gen_id1")
#     # delete!(data_math["gen"], "$gen_id2")
# end


function add_sop_inverter_losses_v2!(data_math, f_bus, t_bus; c_rating_a=0, reconfigurable=false, dc_link=true)
    sbase = data_math["settings"]["sbase"]                          # p.u.
    sbace_factor = data_math["settings"]["power_scale_factor"]      # 
    # vbase = [v for v in values(data_math["settings"]["vbases_default"])][2]
    vbase_factor = data_math["settings"]["voltage_scale_factor"]
    vbase = 0.2309      # [kV]  data_math["settings"]["vbases_default"]["5"]
    Ibase = (sbase * sbace_factor) / (vbase * vbase_factor)  #[kA]
    zbase = (vbase * vbase_factor)^2 / (sbase * sbace_factor)
    vbase_max = vbase*1.1  # [V]
        
    ## add a branch connecting the two inverters
    Rin = 0.015
    Cin = 0.33E-9
    Lf = 0.42E-3
    Vdc = 0.6   # [kV]
    Rin_Vdc = Rin / Vdc^2
    Cin_Vdc = Cin / Vdc^2
    new_branch_id = length(data_math["branch"]) + 1
    branch_data = deepcopy(data_math["branch"]["1"])
    data_math["branch"]["$new_branch_id"] = deepcopy(branch_data)
    data_math["branch"]["$new_branch_id"]["index"] = new_branch_id
    data_math["branch"]["$new_branch_id"]["name"] = "SOP_branch_$(f_bus)_$(t_bus)"
    data_math["branch"]["$new_branch_id"]["source_id"] = "SOP_branch_$(f_bus)_$(t_bus)"
    # data_math["branch"]["$new_branch_id"]["f_connections"] = branch_data["f_connections"][[1,4]]
    # data_math["branch"]["$new_branch_id"]["t_connections"] = branch_data["t_connections"][[1,4]]
    # data_math["branch"]["$new_branch_id"]["rate_a"] = [Inf, Inf] # branch_data["rate_a"][[1,4]]
    # data_math["branch"]["$new_branch_id"]["rate_b"] = [Inf, Inf] # branch_data["rate_b"][[1,4]]
    # data_math["branch"]["$new_branch_id"]["rate_c"] = [Inf, Inf] # branch_data["rate_c"][[1,4]]
    # data_math["branch"]["$new_branch_id"]["c_rating_a"] = branch_data["c_rating_a"][[1,4]]
    # data_math["branch"]["$new_branch_id"]["c_rating_b"] = branch_data["c_rating_a"][[1,4]] #branch_data["c_rating_b"][[1,4]]
    # data_math["branch"]["$new_branch_id"]["c_rating_c"] = branch_data["c_rating_a"][[1,4]] #branch_data["c_rating_c"][[1,4]]
    # data_math["branch"]["$new_branch_id"]["angmin"] = branch_data["angmin"][[1,4]]
    # data_math["branch"]["$new_branch_id"]["angmax"] = branch_data["angmax"][[1,4]]
    data_math["branch"]["$new_branch_id"]["br_r"] = diagm(Rin_Vdc/zbase * ones(4))
    data_math["branch"]["$new_branch_id"]["br_x"] = zeros(4,4) #diagm(Lf*2*pi*50/zbase * ones(2))
    data_math["branch"]["$new_branch_id"]["g_fr"] = zeros(4,4)
    data_math["branch"]["$new_branch_id"]["g_to"] = zeros(4,4)
    data_math["branch"]["$new_branch_id"]["b_fr"] = diagm(Cin_Vdc*2*pi*50*zbase * ones(4))
    data_math["branch"]["$new_branch_id"]["b_to"] = diagm(Cin_Vdc*2*pi*50*zbase * ones(4))
    data_math["branch"]["$new_branch_id"]["f_bus"] = f_bus
    data_math["branch"]["$new_branch_id"]["t_bus"] = t_bus


    ### If the actual c_rating_a is not given, select the higher pmax as the power rating of the SOP
    if c_rating_a == 0
        # pmax1 = data_math["gen"]["$gen_id1"]["pmax"]
        # pmax2 = data_math["gen"]["$gen_id2"]["pmax"]
        # pmax = maximum([pmax1, pmax2]) * (sbase * sbace_factor)
        pmax = 1e3 * ones(3) * (sbase * sbace_factor)
        c_rating_a = pmax / vbase_max #/ Ibase
    end

    c_rating_a_pu = c_rating_a / Ibase
    data_math["branch"]["$new_branch_id"]["c_rating_a"] = [c_rating_a_pu ; c_rating_a_pu[1]]
    
    if reconfigurable
        data_math["branch"]["$new_branch_id"]["m_legs"] = 12
    else
        data_math["branch"]["$new_branch_id"]["m_legs"] = 8
    end

    if dc_link 
        data_math["branch"]["$new_branch_id"]["pdcmin"] = 0
        data_math["branch"]["$new_branch_id"]["pdcmax"] = 10
    end


    if isempty([i for (i, gen) in data_math["gen"] if gen["gen_bus"] == f_bus])
        data_math["bus"]["$f_bus"]["bus_type"] = 1
    end
    if isempty([i for (i, gen) in data_math["gen"] if gen["gen_bus"] == t_bus])
        data_math["bus"]["$t_bus"]["bus_type"] = 1
    end
    
    return new_branch_id
end