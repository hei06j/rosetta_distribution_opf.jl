function add_sop_inverter_losses!(data_math, gen_id1, gen_id2; reconfigurable=false, dc_link=true)

    # ## add a branch for each inverter
    # new_gen1_bus, new_branch_id1 = RPMD.add_inverter_losses!(data_math, gen_id1; reconfigurable=reconfigurable)
    # new_gen2_bus, new_branch_id2 = RPMD.add_inverter_losses!(data_math, gen_id2; reconfigurable=reconfigurable)

    # @show new_gen1_bus, new_branch_id1
    # @show new_gen2_bus, new_branch_id2

    ## add a branch for each inverter
    f_bus = data_math["gen"]["$gen_id1"]["gen_bus"]
    t_bus = data_math["gen"]["$gen_id2"]["gen_bus"]    

    ## add a branch connecting the two inverters
    Rin = 0.015
    Cin = 0.33E-9
    Lf = 0.42E-3
    Vdc = 0.6   # [kV]
    Rin_Vdc = Rin / Vdc^2
    Cin_Vdc = Cin / Vdc^2
    zbase = 230.94^2 / 1000
    new_branch_id = length(data_math["branch"]) + 1
    branch_data = deepcopy(data_math["branch"]["$(new_branch_id-1)"])
    data_math["branch"]["$new_branch_id"] = deepcopy(branch_data)
    data_math["branch"]["$new_branch_id"]["index"] = new_branch_id
    data_math["branch"]["$new_branch_id"]["name"] = "SOP_branch__$(gen_id1)_$(gen_id2)"
    data_math["branch"]["$new_branch_id"]["source_id"] = "SOP_branch__$(gen_id1)_$(gen_id2)"
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

    ## calculate the SOP current rating
    Sbase = data_math["settings"]["sbase"]                       # MVA p.u.
    Sbace_Factor = data_math["settings"]["power_scale_factor"]   # 
    Vbase = 0.2309                                               # [kV]  data_math["settings"]["vbases_default"]["5"]
    Vbase_Factor = data_math["settings"]["voltage_scale_factor"]
    Ibase = (Sbase * Sbace_Factor) / (Vbase * Vbase_Factor)      #[kA]
    vbase_max = 253                                              # [V]

    ## select the higher pmax as the S rating of the SOP
    pmax1 = data_math["gen"]["$gen_id1"]["pmax"]
    pmax2 = data_math["gen"]["$gen_id2"]["pmax"]
    pmax = maximum([pmax1, pmax2]) * (Sbase * Sbace_Factor)
    c_rating_a = pmax / vbase_max / Ibase
    data_math["branch"]["$new_branch_id"]["c_rating_a"] = [c_rating_a ; c_rating_a[1]]
    
    if reconfigurable
        data_math["branch"]["$new_branch_id"]["m_legs"] = 12
    else
        data_math["branch"]["$new_branch_id"]["m_legs"] = 8
    end

    if dc_link 
        data_math["branch"]["$new_branch_id"]["pdcmin"] = -Inf
        data_math["branch"]["$new_branch_id"]["pdcmax"] = Inf
    end

    delete!(data_math["gen"], "$gen_id1")
    delete!(data_math["gen"], "$gen_id2")

    if isempty([i for (i, gen) in data_math["gen"] if gen["gen_bus"] == f_bus])
        data_math["bus"]["$f_bus"]["bus_type"] = 1
    end
    if isempty([i for (i, gen) in data_math["gen"] if gen["gen_bus"] == t_bus])
        data_math["bus"]["$t_bus"]["bus_type"] = 1
    end

    ## set the generations output to zero, so that no exgenous input exists
    # for id in [gen_id1, gen_id2]
    #     # data_math["gen"]["$id"]["pmin"] *= 0
    #     # data_math["gen"]["$id"]["pmax"] *= 0
    #     # data_math["gen"]["$id"]["qmin"] *= 0
    #     # data_math["gen"]["$id"]["qmax"] *= 0
    #     if reconfigurable
    #         data_math["gen"]["$id"]["m_legs"] = 8
    #     end
    # end
    # data_math["gen"]["$gen_id1"]["inverter_branch"] = new_branch_id1
    # data_math["gen"]["$gen_id2"]["inverter_branch"] = new_branch_id2
    # delete!(data_math["gen"], "$gen_id1")
    # delete!(data_math["gen"], "$gen_id2")
end