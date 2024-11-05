function get_coordinates(eng, bus_coordinates_file)
    for line in eachline(bus_coordinates_file)
        words = split(line, (',', ','), keepempty=false)
        eng["bus"][words[1]]["lon"] = parse(Float64, words[2])
        eng["bus"][words[1]]["lat"] = parse(Float64, words[3])
    end
end

function construct_graph_eng(eng)
    # graph = _G.AbstractGraph(length(eng["bus"]));
    graph = _G.SimpleGraph(length(eng["bus"]));
    map = Dict(i => c for (c, (i, bus)) in enumerate(eng["bus"]))

    for (l, line) in eng["line"]
        _G.add_edge!(graph, map[line["f_bus"]], map[line["t_bus"]]);
    end
    #TODO support multiwinding transformers
    if haskey(eng, "transformer")
        for (l, transformer) in eng["transformer"]
            _G.add_edge!(graph, map[transformer["bus"][1]], map[transformer["bus"][2]]);
        end
    end
    if haskey(eng, "switch")
        for (s, switch) in eng["switch"]
            _G.add_edge!(graph, map[switch["f_bus"]], map[switch["t_bus"]]);
        end
    end

    return graph
end

function construct_graph_math(math)
    graph = _G.SimpleGraph(length(math["bus"]));
    for (l, branch) in math["branch"]
        _G.add_edge!(graph, branch["f_bus"], branch["t_bus"]);
    end

    for (l, transformer) in math["transformer"]
        _G.add_edge!(graph, transformer["f_bus"], transformer["t_bus"]);
    end
    return graph
end


function plot_network(graph, eng, math)
    nbuses = (length(eng["bus"]))
    nbranches = (length(eng["line"]))
    locs_x = zeros(nbuses)
    locs_y = zeros(nbuses)
    nodelabel= ["" for i in 1:nbuses]
    edgelabel= ["" for i in 1:nbranches]


    for (i, bus) in math["bus"]
        busname = bus["name"]
        if startswith(busname, "_virtual_bus.voltage_source.")
            voltage_source_name = split(busname, "_virtual_bus.voltage_source.")[2]
            @show voltage_source_name
            voltage_source_bus = eng["voltage_source"][voltage_source_name]["bus"]
            @show voltage_source_bus
            @show eng["bus"][voltage_source_bus]
            locs_x[bus["bus_i"]] = eng["bus"][voltage_source_bus]["lon"]
            locs_y[bus["bus_i"]] = eng["bus"][voltage_source_bus]["lat"]
        else
            locs_x[bus["bus_i"]] = eng["bus"][busname]["lon"]
            locs_y[bus["bus_i"]] = eng["bus"][busname]["lat"]
            nodelabel[bus["bus_i"]] = busname
        end
    end
    # _GP.gplot(graph, edgelabel=edgelabel, nodelabel=nodelabel)
    _GP.gplot(graph)
end


function plot_network_graph(data_eng, data_math, bus_coordinates_file)
    get_coordinates(data_eng, bus_coordinates_file)
    graph = construct_graph_eng(data_eng)
    # graph = RPMD.construct_graph_math(data_math_mx)
    return plot_network(graph, data_eng, data_math)
end