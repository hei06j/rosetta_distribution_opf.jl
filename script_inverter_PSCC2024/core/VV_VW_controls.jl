## Inverter control Volt-var
gen_id = 1
gen = data_math["gen"]["$gen_id"]
bus_id = gen["gen_bus"]
vmin = 0.9; vmax = 1.1;
if 4 in gen["connections"]
    phases = gen["connections"][1:end-1]
    n = gen["connections"][end]
else
    phases = gen["connections"]
end

terminals = data_math["bus"]["$bus_id"]["terminals"]
vm = JuMP.@variable(model, [t in terminals, i=bus_id], base_name="vm", lower_bound=0)
JuMP.@constraint(model, [t in terminals], vm[t,bus_id]^2 == vr[t,bus_id]^2 + vi[t,bus_id]^2)


vmpp = JuMP.@variable(model, [t in 1:3, i=bus_id], base_name="vmpp", lower_bound=0)
JuMP.@constraint(model, vmpp[1,bus_id]^2 == (vr[1,bus_id]-vr[2,bus_id])^2 + (vi[1,bus_id]-vi[2,bus_id])^2)
JuMP.@constraint(model, vmpp[2,bus_id]^2 == (vr[2,bus_id]-vr[3,bus_id])^2 + (vi[2,bus_id]-vi[3,bus_id])^2)
JuMP.@constraint(model, vmpp[3,bus_id]^2 == (vr[3,bus_id]-vr[1,bus_id])^2 + (vi[3,bus_id]-vi[1,bus_id])^2)

phases = collect(1:3)
vmpn = JuMP.@variable(model, [t in 1:3, i=bus_id], base_name="vmpn", lower_bound=0)
JuMP.@constraint(model, [p in phases], vmpn[p,bus_id].^2 == (vr[p,bus_id].-vr[4,bus_id]).^2 + (vi[p,bus_id].-vi[4,bus_id]).^2)
# JuMP.@constraint(model, vmpn[2,bus_id]^2 == (vr[2,bus_id]-vr[4,bus_id])^2 + (vi[2,bus_id]-vi[4,bus_id])^2)
# JuMP.@constraint(model, vmpn[3,bus_id]^2 == (vr[3,bus_id]-vr[4,bus_id])^2 + (vi[3,bus_id]-vi[4,bus_id])^2)


# vv_curve(v, qmax) = v <= 0.95 ? qmax : (v >= 1.05 ? -qmax : -2*qmax/0.1*(v-1))
# vw_curve(v, pmax) = v <= 0.95 ? pmax : (v >= 1.05 ? 0.0 : -pmax/0.1*(v-1.05))
vv_curve(v, qmax) = v <= 0.9 ? qmax : (v >= 1.1 ? -qmax : -2*qmax/0.2*(v-1))
vw_curve(v, pmax) = v <= 0.9 ? pmax : (v >= 1.1 ? 0.0 : -pmax/0.2*(v-1.1))
JuMP.@operator(model, vv, 2, vv_curve)
JuMP.@operator(model, vw, 2, vw_curve)

### control = Volt_Watt_pn, Volt_var_pn, Volt_Watt_pp, Volt_var_pp, Volt_Watt_av, Volt_var_av
if control == "Volt_var_pn"
    JuMP.@constraint(model, [i in phases], qg[i,gen_id] == vv(vmpn[i,bus_id], gen["qmax"][i]) )
    # JuMP.@constraint(model, [i in phases], qg[i,gen_id] == vv(vm[i,bus_id]-vm[n,bus_id], gen["qmax"][i]) )

elseif control == "Volt_Watt_pn"
    # JuMP.@constraint(model, [i in phases], pg[i,gen_id] == vw(vm[i,bus_id], gen["pmax"][i]) )
    # JuMP.@constraint(model, [i in phases], pg[i,gen_id] == vw(vm[i,bus_id]-vm[n,bus_id], gen["pmax"][i]) )
    JuMP.@constraint(model, [i in phases], pg[i,gen_id] == vw(vmpn[i,bus_id], gen["pmax"][i]) )

elseif control == "Volt_var_pp"
    JuMP.@constraint(model, qg[1,gen_id] - qg[2,gen_id] == vv(vmpp[1,bus_id]/sqrt(3), gen["qmax"][1]) )
    JuMP.@constraint(model, qg[2,gen_id] - qg[3,gen_id] == vv(vmpp[2,bus_id]/sqrt(3), gen["qmax"][2]) )
    JuMP.@constraint(model, qg[3,gen_id] - qg[1,gen_id] == vv(vmpp[3,bus_id]/sqrt(3), gen["qmax"][3]) )

elseif control == "Volt_Watt_pp"
    JuMP.@constraint(model, pg[1,gen_id] - pg[2,gen_id] == vw(vmpp[1,bus_id]/sqrt(3), gen["pmax"][1]) )
    JuMP.@constraint(model, pg[2,gen_id] - pg[3,gen_id] == vw(vmpp[2,bus_id]/sqrt(3), gen["pmax"][2]) )
    JuMP.@constraint(model, pg[3,gen_id] - pg[1,gen_id] == vw(vmpp[3,bus_id]/sqrt(3), gen["pmax"][3]) )

elseif control == "Volt_var_av"
    JuMP.@constraint(model, qg[1, gen_id] == vv(sum(vm[i,bus_id] for i in ref[:bus][1]["terminals"][1:3])/3, gen["qmax"][1]) )
    JuMP.@constraint(model, qg[2, gen_id] == qg[1, gen_id])
    JuMP.@constraint(model, qg[3, gen_id] == qg[1, gen_id])

elseif control == "Volt_Watt_av"
    JuMP.@constraint(model, pg[1, gen_id] == vw(sum(vm[i,bus_id] for i in ref[:bus][1]["terminals"][1:3])/3, gen["pmax"][1]) )
    JuMP.@constraint(model, pg[2, gen_id] == pg[1, gen_id])
    JuMP.@constraint(model, pg[3, gen_id] == pg[1, gen_id])
end