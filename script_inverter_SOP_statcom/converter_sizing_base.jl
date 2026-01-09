using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using JuMP
using Ipopt
using LinearAlgebra: diag, diagm
using LaTeXStrings
using DataFrames
using CSV
using Plots
using Printf
using Statistics

const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels

PMD.silence!()

data_path = "./data/inverter_4w_wye_unbalanced_loads_2bus.dss"

function build_inverter_case(data_eng)
    data_math = PMD.transform_data_model(data_eng, kron_reduce=false, phase_project=false)
    # Add inverters but keep ratings fixed (no sizing/optimization)
    pv_gen_ids = [i for (i, gen) in data_math["gen"] if !occursin("source", gen["name"])]
    for gen_id in pv_gen_ids
        RPMD.add_inverter_losses!(data_math, gen_id; c_rating_a=30*ones(3))
        # fixed ratings only
        data_math["gen"][gen_id]["srating_min"] = 50
        data_math["gen"][gen_id]["srating_max"] = 50
        data_math["gen"][gen_id]["pdcrating_min"] = 50
        data_math["gen"][gen_id]["pdcrating_max"] = 50

        #=
        # No converter
        data_math["gen"][gen_id]["srating_min"] = 0
        data_math["gen"][gen_id]["srating_max"] = 0
        data_math["gen"][gen_id]["pdcrating_min"] = 0
        data_math["gen"][gen_id]["pdcrating_max"] = 0

        # Force AC power injection to zero
        data_math["gen"][gen_id]["pmin"] .= 0.0
        data_math["gen"][gen_id]["pmax"] .= 0.0
        data_math["gen"][gen_id]["qmin"] .= 0.0
        data_math["gen"][gen_id]["qmax"] .= 0.0

        # Zero current limits
        data_math["gen"][gen_id]["c_rating_a"] = [0.0, 0.0, 0.0]
        =#
    end
    return data_math
end

function build_data_math(data_path; sbase=1.0)
    data_eng = PMD.parse_file(data_path, transformations=[PMD.transform_loops!])
    data_eng["settings"]["sbase_default"] = sbase
    data_eng["voltage_source"]["source"]["rs"] *= 0
    data_eng["voltage_source"]["source"]["xs"] *= 0

    data_math = build_inverter_case(data_eng)
    # set the source inverter to 0 for forward case
    data_math["gen"]["1"]["pmax"] .= 0
    data_math["gen"]["1"]["qmax"] .= 0
    data_math["gen"]["1"]["qmin"] .= 0

    # fix bus voltage limits (no sizing)
    for (i,bus) in data_math["bus"]
        bus["vmin"] .= 0.9
        bus["vmax"] .= 1.1
    end

    # scale voltage and current bases
    sbase_factor = data_math["settings"]["power_scale_factor"]
    vbase = first(values(data_math["settings"]["vbases_default"]))
    vbase_factor = data_math["settings"]["voltage_scale_factor"]
    Ibase = (sbase * sbase_factor) / (vbase * vbase_factor)

    return data_math, Ibase
end

# -------------------- MAIN BASE CASE --------------------
setting = Dict("dc_link" => true)
sbase = 1.0
data_math, Ibase = build_data_math(data_path; sbase=sbase)

# replicate for 24-hour horizon
PMD.add_start_vrvi!(data_math)
mn_data = IM.replicate(data_math, 24, PMD._pmd_global_keys)
PMD.make_multinetwork(mn_data)
for (n, nw) in mn_data["nw"]
    nw["per_unit"] = true
end

# -------------------- LOAD ASSIGNMENT --------------------
pd_peak = 4.0
qd_peak = 1.0

normalized_load_a = [9.231668355, 8.708293005, 8.284840935, 8.113415685, 8.12055168, 8.388332595, 9.013518645000001, 9.747059055000001, 2.7997820064, 2.93546004, 2.914690196, 3.041727568, 3.169100912, 3.42147692, 3.6439375, 3.86539076, 10.0, 9.95703895, 9.69301121, 9.350277380000001, 9.061107309999999, 8.461982129999999, 7.7090232599999995, 6.80178199]

normalized_load_b = [4.308111899, 4.063870069, 3.866259103, 3.7862606530000003, 3.7895907839999996, 3.9145552109999997, 4.206308701, 4.548627559000001, 6.99945516, 7.338650100000001, 7.28672549, 7.604318920000001, 7.922752280000001, 8.5536923, 9.10984375, 9.6634769, 0.4, 0.39828155800000004, 0.38772044840000003, 0.37401109520000003, 0.3624442924, 0.33847928520000004, 0.3083609304, 0.2720712796]

normalized_load_c = [0.615444557, 0.580552867, 0.552322729, 0.540894379, 0.541370112, 0.559222173, 0.600901243, 0.649803937, 4.199673096, 4.40319006, 4.372035294, 4.562591352, 4.753651368, 5.13221538, 5.46590625, 5.79808614, 15.0, 14.935558425, 14.539516814999999, 14.02541607, 13.591660964999999, 12.692973195, 11.56353489, 10.202672985]

pd_load_a = [pd_peak * val for val in normalized_load_a]
qd_load_a = [qd_peak * val for val in normalized_load_a]

pd_load_b = [pd_peak * val for val in normalized_load_b]
qd_load_b = [qd_peak * val for val in normalized_load_b]

pd_load_c = [pd_peak * val for val in normalized_load_c]
qd_load_c = [qd_peak * val for val in normalized_load_c]

for time_step in 1:24
    nw_key = string(time_step)
    load_data = mn_data["nw"][nw_key]["load"]
    for load_id in keys(load_data)
        if load_id == "1"
            load_data[load_id]["pd"] = pd_load_a[time_step]
            load_data[load_id]["qd"] = qd_load_a[time_step]
        elseif load_id == "2"
            load_data[load_id]["pd"] = pd_load_b[time_step]
            load_data[load_id]["qd"] = qd_load_b[time_step]
        elseif load_id == "3"
            load_data[load_id]["pd"] = pd_load_c[time_step]
            load_data[load_id]["qd"] = qd_load_c[time_step]
        end
    end
end


model = PMD.instantiate_mc_model(mn_data, PMD.IVRENPowerModel, RPMD.build_mn_mc_opf_sizing; setting=setting)

# objective_mc_min_sizing(model; a_s=339.96, a_pdc=69.72)

result = PMD.optimize_model!(model, optimizer=ipopt_solver)

Vdc_nominal = 1100  

for n in 1:24
    gen_dict = result["solution"]["nw"][string(n)]["gen"]["1"]

    # Only compute if pdcrating exists
    if haskey(gen_dict, "pdcrating")
        I_dc_2w = gen_dict["pdcrating"] / Vdc_nominal
        gen_dict["idc2w"] = I_dc_2w   # Add to the dictionary
    end

    println("Network $n: ", gen_dict)
end


# ====================== Plots ========================

# Compute phase currents and sequence components over 24 hours

# --- assume pd_load_a/b/c and qd_load_a/b/c are already defined ---
n_hours = length(pd_load_a)
V_nom = 1.0
a = cis(2*pi/3)

# phase voltage angles (A = 0, B = -120°, C = +120°)
θa, θb, θc = 0.0, -2*pi/3, 2*pi/3

# arrays
Ia = zeros(ComplexF64, n_hours)
Ib = zeros(ComplexF64, n_hours)
Ic = zeros(ComplexF64, n_hours)
I0 = zeros(ComplexF64, n_hours)
Ipos = zeros(ComplexF64, n_hours)   # positive sequence (I1)
Ineg = zeros(ComplexF64, n_hours)   # negative sequence (I2)
CUF = zeros(n_hours)

# compute currents and sequences
for t in 1:n_hours
    Sa = pd_load_a[t] + im*qd_load_a[t]
    Sb = pd_load_b[t] + im*qd_load_b[t]
    Sc = pd_load_c[t] + im*qd_load_c[t]

    Va = V_nom * cis(θa)
    Vb = V_nom * cis(θb)
    Vc = V_nom * cis(θc)

    Ia[t] = conj(Sa / Va)
    Ib[t] = conj(Sb / Vb)
    Ic[t] = conj(Sc / Vc)

    # use the standard Fortescue ordering: positive = (Ia + a*Ib + a^2*Ic)/3
    I0[t]   = (Ia[t] + Ib[t] + Ic[t]) / 3
    Ipos[t] = (Ia[t] + a*Ib[t] + a^2*Ic[t]) / 3   # positive sequence
    Ineg[t] = (Ia[t] + a^2*Ib[t] + a*Ic[t]) / 3   # negative sequence

    CUF[t] = abs(Ineg[t]) / (abs(Ipos[t]) + 1e-12)
end

# quick diagnostic prints
@printf("\nBalanced test (sanity check) with Ia=1∠0, Ib=1∠-120, Ic=1∠+120:\n")
Ia_test, Ib_test, Ic_test = 1.0*cis(0.0), 1.0*cis(-2*pi/3), 1.0*cis(2*pi/3)
I0_t = (Ia_test + Ib_test + Ic_test)/3
Ipos_t = (Ia_test + a*Ib_test + a^2*Ic_test)/3
Ineg_t = (Ia_test + a^2*Ib_test + a*Ic_test)/3
@printf(" |I0|=%.6f |Ipos|=%.6f |Ineg|=%.6f   (expected Ipos≈1)\n\n", abs(I0_t), abs(Ipos_t), abs(Ineg_t))

@printf("Sample diagnostics (first 6 hours):\n")
@printf("t |   |Ia|   |Ib|   |Ic|   |Ipos|   |Ineg|   CUF\n")
for t in 1:min(6,n_hours)
    @printf("%2d | %7.4f %7.4f %7.4f %9.6f %9.6f %8.6f\n",
            t, abs(Ia[t]), abs(Ib[t]), abs(Ic[t]), abs(Ipos[t]), abs(Ineg[t]), CUF[t])
end
@printf("\nCUF range: (%.6g, %.6g)\n\n", minimum(CUF), maximum(CUF))

# --- Plot 1: Phase currents ---
p1 = plot(1:n_hours, abs.(Ia), label="|Ia|", lw=2)
plot!(1:n_hours, abs.(Ib), label="|Ib|", lw=2)
plot!(1:n_hours, abs.(Ic), label="|Ic|", lw=2)
xlabel!("Hour"); ylabel!("Phase Current (p.u.)")
title!("Phase Currents Over 24 Hours")
savefig(p1, "phase_currents.png")
println("Saved phase_currents.png")

# --- Plot 2: Sequence components (magnitudes) ---
p2 = plot(1:n_hours, abs.(Ipos), label="Positive Seq (Ipos)", lw=2)
plot!(1:n_hours, abs.(Ineg), label="Negative Seq (Ineg)", lw=2)
plot!(1:n_hours, abs.(I0), label="Zero Seq (I0)", lw=2)
xlabel!("Hour"); ylabel!("Sequence Magnitude (p.u.)")
title!("Sequence Components Over 24 Hours (corrected)")
savefig(p2, "sequence_components.png")
println("Saved sequence_components.png")

# --- Plot 3: CUF (negative / positive) ---
ymin, ymax = minimum(CUF), maximum(CUF)
# handle near-constant case gracefully
if isapprox(ymin, ymax; atol=1e-12)
    center = (ymin + ymax)/2
    span = max(1e-4, center*0.1 + 1e-6)
    ylims = (max(0.0, center - span), center + span)
else
    padding = 0.1*(ymax - ymin)
    ylims = (max(0.0, ymin - padding), ymax + padding)
end

p3 = plot(1:n_hours, CUF,
          lw=2,
          xlabel="Hour", ylabel="Current Unbalance Factor (|Ineg|/|Ipos|)",
          title="Current Unbalance Factor Over 24 Hours",
          legend=false, grid=true, ylims=ylims,
          yformatter = y -> string(round(y, digits=6))
)
savefig(p3, "current_unbalance_factor.png")
println("Saved current_unbalance_factor.png")

println("\nDone. Inspect the three PNG files and the printed diagnostics above.")


# ====================== Results ========================
# Full 24-hour extraction, CSV export, plotting and sequence/unbalance analysis
# Assumes `data_math` and `result` are already in the workspace (from your earlier parse/solve).
#
# Requirements: Plots, DataFrames, CSV, Statistics

# ---------------------------
# Settings
# ---------------------------
n_hours = 24                 # number of timesteps in result (1..24)
out_folder = "."             # change if you want to save somewhere else
eps = 1e-9                   # small number to protect divisions

# ---------------------------
# 0. Discover keys / maps (automatic)
# ---------------------------
@assert haskey(result, "solution") && haskey(result["solution"], "nw") "result structure unexpected"

# Build mapping: bus_name => bus_key (string) and key => name
bus_name_to_key = Dict{String,String}()
key_to_bus_name = Dict{String,String}()
if haskey(data_math, "bus")
    for (k, v) in data_math["bus"]
        if haskey(v, "name")
            bus_name_to_key[v["name"]] = k
            key_to_bus_name[k] = v["name"]
        else
            # fallback: use numeric string as name
            key_to_bus_name[k] = "bus_$k"
            bus_name_to_key["bus_$k"] = k
        end
    end
else
    error("data_math does not contain 'bus' key.")
end

# Branch labels from data_math (if available) to improve plot legends
branch_key_to_label = Dict{String,String}()
if haskey(data_math, "line")
    # try to read bus1/bus2 fields from lines
    for (k, ln) in data_math["line"]
        b1 = haskey(ln, "bus1") ? string(ln["bus1"]) : nothing
        b2 = haskey(ln, "bus2") ? string(ln["bus2"]) : nothing
        if b1 !== nothing && b2 !== nothing && haskey(key_to_bus_name, b1) && haskey(key_to_bus_name, b2)
            branch_key_to_label[k] = "$(key_to_bus_name[b1]) → $(key_to_bus_name[b2])"
        else
            branch_key_to_label[k] = "branch_$k"
        end
    end
else
    # fallback: just use numeric keys
    println("Warning: data_math['line'] not found; branch labels will use numeric keys.")
end

# Determine which keys exist in result for a single timestep (to validate)
sample_nw = result["solution"]["nw"][string(1)]
@show keys(sample_nw)     # helpful debug: shows what components exist (bus, branch, load, gen, ...)

# Branch keys available in the solution (strings)
branch_keys = collect(keys(sample_nw["branch"]))   # e.g. ["1","2"]
bus_keys_in_result = collect(keys(sample_nw["bus"])) # e.g. ["1","2","3"]

# helper to pick the "load" bus (prefer bus_type==1 in data_math if present)
load_bus_key = nothing
for (k, v) in data_math["bus"]
    if haskey(v, "bus_type") && v["bus_type"] == 1
        load_bus_key = k
        break
    end
end
if load_bus_key === nothing
    # fallback: pick key corresponding to name "b2" if exists
    if haskey(bus_name_to_key, "b2")
        load_bus_key = bus_name_to_key["b2"]
    else
        # fallback to first bus key in result
        load_bus_key = first(bus_keys_in_result)
    end
end
println("Using load bus key: ", load_bus_key, " (", key_to_bus_name[load_bus_key], ")")

# find branch supplying the load: branch where one endpoint equals load bus (look in data_math.line)
supply_branch_key = nothing
if haskey(data_math, "line")
    for (k, ln) in data_math["line"]
        # ln["bus1"], ln["bus2"] might be numeric indices or strings; coerce to string
        b1 = haskey(ln, "bus1") ? string(ln["bus1"]) : nothing
        b2 = haskey(ln, "bus2") ? string(ln["bus2"]) : nothing
        if b1 == load_bus_key || b2 == load_bus_key
            supply_branch_key = k
            break
        end
    end
end
if supply_branch_key === nothing
    # fallback to first branch key
    supply_branch_key = branch_keys[1]
end
println("Using supply branch key: ", supply_branch_key, " label=", get(branch_key_to_label, supply_branch_key, supply_branch_key))

# ---------------------------
# 1. Containers for extracted data
# ---------------------------
# We'll iterate over all bus keys present in data_math (numeric string keys)
all_bus_keys = collect(keys(data_math["bus"]))  # e.g. ["1","2","3"]
all_branch_keys = branch_keys                   # keys available in results (e.g. ["1","2"])

# Voltage complex (phase-to-neutral) as arrays per bus
V_a = Dict{String, Vector{ComplexF64}}()
V_b = Dict{String, Vector{ComplexF64}}()
V_c = Dict{String, Vector{ComplexF64}}()
for bk in all_bus_keys
    V_a[bk] = zeros(ComplexF64, n_hours)
    V_b[bk] = zeros(ComplexF64, n_hours)
    V_c[bk] = zeros(ComplexF64, n_hours)
end

# Branch currents (from-end)
I_a = Dict{String, Vector{ComplexF64}}()
I_b = Dict{String, Vector{ComplexF64}}()
I_c = Dict{String, Vector{ComplexF64}}()
for br in all_branch_keys
    I_a[br] = zeros(ComplexF64, n_hours)
    I_b[br] = zeros(ComplexF64, n_hours)
    I_c[br] = zeros(ComplexF64, n_hours)
end

# Reactive power qf per branch per phase (store as Float64)
Qf_a = Dict{String, Vector{Float64}}()
Qf_b = Dict{String, Vector{Float64}}()
Qf_c = Dict{String, Vector{Float64}}()
for br in all_branch_keys
    Qf_a[br] = zeros(n_hours); Qf_b[br] = zeros(n_hours); Qf_c[br] = zeros(n_hours)
end

# ---------------------------
# 2. Extract per-hour values from result
# ---------------------------
for t in 1:n_hours
    nw = result["solution"]["nw"][string(t)]

    # BUSES: construct phase-to-neutral phasors using neutral (4th element)
    for bk in all_bus_keys
        bus_res = nw["bus"][bk]      # structure with "vr" and "vi" arrays
        vr = bus_res["vr"]
        vi = bus_res["vi"]

        # safety check lengths
        if length(vr) < 4 || length(vi) < 4
            error("bus $bk: expected vr/vi length >=4 (A,B,C,N), got $(length(vr))/$(length(vi))")
        end

        Vn = vr[4] + im*vi[4]  # neutral phasor (phase-to-neutral reference)
        Va_ph = (vr[1] + im*vi[1]) - Vn
        Vb_ph = (vr[2] + im*vi[2]) - Vn
        Vc_ph = (vr[3] + im*vi[3]) - Vn

        V_a[bk][t] = Va_ph
        V_b[bk][t] = Vb_ph
        V_c[bk][t] = Vc_ph
    end

    # BRANCHES: use cr_fr and ci_fr for 'from' complex currents (per-phase)
    for br in all_branch_keys
        br_res = nw["branch"][br]

        # prefer "cr_fr" and "ci_fr"; some variants might be named differently -> check alternatives
        if haskey(br_res, "cr_fr") && haskey(br_res, "ci_fr")
            cr = br_res["cr_fr"]
            ci = br_res["ci_fr"]
            # check lengths
            if length(cr) < 3 || length(ci) < 3
                error("branch $br: expected cr_fr/ci_fr length >=3, got $(length(cr))/$(length(ci))")
            end
            I_a[br][t] = cr[1] + im*ci[1]
            I_b[br][t] = cr[2] + im*ci[2]
            I_c[br][t] = cr[3] + im*ci[3]
        elseif haskey(br_res, "cr_to") && haskey(br_res, "ci_to")
            # fallback using to-end if from-end not present
            cr = br_res["cr_to"]; ci = br_res["ci_to"]
            I_a[br][t] = cr[1] + im*ci[1]
            I_b[br][t] = cr[2] + im*ci[2]
            I_c[br][t] = cr[3] + im*ci[3]
        else
            error("branch $br: expected cr_fr/ci_fr or cr_to/ci_to fields")
        end

        # Reactive power qf: if qf present (vector length >=3), store per phase; else attempt mean fallback
        if haskey(br_res, "qf") && length(br_res["qf"]) >= 3
            qf = br_res["qf"]
            Qf_a[br][t] = qf[1]
            Qf_b[br][t] = qf[2]
            Qf_c[br][t] = qf[3]
        else
            Qf_a[br][t] = 0.0
            Qf_b[br][t] = 0.0
            Qf_c[br][t] = 0.0
        end
    end
end

# ---------------------------
# 3. Save CSVs (voltages, currents, reactive power) and print head
# ---------------------------
# Voltages CSV: magnitude and real/imag parts per phase, per bus
dfV = DataFrame(Hour = 1:n_hours)
for bk in all_bus_keys
    name = get(key_to_bus_name, bk, "bus_$bk")
    dfV[!, "$(name)_Va_r"] = real.(V_a[bk])
    dfV[!, "$(name)_Va_i"] = imag.(V_a[bk])
    dfV[!, "$(name)_Va_abs"] = abs.(V_a[bk])
    dfV[!, "$(name)_Vb_r"] = real.(V_b[bk])
    dfV[!, "$(name)_Vb_i"] = imag.(V_b[bk])
    dfV[!, "$(name)_Vb_abs"] = abs.(V_b[bk])
    dfV[!, "$(name)_Vc_r"] = real.(V_c[bk])
    dfV[!, "$(name)_Vc_i"] = imag.(V_c[bk])
    dfV[!, "$(name)_Vc_abs"] = abs.(V_c[bk])
end
CSV.write(joinpath(out_folder, "bus_voltages_phasors.csv"), dfV)
println("\nSaved bus_voltages_phasors.csv — head:")
display(first(dfV, 5))

# Branch currents CSV
dfI = DataFrame(Hour = 1:n_hours)
for br in all_branch_keys
    label = get(branch_key_to_label, br, "branch_$br")
    dfI[!, "$(label)_Ia_r"] = real.(I_a[br])
    dfI[!, "$(label)_Ia_i"] = imag.(I_a[br])
    dfI[!, "$(label)_Ia_abs"] = abs.(I_a[br])
    dfI[!, "$(label)_Ib_r"] = real.(I_b[br])
    dfI[!, "$(label)_Ib_i"] = imag.(I_b[br])
    dfI[!, "$(label)_Ib_abs"] = abs.(I_b[br])
    dfI[!, "$(label)_Ic_r"] = real.(I_c[br])
    dfI[!, "$(label)_Ic_i"] = imag.(I_c[br])
    dfI[!, "$(label)_Ic_abs"] = abs.(I_c[br])
end
CSV.write(joinpath(out_folder, "branch_currents_phasors.csv"), dfI)
println("\nSaved branch_currents_phasors.csv — head:")
display(first(dfI, 5))

# Reactive power CSV (per branch)
dfQ = DataFrame(Hour = 1:n_hours)
for br in all_branch_keys
    label = get(branch_key_to_label, br, "branch_$br")
    dfQ[!, "$(label)_Qf_a"] = Qf_a[br]
    dfQ[!, "$(label)_Qf_b"] = Qf_b[br]
    dfQ[!, "$(label)_Qf_c"] = Qf_c[br]
end
CSV.write(joinpath(out_folder, "branch_reactive_power.csv"), dfQ)
println("\nSaved branch_reactive_power.csv — head:")
display(first(dfQ, 5))

# Assume branch 2 key is known
branch2_key = "2"  # adjust if your branch key is different
branch2_label = get(branch_key_to_label, branch2_key, "branch_2")

# Create DataFrame with Hour column
dfQ_branch2 = DataFrame(Hour = 1:n_hours)

# Compute mean reactive power across phases for branch 2
dfQ_branch2[!, "$(branch2_label)_Q_mean"] = (Qf_a[branch2_key] .+ Qf_b[branch2_key] .+ Qf_c[branch2_key]) ./ 3

# Save CSV
CSV.write(joinpath(out_folder, "branch2_mean_reactive_power.csv"), dfQ_branch2)
println("\nSaved branch2_mean_reactive_power.csv — head:")
first(dfQ_branch2, 5)  # show first 5 rows


# Voltage plot (all buses, three phases)
pV = plot(title = "Phase-to-Neutral Voltages (All Buses)", xlabel = "Hour", ylabel = "Voltage (p.u.)", legend = :outerright)

# Line styles per bus
line_styles = Dict(1 => :solid, 2 => :dash, 3 => :dot)

# Phase labels
phases = ["a", "b", "c"]

for (i, bk) in enumerate(all_bus_keys)
    linestyle = get(line_styles, i, :dashdot)  # default style if more buses exist
    name = get(key_to_bus_name, bk, "bus_$bk")
    for ph in phases
        # Dynamic access to phase voltages
        V_arr = eval(Symbol("V_$ph"))[bk]  # V_a[bk], V_b[bk], V_c[bk]
        label = "$(name) Phase $(uppercase(ph))"
        plot!(1:n_hours, abs.(V_arr), label=label, lw=2, ls=linestyle)
    end
end

savefig(joinpath(out_folder, "voltages_all_buses.png"))
println("\nSaved voltages_all_buses.png")

# Branch current plot (all branches, three phases)
pI = plot(title = "Branch Currents (from-end) - All Branches", xlabel = "Hour", ylabel = "Current (p.u.)", legend = :outerright)

# Define line styles per branch
line_styles = Dict(1 => :solid, 2 => :dash)

# Phase labels
phases = ["a", "b", "c"]

for (i, br) in enumerate(all_branch_keys)
    linestyle = get(line_styles, i, :dot)  # default dotted if more branches exist
    for ph in phases
        # Dynamic access to phase currents
        I_arr = eval(Symbol("I_$ph"))[br]  # I_a[br], I_b[br], I_c[br]
        label = "$(get(branch_key_to_label, br, "branch_$br")) Phase $(uppercase(ph))"
        plot!(1:n_hours, abs.(I_arr), label=label, lw=2, ls=linestyle)
    end
end

savefig(joinpath(out_folder, "branch_currents_per_phase.png"))
println("Saved branch_currents_per_phase.png")

pQ = plot(title = "Branch Reactive Power per Phase", xlabel = "Hour", ylabel = "Reactive Power (p.u.)", legend = :outerright)

# Define line styles per branch
line_styles = Dict(1 => :solid, 2 => :dash)

# Phase labels
phases = ["a", "b", "c"]

for (i, br) in enumerate(all_branch_keys)
    linestyle = get(line_styles, i, :dot)  # default to dotted if more branches exist
    for ph in phases
        # Get reactive power array dynamically
        Q_arr = eval(Symbol("Qf_$ph"))[br]  # Qf_a[br], Qf_b[br], Qf_c[br]
        label = "$(get(branch_key_to_label, br, "branch_$br")) Phase $(uppercase(ph))"
        plot!(1:n_hours, Q_arr, label=label, lw=2, ls=linestyle)
    end
end

savefig(joinpath(out_folder, "branch_reactive_power_per_phase.png"))
println("Saved branch_reactive_power_per_phase.png")

# Reactive power plot (per branch, single series per branch as mean Qf)
pQ = plot(title="Branch Reactive Power (mean per-branch)", xlabel="Hour", ylabel="Reactive Power (p.u.)", legend=:outerright)
for br in all_branch_keys
    label = get(branch_key_to_label, br, "branch_$br")
    # use mean of phases or sum depending on your preference
    plot!(1:n_hours, (Qf_a[br] .+ Qf_b[br] .+ Qf_c[br]) ./ 3, label="$(label) Qf_mean", lw=2)
end
savefig(joinpath(out_folder, "branch_reactive_power_mean.png"))
println("Saved branch_reactive_power_mean.png")

# ---------------------------
# 5. Sequence components and Unbalance (VUF & CUF) — robust computation
# ---------------------------
# sequence operator matrix (Fortescue)
a = cis(2*pi/3)
A = (1/3) * [1 1 1; 1 a a^2; 1 a^2 a]

# choose target for sequence analysis:
target_bus_key = load_bus_key      # load bus (e.g. "1")
target_branch_key = supply_branch_key

println("\nSequence analysis: using bus ", target_bus_key, " (", key_to_bus_name[target_bus_key], ") and branch ", target_branch_key, " (", get(branch_key_to_label, target_branch_key, target_branch_key), ")")

V0 = zeros(ComplexF64, n_hours); V1 = zeros(ComplexF64, n_hours); V2 = zeros(ComplexF64, n_hours)
I0 = zeros(ComplexF64, n_hours); I1 = zeros(ComplexF64, n_hours); I2 = zeros(ComplexF64, n_hours)
VUF = zeros(Float64, n_hours); CUF = zeros(Float64, n_hours)

for t in 1:n_hours
    # phasor vector (Va, Vb, Vc) at target bus (phase-to-neutral, complex)
    vvec = [V_a[target_bus_key][t], V_b[target_bus_key][t], V_c[target_bus_key][t]]
    # remove mean (small DC offset) to reduce numerical artifacts
    vvec .-= mean(vvec)

    seqv = A * vvec
    V0[t] = seqv[1]; V1[t] = seqv[2]; V2[t] = seqv[3]

    # branch currents
    ivec = [I_a[target_branch_key][t], I_b[target_branch_key][t], I_c[target_branch_key][t]]
    ivec .-= mean(ivec)
    seqi = A * ivec
    I0[t] = seqi[1]; I1[t] = seqi[2]; I2[t] = seqi[3]

    # Compute unbalance factors as percentages
    # Guard against tiny denominators
    V0mag = abs(V0[t]); I0mag = abs(I0[t])
    V1mag = abs(V1[t]); I1mag = abs(I1[t])
    V2mag = abs(V2[t]); I2mag = abs(I2[t])

    VUF[t] = sqrt(V0mag^2 + V2mag^2) / (V1mag + eps) * 100
    CUF[t] = sqrt(I0mag^2 + I2mag^2) / (I1mag + eps) * 100
end

# Save sequence and unbalance CSV
dfSeq = DataFrame(Hour = 1:n_hours,
    V1_abs = abs.(V1), V2_abs = abs.(V2), V0_abs = abs.(V0),
    I1_abs = abs.(I1), I2_abs = abs.(I2), I0_abs = abs.(I0),
    VUF_pct = VUF, CUF_pct = CUF)
CSV.write(joinpath(out_folder, "sequence_unbalance.csv"), dfSeq)
println("\nSaved sequence_unbalance.csv — head:")
display(first(dfSeq, 6))

# Plots: sequence components
p_vs = plot(title="Voltage sequence components (target bus)", xlabel="Hour", ylabel="Magnitude (p.u.)")
plot!(1:n_hours, abs.(V1), label="V1 (pos)", lw=2)
plot!(1:n_hours, abs.(V2), label="V2 (neg)", lw=2)
plot!(1:n_hours, abs.(V0), label="V0 (zero)", lw=2)
savefig(joinpath(out_folder, "voltage_sequence_components.png"))
println("Saved voltage_sequence_components.png")

p_is = plot(title="Current sequence components (target branch)", xlabel="Hour", ylabel="Magnitude (p.u.)")
plot!(1:n_hours, abs.(I1), label="I1 (pos)", lw=2)
plot!(1:n_hours, abs.(I2), label="I2 (neg)", lw=2)
plot!(1:n_hours, abs.(I0), label="I0 (zero)", lw=2)
savefig(joinpath(out_folder, "current_sequence_components.png"))
println("Saved current_sequence_components.png")

# Plots: VUF and CUF
p_vuf = plot(1:n_hours, VUF, lw=2, xlabel="Hour", ylabel="VUF", title="Voltage Unbalance Factor (%)", legend=false, grid=true)
savefig(joinpath(out_folder, "VUF_plot.png"))
println("Saved VUF_plot.png")

p_cuf = plot(1:n_hours, CUF, lw=2, xlabel="Hour", ylabel="CUF", title="Current Unbalance Factor (%)", legend=false, grid=true)
savefig(joinpath(out_folder, "CUF_plot.png"))
println("Saved CUF_plot.png")

CSV.write("VUF_raw.csv", DataFrame(hour=1:n_hours, VUF=VUF))
CSV.write("CUF_raw.csv", DataFrame(hour=1:n_hours, CUF=CUF))
println("Raw VUF range: ", extrema(VUF))
println("Raw CUF range: ", extrema(CUF))
println("VUF: min=", minimum(VUF), " mean=", mean(VUF), " median=", median(VUF), " p90=", quantile(VUF, 0.9))
println("CUF: min=", minimum(CUF), " mean=", mean(CUF), " median=", median(CUF), " p90=", quantile(CUF, 0.9))

th = 0.20   # 20% in per-unit
high_idx = findall(x -> x > th, CUF)
println("Hours with CUF > 20%: ", high_idx)

plot(1:n_hours, CUF, lw=2, xlabel="Hour", ylabel="CUF (p.u.)", title="CUF (raw)")
scatter!(high_idx, CUF[high_idx], color=:red, label="CUF > 20%")
savefig("CUF_flagged.png")


println("\nAll done — CSVs and PNGs are in: ", abspath(out_folder))