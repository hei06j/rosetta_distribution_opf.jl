using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP  # bl/array_nl
import LinearAlgebra: diag, diagm
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels


##
include("./ENWL_OPF_NoInv_4w.jl")
include("./ENWL_OPF_GFL_4w.jl")
include("./ENWL_OPF_GFM_4w.jl")

##
using Plots

histogram(c2m.*100, label=false, xlabel="Current Negative Sequence (%)")

v2m_GFLs_plot = histogram(v2m_GEN.*100, alpha=0.5, bins=range(0,1.5, step = 0.1), label=false)
histogram!(v2m_GFL.*100, alpha=0.5, bins=range(0,1.5, step = 0.1), label=false)
histogram!(v2m_GFM.*100, alpha=0.5, bins=range(0,1.5, step = 0.1), label=false, xlabel="Voltage Negative Sequence (%)")
# vline!([2], label="Grid Code Limit", linewidth=2, legend=:topright)
savefig(v2m_GFLs_plot, "./Figures/v2m_GFM_plot.pdf")

v20_GFLs_plot = histogram(v0m_GEN.*100, alpha=0.5, bins=range(0,1.2, step = 0.1), label="No inverter")
histogram!(v0m_GFL.*100, alpha=0.5, bins=range(0,1.2, step = 0.1), label="GFLs 4-w")
histogram!(v0m_GFM.*100, alpha=0.5, bins=range(0,1.2, step = 0.1), label="GFM/GFL 3/4-w", xlabel="Voltage Zero Sequence (%)")
savefig(v20_GFLs_plot, "./Figures/v20_GFLs_plot.pdf")

v0_v2_plot = plot(v2m_GFLs_plot, v20_GFLs_plot, layout=(2,1), size=(600, 300))
savefig(v0_v2_plot, "./Figures/v0_v2_plot.pdf")
