using Pkg
Pkg.activate("./")
using rosetta_distribution_opf
import PowerModelsDistribution
import InfrastructureModels
using Ipopt
using JuMP
const PMD = PowerModelsDistribution
const RPMD = rosetta_distribution_opf
const IM = InfrastructureModels


# data_path = "./data/test_gen_3ph_wye.dss"

##
using Plots
f = 50
phi = [0, -20 , 20] * pi/180 # rad
Vm = [1, 1, 1]
t = 0:.001:.03

function vabc(f, phi, Vm, t)
    w = 2*pi*f
    vat = Vm[1] * cos.(w * t .+ phi[1])
    vbt = Vm[2] * cos.(w * t .+ phi[2] .- 2*pi/3)
    vct = Vm[3] * cos.(w * t .+ phi[3] .+ 2*pi/3)
    vabct = [vat' ; vbt' ; vct']
    return vabct
end

function clarke_transform(vabct)
    Tc = 2/3 * [1 -1/2 -1/2 ; 0 sqrt(3)/2 -sqrt(3)/2 ; 1/2 1/2 1/2]
    vabgt = zeros(size(vabct))
    for i in 1:size(vabct,2)
        vabgt[:,i] = Tc * vabct[:,i]
    end
    return vabgt
end


function plot_abc_abg(vabct, vabgt)
    vabc_plot = plot(t, vabct[1,:], label="va")
    plot!(t, vabct[2,:], label="vb")
    plot!(t, vabct[3,:], label="vc")

    vabg_plot = plot(t, vabgt[1,:], label="v_alpha")
    plot!(t, vabgt[2,:], label="v_beta")
    plot!(t, vabgt[3,:], label="v_gamma")

    plot(vabc_plot, vabg_plot, layout=(2,1))
end

vabct = vabc(f, phi, Vm, t)
vabgt = clarke_transform(vabct)
plot_abc_abg(vabct, vabgt)


