using Plots
using LaTeXStrings

w = 2*pi*50
t = 0:0.001:0.051
l = Int(length(t))

Ima = [1*ones(Int(l/2)) ;  1*ones(Int(l/2))]
Imb = [1*ones(Int(l/2)) ;  0.8*ones(Int(l/2))]
Imc = [1*ones(Int(l/2)) ;  0.6*ones(Int(l/2))]

Iaa = pi/180 * [   0*ones(Int(l/2)) ;  40*ones(Int(l/2))]
Iab = pi/180 * [ 120*ones(Int(l/2)) ; 130*ones(Int(l/2))]
Iac = pi/180 * [-120*ones(Int(l/2)) ;-140*ones(Int(l/2))]


Ia = Ima .* sin.(w*t .+ Iaa)
Ib = Imb .* sin.(w*t .+ Iab)
Ic = Imc .* sin.(w*t .+ Iac)
In = Ia .+ Ib .+ Ic
I = [Ia Ib Ic In]'

unbalance_ac_signals = plot(t, Ia, color=:blue, linewidth=3, label="phase a", ylabel="Current (A pu)", size=(800,400))
plot!(t, Ib, label="phase b", color=:red, linewidth=3, xticks=(0:0.01:0.05, 0:10:50), xlabel="Time (ms)")
plot!(t, Ic, label="phase c", color=:green, linewidth=3)
plot!(t, In, label="neutral", color=:black, linewidth=3)
vline!([t[Int(l/2)]], color=:black, label=false)
Plots.savefig(unbalance_ac_signals, "./Figures/unbalance_ac_signals.png")

##
alpha = exp(im*2/3*pi)
T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
Tre = real.(T)
Tim = imag.(T)

Ima = 1;  Imb = 1; Imc = 1;
Iaa = 0*pi/180;  Iab = 120*pi/180; Iac = -120*pi/180;
Ia = Ima .* exp.(im*(w*t .+ Iaa))
Ib = Imb .* exp.(im*(w*t .+ Iab))
Ic = Imc .* exp.(im*(w*t .+ Iac))
In = Ia .+ Ib .+ Ic
I = [Ia Ib Ic In]'
I012 = T * [Ia Ib Ic]'
I0, I1, I2 = abs.(I012[:,1])
plt = plot_phasors(I, labeled=true, I1=I1, I2=I2, I0=I0)
Plots.savefig(plt, "./Figures/unbalance_phasors1.png")


Ima = 1;  Imb = 0.8; Imc = 0.6;
Iaa = 40*pi/180;  Iab = 130*pi/180; Iac = -140*pi/180;
Ia = Ima .* exp.(im*(w*t .+ Iaa))
Ib = Imb .* exp.(im*(w*t .+ Iab))
Ic = Imc .* exp.(im*(w*t .+ Iac))
In = Ia .+ Ib .+ Ic
I = [Ia Ib Ic In]'
I012 = T * [Ia Ib Ic]'
I0, I1, I2 = abs.(I012[:,1])
plt = plot_phasors(I, labeled=true, I1=I1, I2=I2, I0=I0)
Plots.savefig(plt, "./Figures/unbalance_phasors2.png")



function plot_phasors(phasor; Imax=1, labeled=false, I1=[], I2=[], I0=[])
    plt = Plots.plot([0,imag.(phasor[1])], [0,real.(phasor[1])], arrow=true, color=:blue, linewidth=3, linestyle=:solid, label="a", border=:none)
    Plots.plot!([0,imag.(phasor[2])], [0,real.(phasor[2])], arrow=true, color=:red, linewidth=3, linestyle=:solid, label="b", border=:none)
    Plots.plot!([0,imag.(phasor[3])], [0,real.(phasor[3])], arrow=true, color=:green, linewidth=3, linestyle=:solid, label="c", border=:none)
    Plots.plot!([0,imag.(phasor[4])], [0,real.(phasor[4])], arrow=true, color=:black, linewidth=4, linestyle=:solid, label="n", border=:none)
    Plots.plot!([0,0], [0,1.1*Imax], arrow=true, color=:grey, linestyle=:dot, label=false)
    Plots.plot!([0,1.1*Imax*real(exp(im*210/180*pi))], [0,1.1*Imax*imag(exp(im*210/180*pi))], arrow=true, color=:grey, linestyle=:dot, label=false)
    Plots.plot!([0,1.1*Imax*real(exp(im*330/180*pi))], [0,1.1*Imax*imag(exp(im*330/180*pi))], arrow=true, color=:grey, linestyle=:dot, label=false)
    if labeled
        Plots.plot!(Imax*exp.(im*(0:0.01:2pi)), color=:black, border=:none, label=false, markersize=5, legend=:topright, legendcolumns=4, legendfontsize=10)
    else
        Plots.plot!(Imax*exp.(im*(0:0.01:2pi)), color=:black, border=:none, label=false, markersize=5, legend=false)
    end
    annotate!([-0.8], [-Imax+0.2], text(latexstring("Pos Seq= $(round(I1, digits=2))"), :black, 20))
    annotate!([-0.8], [-Imax], text(latexstring("Neg Seq= $(round(I2, digits=2))"), :black, 20))
    annotate!([-0.8], [-Imax-0.2], text(latexstring("Zero Seq= $(round(I0, digits=2))"), :black, 20))
    return plt
end
plt = plot_phasors(I, labeled=true, I1=I1, I2=I2, I0=I0)
