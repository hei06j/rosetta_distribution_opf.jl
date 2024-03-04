using Plots
using LaTeXStrings

w = 2*pi*50
t = 0:0.001:0.103
l = Int(length(t))
theta = w * t *180/pi

Ima = [1*ones(Int(l/4)) ;  1*ones(Int(l/4)) ;  1*ones(Int(l/4)) ;  1*ones(Int(l/4))]
Imb = [1*ones(Int(l/4)) ;  1*ones(Int(l/4)) ;  0.9*ones(Int(l/4)) ;  1*ones(Int(l/4))]
Imc = [1*ones(Int(l/4)) ;  1*ones(Int(l/4)) ;  0.8*ones(Int(l/4)) ;  1*ones(Int(l/4))]

Iaa = pi/180 * [   0*ones(Int(l/4)) ;  20*ones(Int(l/4)) ;   0*ones(Int(l/4)) ;  20*ones(Int(l/4))]
Iab = pi/180 * [ 120*ones(Int(l/4)) ; 130*ones(Int(l/4)) ; 120*ones(Int(l/4)) ; 130*ones(Int(l/4))]
Iac = pi/180 * [-120*ones(Int(l/4)) ;-140*ones(Int(l/4)) ;-120*ones(Int(l/4)) ;-140*ones(Int(l/4))]


Ia = Ima .* sin.(w*t .+ Iaa)
Ib = Imb .* sin.(w*t .+ Iab)
Ic = Imc .* sin.(w*t .+ Iac)
In = Ia .+ Ib .+ Ic
Iabc = [Ia Ib Ic]'
I = [Ia Ib Ic In]'

unbalance_abc_signals = plot(t, Ia, color=:blue, linewidth=3, label="phase a", ylabel="Current (A pu)", size=(1200,400))
plot!(t, Ib, label="phase b", color=:red, linewidth=3, xticks=(0:0.01:0.100, 0:10:100), xlabel="Time (ms)")
plot!(t, Ic, label="phase c", color=:green, linewidth=3)
plot!(t, In, label="neutral", color=:black, linewidth=3)
vline!([t[Int(l/4)]], color=:black, label=false)
vline!([t[2*Int(l/4)]], color=:black, label=false)
vline!([t[3*Int(l/4)]], color=:black, label=false)
Plots.savefig(unbalance_abc_signals, "./Figures/unbalance_abc_signals.png")

Tab0 = 2/3 * [1 -1/2 -1/2 ; 0 sqrt(3)/2 -sqrt(3)/2 ; 1/2 1/2 1/2]
Iab0 = Tab0 * Iabc
unbalance_ab0_signals = plot(t, Iab0[1,:], color=:blue, linewidth=3, label=L"\alpha", ylabel="Current (A pu)", size=(800,400))
plot!(t, Iab0[2,:], label=L"\beta", color=:red, linewidth=3, xticks=(0:0.01:0.100, 0:10:100), xlabel="Time (ms)")
plot!(t, Iab0[3,:], label=L"0", color=:green, linewidth=3)
# plot!(t, In, label="neutral", color=:black, linewidth=3)
vline!([t[Int(l/4)]], color=:black, label=false)
vline!([t[2*Int(l/4)]], color=:black, label=false)
vline!([t[3*Int(l/4)]], color=:black, label=false)
Plots.savefig(unbalance_ab0_signals, "./Figures/unbalance_ab0_signals.png")


# theta = pi/2
# Tdq00 = sqrt(2/3) * [cos(theta)   cos(theta-2*pi/3)  cos(theta+2*pi/3) ; -sin(theta) -sin(theta-2*pi/3) -sin(theta+2*pi/3) ; sqrt(1/2) sqrt(1/2) sqrt(1/2)]
# Idq0 = Tdq00 * Iabc

Tdq0(theta) = sqrt(2/3) * [cos(theta)   cos(theta-2*pi/3)  cos(theta+2*pi/3) ; -sin(theta) -sin(theta-2*pi/3) -sin(theta+2*pi/3) ; sqrt(1/2) sqrt(1/2) sqrt(1/2)]
I_abc2dq0(theta, Iabc) = sqrt(2/3) * [cos(theta)   cos(theta-2*pi/3)  cos(theta+2*pi/3) ; -sin(theta) -sin(theta-2*pi/3) -sin(theta+2*pi/3) ; sqrt(1/2) sqrt(1/2) sqrt(1/2)] * Iabc
I_ab02dq0(theta, Iab0) = sqrt(2/3) * [cos(theta)   sin(theta)  0 ; -sin(theta) cos(theta) 0 ; 0 0 1] * Iab0

Idq0 = []
for i in 1:length(t)
    # push!(Idq0, Tdq0(w*t[i])*Iabc[:,i])
    # push!(Idq0, I_abc2dq0(w*t[i], Iabc[:,i]))
    push!(Idq0, I_ab02dq0(w*t[i], Iab0[:,i]))
end
Idq0 = permutedims(hcat(Idq0...))'

##
unbalance_dq0_signals = plot(t, Idq0[1,:], color=:blue, linewidth=3, label=L"d", ylabel="Current (A pu)", size=(800,400))
plot!(t, Idq0[2,:], label=L"d", color=:red, linewidth=3, xticks=(0:0.01:0.100, 0:10:100), xlabel="Time (ms)")
plot!(t, Idq0[3,:], label=L"0", color=:green, linewidth=3)
# plot!(t, In, label="neutral", color=:black, linewidth=3)
vline!([t[Int(l/4)]], color=:black, label=false)
vline!([t[2*Int(l/4)]], color=:black, label=false)
vline!([t[3*Int(l/4)]], color=:black, label=false)
Plots.savefig(unbalance_dq0_signals, "./Figures/unbalance_dq0_signals.png")

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

Ima = 1;  Imb = 1; Imc = 1;
Iaa = 20*pi/180;  Iab = 130*pi/180; Iac = -140*pi/180;
Ia = Ima .* exp.(im*(w*t .+ Iaa))
Ib = Imb .* exp.(im*(w*t .+ Iab))
Ic = Imc .* exp.(im*(w*t .+ Iac))
In = Ia .+ Ib .+ Ic
I = [Ia Ib Ic In]'
I012 = T * [Ia Ib Ic]'
I0, I1, I2 = abs.(I012[:,1])
plt = plot_phasors(I, labeled=true, I1=I1, I2=I2, I0=I0)
Plots.savefig(plt, "./Figures/unbalance_phasors2.png")

Ima = 1;  Imb = 0.9; Imc = 0.8;
Iaa = 0*pi/180;  Iab = 120*pi/180; Iac = -120*pi/180;
Ia = Ima .* exp.(im*(w*t .+ Iaa))
Ib = Imb .* exp.(im*(w*t .+ Iab))
Ic = Imc .* exp.(im*(w*t .+ Iac))
In = Ia .+ Ib .+ Ic
I = [Ia Ib Ic In]'
I012 = T * [Ia Ib Ic]'
I0, I1, I2 = abs.(I012[:,1])
plt = plot_phasors(I, labeled=true, I1=I1, I2=I2, I0=I0)
Plots.savefig(plt, "./Figures/unbalance_phasors3.png")

Ima = 1;  Imb = 0.9; Imc = 0.8;
Iaa = 20*pi/180;  Iab = 130*pi/180; Iac = -140*pi/180;
Ia = Ima .* exp.(im*(w*t .+ Iaa))
Ib = Imb .* exp.(im*(w*t .+ Iab))
Ic = Imc .* exp.(im*(w*t .+ Iac))
In = Ia .+ Ib .+ Ic
I = [Ia Ib Ic In]'
I012 = T * [Ia Ib Ic]'
I0, I1, I2 = abs.(I012[:,1])
plt = plot_phasors(I, labeled=true, I1=I1, I2=I2, I0=I0)
Plots.savefig(plt, "./Figures/unbalance_phasors4.png")



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
