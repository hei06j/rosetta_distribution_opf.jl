function plot_phasors(phasor, Imax; labeled=false, I2=[], I0=[])
    plt = Plots.plot([0,imag.(phasor[1])], [0,real.(phasor[1])], arrow=true, color=:blue, linewidth=3, linestyle=:solid, label="a", border=:none)
    Plots.plot!([0,imag.(phasor[2])], [0,real.(phasor[2])], arrow=true, color=:red, linewidth=3, linestyle=:solid, label="b", border=:none)
    Plots.plot!([0,imag.(phasor[3])], [0,real.(phasor[3])], arrow=true, color=:green, linewidth=3, linestyle=:solid, label="c", border=:none)
    if phasor[4] !==  0 + 0im
        Plots.plot!([0,imag.(phasor[4])], [0,real.(phasor[4])], arrow=true, color=:black, linewidth=3, linestyle=:solid, label="n", border=:none)
    end
    Plots.plot!([0,0], [0,1.1*Imax], arrow=true, color=:grey, linestyle=:dot, label=false)
    Plots.plot!([0,1.1*Imax*real(exp(im*210/180*pi))], [0,1.1*Imax*imag(exp(im*210/180*pi))], arrow=true, color=:grey, linestyle=:dot, label=false)
    Plots.plot!([0,1.1*Imax*real(exp(im*330/180*pi))], [0,1.1*Imax*imag(exp(im*330/180*pi))], arrow=true, color=:grey, linestyle=:dot, label=false)
    if labeled
        Plots.plot!(Imax*exp.(im*(0:0.01:2pi)), color=:black, border=:none, label=false, markersize=10, legend=:bottom, legendcolumns=4, legendfontsize=30)
    else
        Plots.plot!(Imax*exp.(im*(0:0.01:2pi)), color=:black, border=:none, label=false, markersize=10, legend=false)
    end
    if !isempty(I2)
        Plots.annotate!([-7], [-Imax], Plots.text(latexstring("I_2= $(I2)"), :black, 40))
    end
    if !isempty(I0)
        Plots.annotate!([-7], [-Imax+4], Plots.text(latexstring("I_0= $(I0)"), :black, 40))
    end
    return plt
end
