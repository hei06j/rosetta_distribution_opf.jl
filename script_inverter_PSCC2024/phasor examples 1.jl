using Unitful, Unitful.DefaultSymbols, PyPlot, ElectricalEngineering

cd("Figures")

## Currents
# balance
Ia = cis(0)
Ib = cis(deg2rad(-120))
Ic = cis(deg2rad(120))
In = -(Ia+Ib+Ic)
a = 2.5 # plot scale

rc("text", usetex=true); rc("font", family="sans-serif", size=16)
phasorcosine(abs(Ia),angle(Ia), ylabel=L"$i$", maglabel=L"$\hat{I}_a$", 
    labelrsep=0.5,
    figsize=(7*a,2.5*a),
    color="blue", linestyle="-", 
    )
phasorcosine(abs(Ib),angle(Ib), ylabel=L"$u$", 
    labelrsep=0.5,
    color="green", linestyle="-", add=true)
phasorcosine(abs(Ic), angle(Ic), ylabel=L"$i$", 
    labelrsep=0.5,
    color="red", linestyle="-", add=true)
phasorcosine(abs(In), angle(In), ylabel=L"$i$", # maglabel=L"$\hat{I}_n$", 
    labelrsep=0.5, color="black", linestyle="--", add=true)
gcf()
save2fig("Balance", dpi=300, crop=true);

# Inbalance
Ia = cis(deg2rad(-25))
Ib = 0.6*cis(deg2rad(-130))
Ic = 0.75*cis(deg2rad(125))
In = -(Ia+Ib+Ic)

rc("text", usetex=true); rc("font", family="sans-serif", size=16)
phasorcosine(abs(Ia),angle(Ia), ylabel=L"$i$", maglabel=L"$\hat{I}_a$", 
    labelrsep=0.5,
    figsize=(7*a,2.5*a),
    color="blue", linestyle="-")
phasorcosine(abs(Ib),angle(Ib), ylabel=L"$i$", maglabel=L"$\hat{I}_b$", 
    labelrsep=0.5,
    color="green", linestyle="-", add=true)
phasorcosine(abs(Ic), angle(Ic), ylabel=L"$u$", maglabel=L"$\hat{I}_c$", 
    labelrsep=0.5,labeltsep=-0.1, labelrelrot=true, labelrelangle=deg2rad(180),
    color="red", linestyle="-", add=true)
phasorcosine(abs(In), angle(In), ylabel=L"$u$", maglabel=L"$\hat{I}_n$", 
    labelrsep=0.5,labeltsep=-0.1, labelrelrot=true, labelrelangle=deg2rad(180),
    color="black", linestyle="--", add=true)
gcf()
save2fig("Unbalance", dpi=300, crop=true);

# Inbalance reverse power flow
Ia = cis(deg2rad(-25))
Ib = 0.15*cis(deg2rad(65))
Ic = 0.75*cis(deg2rad(125))
In = -(Ia+Ib+Ic)

rc("text", usetex=true); rc("font", family="sans-serif", size=16)
phasorcosine(abs(Ia),angle(Ia), ylabel=L"$i$", maglabel=L"$\hat{I}_a$", 
    labelrsep=0.5,
    figsize=(7*a,2.5*a),
    color="blue", linestyle="-")
phasorcosine(abs(Ib),angle(Ib), ylabel=L"$i$", maglabel=L"$\hat{I}_b$", 
    labelrsep=0.5, labeltsep=-0.1, labelrelrot=true, labelrelangle=deg2rad(180), 
    color="green", linestyle="-", add=true)
phasorcosine(abs(Ic), angle(Ic), ylabel=L"$u$", maglabel=L"$\hat{I}_c$", 
    labelrsep=0.5, labeltsep=-0.1, labelrelrot=true, labelrelangle=deg2rad(180), 
    color="red", linestyle="-", add=true)
phasorcosine(abs(In), angle(In), ylabel=L"$u$", maglabel=L"$\hat{I}_n$", 
    labelrsep=0.5,
    color="black", linestyle="--", add=true)
gcf()
save2fig("UnbalanceReverse", dpi=300, crop=true);

    # phasorcosine(0.55, 20, add=true, maglabel=L"$\hat{I}_a$")
# phasorcosine(0.55, -100, add=true, maglabel=L"$\hat{I}_b$")
# phasorcosine(0.55, 140, add=true, maglabel=L"$\hat{I}_c$")

# phasorcosine(1, 45°, ylabel=L"$u$", maglabel=L"$\hat{U}_a$", labelrsep=0.3,
#     color="blue", linestyle="-",figsize=(15,7))
# phasorcosine(0.6, -135°, ylabel=L"$u$", maglabel=L"$\hat{U}_b$", labelrsep=0.3,
#     color="green", linestyle="-", add=true)
# phasorcosine(0.8, 125°, ylabel=L"$u$", maglabel=L"$\hat{U}_c$", labelrsep=0.3,
#     color="red", linestyle="-", add=true)
# phasorcosine(0.95, 110°, ylabel=L"$u$", maglabel=L"$\hat{U}_n$", labelrsep=0.3,
#     color="black", linestyle="-", add=true)



# phasorsine(1, 45°, ylabel=L"$u,i$", maglabel=L"$\hat{U}$", labelrsep=0.3,
#     color="gray", linestyle="--")
# phasorsine(0.55, 0, add=true, maglabel=L"$\hat{I}$")
# save2fig("phasorsine",crop=true);

## Voltage phasors

alpha = exp(im*2/3*pi)
T = 1/3 * [1 1 1 ; 1 alpha alpha^2 ; 1 alpha^2 alpha]
Tre = real.(T)
Tim = imag.(T)


# balance
Va = cis(0)
Vb = cis(deg2rad(120))
Vc = cis(deg2rad(-120))
a = 2.5 # plot scale

rc("text", usetex=true); rc("font", family="sans-serif", size=16)
phasorcosine(abs(Va),angle(Va), ylabel=L"$v$", maglabel=L"$\hat{V}_a$", 
    labelrsep=0.5,
    figsize=(7*a,2.5*a),
    color="blue", linestyle="-", 
    )
phasorcosine(abs(Vb),angle(Vb), ylabel=L"$v$", 
    labelrsep=0.5,
    color="green", linestyle="-", add=true)
phasorcosine(abs(Vc), angle(Vc), ylabel=L"$v$", 
    labelrsep=0.5,
    color="red", linestyle="-", add=true)
# phasorcosine(abs(In), angle(In), ylabel=L"$v$", # maglabel=L"$\hat{V}_n$", 
#     labelrsep=0.5, color="black", linestyle="--", add=true)
gcf()
ElectricalEngineering.save2fig("Balanced_voltage", dpi=300, crop=true);


V012 = T * [Va Vb Vc]'
round.(V012, digits=4)
round.(abs.(V012), digits=4)

rc("text", usetex=true); rc("font", family="sans-serif", size=16)
phasorcosine(abs(V012[2]),angle(V012[2]), ylabel=L"$v$", maglabel=L"$\hat{V}_1$", 
    labelrsep=0.5,
    figsize=(7*a,2.5*a),
    color="blue", linestyle="-", 
    )
phasorcosine(abs(V012[3]),angle(V012[3]), ylabel=L"$v$", maglabel=L"$\hat{V}_2$",
    labelrsep=0.5,
    color="green", linestyle="-", add=true)
phasorcosine(abs(V012[1]), angle(V012[1]), ylabel=L"$v$", maglabel=L"$\hat{V}_0$",
    labelrsep=0.5,
    color="red", linestyle="-", add=true)
gcf()
ElectricalEngineering.save2fig("Balanced_voltage_sequence", dpi=300, crop=true);


# Unbalance
Va = cis(deg2rad(-25))
Vb = 0.6*cis(deg2rad(140))
Vc = 0.75*cis(deg2rad(-105))

rc("text", usetex=true); rc("font", family="sans-serif", size=16)
phasorcosine(abs(Va),angle(Va), ylabel=L"$v$", maglabel=L"$\hat{V}_a$", 
    labelrsep=0.5,
    figsize=(7*a,2.5*a),
    color="blue", linestyle="-")
phasorcosine(abs(Vb),angle(Vb), ylabel=L"$v$", maglabel=L"$\hat{V}_b$", 
    labelrsep=0.5,
    color="green", linestyle="-", add=true)
phasorcosine(abs(Vc), angle(Vc), ylabel=L"$u$", maglabel=L"$\hat{V}_c$", 
    labelrsep=0.5,labeltsep=-0.1, labelrelrot=true, labelrelangle=deg2rad(180),
    color="red", linestyle="-", add=true)
gcf()
save2fig("Unbalanced_voltage", dpi=300, crop=true);


V012 = T * [Va Vb Vc]'
round.(V012, digits=4)
round.(abs.(V012), digits=4)


rc("text", usetex=true); rc("font", family="sans-serif", size=16)
phasorcosine1(abs(V012[2]),angle(V012[2]), ylabel=L"$v$", maglabel=L"$\hat{V}_1$", 
    labelrsep=0.5,
    figsize=(7*a,2.5*a),
    color="blue", linestyle="-", 
    )
phasorcosine1(abs(V012[3]),angle(V012[3]), ylabel=L"$v$", maglabel=L"$\hat{V}_2$",
    labelrsep=0.5,
    color="green", linestyle="-", add=true)
phasorcosine1(abs(V012[1]), angle(V012[1]), ylabel=L"$v$", maglabel=L"$\hat{V}_0$",
    labelrsep=0.5,
    color="red", linestyle="-", add=true)
gcf()
ElectricalEngineering.save2fig("Unbalanced_voltage_sequence", dpi=300, crop=true);

# # Inbalance reverse power flow
function phasorcosine1(mag = 1,
    phi = 0;
    add = false,
    figsize = (3.3,1.5),
    xlabel = L"$\omega\!\cdot\!t/^\circ $",
    ylabel = "",
    maglabel = "",
    phasorlabel = maglabel,
    color = "black",
    backgroundcolor = "none",
    linewidth = 1,
    linestyle = "-",
    labeltsep = 0.1,
    labelrsep = 0.5,
    labelrelrot = true,
    labelrelangle = 0,
    colorDash="gray",
    left=0.20,
    right=0.80,
    bottom=0.20,
    top=0.80,
    showcosine = true,
    showdashline = true,
    shift = true,
    marker = "")
    # https://matplotlib.org/tutorials/text/annotations.html#plotting-guide-annotation
    # https://matplotlib.org/users/annotations.html
    # https://stackoverflow.com/questions/17543359/drawing-lines-between-two-plots-in-matplotlib

    # Create figure
    if !add
        # Create new figure
        fig = figure(figsize=figsize)
        subplots_adjust(left=left, right=right, bottom=bottom, top=top)
    end
    # Create left subplot
    # subplot(121)
    # Angle vector to draw circle
    psi = collect(0:pi/500:2*pi)
    # Coordinates of circle
    x = mag*cos.(psi)
    y = mag*sin.(psi)
    # Plot circle
    plot(x, y, color=colorDash, linewidth=1, linestyle=":",
        dash_capstyle="round")
    # Plot phasor
    phasor(pol(mag,phi.+pi/2), ref=1,
        label=phasorlabel, labelrsep=labelrsep, labeltsep=labeltsep,
        labelrelrot=labelrelrot, labelrelangle=labelrelangle,
        color=color, backgroundcolor=backgroundcolor,
        linestyle=linestyle, linewidth=linewidth)
    axis("square")
    ax1 = gca()
    xlim(-1.5,1.1)
    ylim(-1.1,1.43)
    if !add
        ax1.spines["left"].set_visible(false)
        ax1.spines["right"].set_visible(false)
        ax1.spines["bottom"].set_visible(false)
        ax1.spines["top"].set_visible(false)
        # Remove ticks
        xticks([])
        yticks([])
    end
end
