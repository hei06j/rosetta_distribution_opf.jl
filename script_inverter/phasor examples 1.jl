using Unitful, Unitful.DefaultSymbols, PyPlot, ElectricalEngineering

cd("Figures")
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