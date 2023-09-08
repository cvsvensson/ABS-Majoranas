using DrWatson
@quickactivate "Majorana sweet spot"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))
using DataFrames
resultsUVp = collect_results(datadir("UV-scan", "parallel"))
resultsUVap = collect_results(datadir("UV-scan", "anti_parallel"))
resultsUhp = collect_results(datadir("Uh-scan", "parallel"))
resultsUhap = collect_results(datadir("Uh-scan", "anti_parallel"))


## UV
UVpfig, _, _ = plot_sweet_scan(resultsUVp[1, :])
display(UVpfig)
UVapfig, _, _ = plot_sweet_scan(resultsUVap[1, :])
display(UVapfig)
## Uh
Uhpfig, _, _ = plot_sweet_scan(resultsUhp[1, :])
display(Uhpfig)
Uhapfig, _, _ = plot_sweet_scan(resultsUhap[1, :])
display(Uhapfig)
##

fig, _, _ = plot_sweet_scan(resultsUhap[1, :]; datamap=x -> x.parameters.μ1 - x.parameters.U - x.parameters.h, colorscale=identity, colorrange=(-2, 2), colorbar=false)
display(fig)
fig, _, _ = plot_sweet_scan(resultsUhap[1, :]; datamap=x -> x.parameters.μ2 + x.parameters.h, colorscale=identity, colorrange=(-2, 2), colorbar=false)
display(fig)
fig, _, _ = plot_sweet_scan(resultsUhap[1, :]; datamap=x -> x.parameters.ϕ, colorscale=identity, colorrange=(0, pi), colorbar=false)
display(fig)

##
scanfig = let data1 = resultsUhap[3, :], data2 = resultsUVap[3, :]
    fig = Figure(; resolution=400 .* (2, 1), fontsize=20, backgroundcolor=:transparent)
    xlabel1 = paramstyle[data1[:xlabel]]
    xlabel2 = paramstyle[data2[:xlabel]]
    ylabel1 = paramstyle[data1[:ylabel]]
    ylabel2 = paramstyle[data2[:ylabel]]
    ax1 = Axis(fig[1, 1]; xlabel=xlabel1, ylabel=ylabel1)
    ax2 = Axis(fig[1, 2]; xlabel=xlabel2, ylabel=ylabel2)
    g = fig[1, 1:3] = GridLayout()
    hm1 = plot_sweet_scan!(ax1, data1)
    hm2 = plot_sweet_scan!(ax2, data2)
    label = L"1-\text{MP}"
    add_exp_colorbar!(fig[1, 3], hm1;)
    paramstring1 = map(x -> round(x, digits=2), data1[:fixedparams])
    paramstring2 = map(x -> round(x, digits=2), data2[:fixedparams])

    supertitle = Label(fig[1, 3, Top()], label, fontsize=25, padding=(0, 0, 8, 0))
    colgap!(fig.layout, 1, 16)
    colgap!(fig.layout, 2, 10)
    #save(plotsdir(string("uhvplot_transparent", paramstring1, "_", paramstring2, ".png")), fig, px_per_unit=4)
    fig
end


##
ssdata = resultsUhap[1, :]
nx, ny = size(ssdata[:sweet_spots])
smallres = 200
# positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.55, 0.75), (0.25, 0.5), (0.22, 0.1), (0.25, 0.1)])
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.5, 0.8), (0.3, 0.4), (0.1, 0.1)])
sweet_spots = [ssdata[:sweet_spots][pos...] for pos in positions]
paramstring = map(ss -> map(x -> round(x, digits=2), ss.parameters), sweet_spots)
transport = Transport(QuantumDots.Pauli(), (; T=1 / 20, μ=(0.0, 0.0)))
csdata_small = [charge_stability_scan((; ss.parameters...), 0.8, 0.8, smallres; transport) for ss in sweet_spots]

##
colors = cgrad(:rainbow, categorical=true)[[1, 2, 5, 4][1:length(csdata_small)]]
fixedparamstring = map(x -> round(x, digits=2), ssdata[:fixedparams])
f = Figure(resolution=400 .* (1.4, 1), fontsize=20, backgroundcolor=:white);
f = Figure(resolution=400 .* (1.4, 1), fontsize=20, backgroundcolor=:transparent);
g = f[1, 1] = GridLayout()
gb = g[1, 1] = GridLayout()
gs = g[1, 2] = GridLayout()

xlabel = paramstyle[ssdata[:xlabel]]
ylabel = paramstyle[ssdata[:ylabel]]
ax = Axis(gb[1, 1]; xlabel, ylabel)
hm = plot_sweet_scan!(ax, ssdata)

cb = add_exp_colorbar!(gb[2, 1], hm; vertical=false, label=L"1-\mathrm{MP}", flipaxis=false, labelpadding=-5, height=12)

scatter!(ax, map(ss -> ss.parameters[ssdata.xlabel], sweet_spots), map(ss -> ss.parameters[ssdata.ylabel], sweet_spots),
    color=colors, markersize=20, marker=:xcross, strokewidth=2)

datamap = x -> sign(x.gap)
datamap2 = x -> real(x.conductance)[1, 2] * real(x.conductance)[1, 1]
datamap2 = x -> real(x.conductance)[1, 2] - real(x.conductance)[2, 1]
datamap2 = x -> real(x.conductance)[1, 2]
spinecolor = [NamedTuple(map(x -> x => colors[n], [:bottomspinecolor, :leftspinecolor, :topspinecolor, :rightspinecolor])) for n in eachindex(csdata_small)]
small_axes = [Axis(gs[n, 1]; spinecolor[n]..., spinewidth=4, aspect=1) for (n, data) in enumerate(csdata_small)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
small_axes2 = [Axis(gs[n, 2]; spinecolor[n]..., spinewidth=4, aspect=1) for (n, data) in enumerate(csdata_small)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
foreach((ax, data) -> plot_charge_stability!(ax, data; datamap), small_axes, csdata_small)
foreach((ax, data) -> plot_charge_stability!(ax, data; datamap=datamap2, colormap=:vik), small_axes2, csdata_small)
hidedecorations!.(small_axes)
hidedecorations!.(small_axes2)
labels = [Label(gs[n, 1:2, Bottom()], L"MP ≈ %$(round(1-MPU(ss),digits=2))", padding=(0, 0, -15, 4)) for (n, ss) in enumerate(sweet_spots)]
Label(gs[1, 1, Top()], L"\text{Parity}", padding=(0, 0, 2, -10))
Label(gs[1, 2, Top()], L"G_{\text{nl}}", padding=(0, 0, 2, -10))


colsize!(g, 2, Auto(0.4))
colsize!(g, 2, Auto(0.9))
cb.alignmode = Mixed(top=-15, bottom=0)
f |> display

##
save(plotsdir(string("uhplot_transparent", fixedparamstring, ".png")), f, px_per_unit=4)

######
##
ssdata = resultsUVap[1, :]
nx, ny = size(ssdata[:sweet_spots])
smallres = 200
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.55, 0.75), (0.25, 0.5), (0.22, 0.1), (0.25, 0.1)])
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.1, 0.6), (0.4, 0.2), (0.4, 0.4)])
sweet_spots = [ssdata[:sweet_spots][pos...] for pos in positions]
transport = Transport(QuantumDots.Pauli(), (; T=1 / 20, μ=(0.0, 0.0)))
paramstring = map(ss -> map(x -> round(x, digits=2), ss.parameters), sweet_spots)

csdata_small = [charge_stability_scan((; ss.parameters...), 1, 1, smallres; transport) for ss in sweet_spots]

##
colors = cgrad(:rainbow, categorical=true)[[1, 2, 5, 4][1:length(csdata_small)]]
fixedparamstring = map(x -> round(x, digits=2), ssdata[:fixedparams])
f = Figure(resolution=400 .* (1.4, 1), fontsize=20, backgroundcolor=:white);
f = Figure(resolution=400 .* (1.4, 1), fontsize=20, backgroundcolor=:transparent);
g = f[1, 1] = GridLayout()
gb = g[1, 1] = GridLayout()
gs = g[1, 2] = GridLayout()

xlabel = paramstyle[ssdata[:xlabel]]
ylabel = paramstyle[ssdata[:ylabel]]
ax = Axis(gb[1, 1]; xlabel, ylabel)
hm = plot_sweet_scan!(ax, ssdata)

cb = add_exp_colorbar!(gb[2, 1], hm; vertical=false, label=L"1-\mathrm{MP}", flipaxis=false, labelpadding=-5, height=12)

scatter!(ax, map(ss -> ss.parameters[ssdata.xlabel], sweet_spots), map(ss -> ss.parameters[ssdata.ylabel], sweet_spots),
    color=colors, markersize=20, marker=:xcross, strokewidth=2)

datamap = x -> sign(x.gap)
datamap2 = x -> real(x.conductance)[1, 2] * real(x.conductance)[1, 1]
datamap2 = x -> real(x.conductance)[1, 2] - real(x.conductance)[2, 1]
datamap2 = x -> real(x.conductance)[1, 2]
spinecolor = [NamedTuple(map(x -> x => colors[n], [:bottomspinecolor, :leftspinecolor, :topspinecolor, :rightspinecolor])) for n in eachindex(csdata_small)]
small_axes = [Axis(gs[n, 1]; spinecolor[n]..., spinewidth=4, aspect=1) for (n, data) in enumerate(csdata_small)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
small_axes2 = [Axis(gs[n, 2]; spinecolor[n]..., spinewidth=4, aspect=1) for (n, data) in enumerate(csdata_small)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
foreach((ax, data) -> plot_charge_stability!(ax, data; datamap), small_axes, csdata_small)
foreach((ax, data) -> plot_charge_stability!(ax, data; datamap=datamap2, colormap=:vik, colorrange=2), small_axes2, csdata_small)
hidedecorations!.(small_axes)
hidedecorations!.(small_axes2)
labels = [Label(gs[n, 1:2, Bottom()], L"MP ≈ %$(round(1-MPU(ss),digits=2))", padding=(0, 0, -15, 4)) for (n, ss) in enumerate(sweet_spots)]
Label(gs[1, 1, Top()], L"\text{Parity}", padding=(0, 0, 2, -10))
Label(gs[1, 2, Top()], L"G_{\text{nl}}", padding=(0, 0, 2, -10))

colsize!(g, 2, Auto(0.4))
colsize!(g, 2, Auto(0.9))
cb.alignmode = Mixed(top=-15, bottom=0)

f |> display

##
save(plotsdir(string("_uvplot_transparent", fixedparamstring, ".png")), f, px_per_unit=4)

