using DrWatson
@quickactivate "Majorana sweet spot"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))
using DataFrames
# resultsUVp = collect_results(datadir("UV-scan", "parallel"))
# resultsUhp = collect_results(datadir("Uh-scan", "parallel"))
# Methods = [:generating_set_search, :probabilistic_descent, :adaptive_de_rand_1_bin_radiuslimited]
resultsUVap = collect_results(datadir("UV-scan", "anti_parallel", "final"))
resultsUhap = collect_results(datadir("Uh-scan", "anti_parallel", "final"))

resultsUh = combine_results(resultsUhap)
resultsUV = combine_results(resultsUVap)

##
for n in 1:4
    scanfig = let data1 = resultsUhap[n, :], data2 = resultsUVap[n, :]
        fig = Figure(; resolution=400 .* (2, 1), fontsize=20, backgroundcolor=:white)
        xlabel1 = paramstyle[data1[:xlabel]]
        xlabel2 = paramstyle[data2[:xlabel]]
        ylabel1 = paramstyle[data1[:ylabel]]
        ylabel2 = paramstyle[data2[:ylabel]]
        ax1 = Axis(fig[1, 1]; xlabel=xlabel1, ylabel=ylabel1)
        ax2 = Axis(fig[1, 2]; xlabel=xlabel2, ylabel=ylabel2)
        g = fig[1, 1:3] = GridLayout()
        hm1 = plot_sweet_scan!(ax1, data1)#; datamap = x->abs(x.gap), colorrange = (1e-5,1))
        hm2 = plot_sweet_scan!(ax2, data2)#; datamap = x->abs(x.gap), colorrange = (1e-5,1))
        label = L"1-\text{MP}"
        add_exp_colorbar!(fig[1, 3], hm1;)
        paramstring1 = map(x -> round.(x, digits=2), data1[:fixedparams])
        paramstring2 = map(x -> round.(x, digits=2), data2[:fixedparams])

        supertitle = Label(fig[1, 3, Top()], label, fontsize=25, padding=(0, 0, 8, 0))
        supertitle2 = Label(fig[1, 1, Top()], string(data1[:sweet_spots][1].optimization.Method), fontsize=20, padding=(0, 0, 8, 0))
        supertitle2 = Label(fig[1, 2, Top()], string(data2[:sweet_spots][1].optimization.Method), fontsize=20, padding=(0, 0, 8, 0))
        colgap!(fig.layout, 1, 16)
        colgap!(fig.layout, 2, 10)
        #save(plotsdir(string("uhvplot_transparent", paramstring1, "_", paramstring2, ".png")), fig, px_per_unit=4)
        fig |> display
    end
end

##
ssdata = resultsUh
nx, ny = size(ssdata[:sweet_spots])
smallres = 100
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.4, 0.7), (0.25, 0.5), (0.25, 0.38), (0.25, 0.33)])
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.5, 0.8), (0.3, 0.4), (0.1, 0.1)])
sweet_spots = [ssdata[:sweet_spots][pos...] for pos in positions]
paramstring = map(ss -> map(x -> round(x, digits=2), ss.parameters), sweet_spots)
# transport = Transport(QuantumDots.Pauli(), (; T=1 / 20, μ=(0.0, 0.0)))
transport = missing
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
spinecolor = [NamedTuple(map(x -> x => colors[n], [:bottomspinecolor, :leftspinecolor, :topspinecolor, :rightspinecolor])) for n in eachindex(csdata_small)]
small_axes = [Axis(gs[n, 1]; spinecolor[n]..., spinewidth=4, aspect=1, xlabel= paramstyle[:μ1]) for (n, data) in enumerate(csdata_small)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
foreach((ax, data) -> plot_charge_stability!(ax, data; datamap), small_axes, csdata_small)

hidedecorations!.(small_axes)
labels = [Label(gs[n, 1, Bottom()], L"MP ≈ %$(round(1-MPU(ss),digits=2))", padding=(0, 0, -15, 4)) for (n, ss) in enumerate(sweet_spots)]
Label(gb[1, 1, TopLeft()], L"a)", padding=(0, 25, 0, -5), tellheight=false, tellwidth=false)
Label(gs[1, 1, TopLeft()], L"b)", padding=(0, 10, 0, -5), tellheight=false, tellwidth=false)

colsize!(g, 2, Auto(0.4))
cb.alignmode = Mixed(top=-15, bottom=0)
f |> display

##
save(plotsdir(string("uhplot_transparent", fixedparamstring, ".png")), f, px_per_unit=4)
######
##
ssdata = resultsUV
nx, ny = size(ssdata[:sweet_spots])
smallres = 100
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.55, 0.75), (0.25, 0.5), (0.22, 0.1), (0.25, 0.1)])
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.1, 0.6), (0.4, 0.2), (0.6, 0.6)])
sweet_spots = [ssdata[:sweet_spots][pos...] for pos in positions]
# transport = Transport(QuantumDots.Pauli(), (; T=1 / 20, μ=(0.0, 0.0)))
transport = missing
paramstring = map(ss -> map(x -> round(x, digits=2), ss.parameters), sweet_spots)

csdata_small = [charge_stability_scan((; ss.parameters..., ϕ=ss.parameters.ϕ, μ1=ss.parameters.μ1, μ2=ss.parameters.μ2 ), 1.5, 1.5, smallres; transport) for ss in sweet_spots]

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
spinecolor = [NamedTuple(map(x -> x => colors[n], [:bottomspinecolor, :leftspinecolor, :topspinecolor, :rightspinecolor])) for n in eachindex(csdata_small)]
small_axes = [Axis(gs[n, 1]; spinecolor[n]..., spinewidth=4, aspect=1, xlabel= paramstyle[:μ1]) for (n, data) in enumerate(csdata_small)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
foreach((ax, data) -> plot_charge_stability!(ax, data; datamap), small_axes, csdata_small)

hidedecorations!.(small_axes)
labels = [Label(gs[n, 1, Bottom()], L"MP ≈ %$(round(1-MPU(ss),digits=2))", padding=(0, 0, -15, 4)) for (n, ss) in enumerate(sweet_spots)]
Label(gb[1, 1, TopLeft()], L"a)", padding=(0, 45, 0, -5), tellheight=false, tellwidth=false)
Label(gs[1, 1, TopLeft()], L"b)", padding=(0, 10, 0, -5), tellheight=false, tellwidth=false)

colsize!(g, 2, Auto(0.4))
cb.alignmode = Mixed(top=-15, bottom=0)
f |> display


##
save(plotsdir(string("uvplot_transparent", fixedparamstring, ".png")), f, px_per_unit=4)

