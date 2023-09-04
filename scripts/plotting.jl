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
scanfig = let data1 = resultsUhap[1, :], data2 = resultsUVap[1, :]
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
    # add_exp_colorbar!(fig[1,3], hm)
    # ax, hm = plot_sweet_scan!(fig[1, 1], data1; colorrange)
    # ax, hm = plot_sweet_scan!(fig[1, 2], data2; colorrange)
    # add_colorbar!(fig[1, 1+end], hm; colorrange)
    paramstring1 = map(x -> round(x, digits=2), data1[:fixedparams])
    paramstring2 = map(x -> round(x, digits=2), data2[:fixedparams])

    supertitle = Label(fig[1, 3, Top()], label, fontsize=25, padding=(0, 0, 8, 0))
    colgap!(fig.layout, 1, 16)
    colgap!(fig.layout, 2, 10)
    save(plotsdir(string("____uhvplot_transparent", paramstring1, "_", paramstring2, ".png")), fig, px_per_unit=4)
    fig
end


##
nx, ny = size(resultsUhap[1, :sweet_spots])
smallres = 200
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.5, 0.8), (0.25, 0.4), (0.1, 0.1)])
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.55, 0.75), (0.25, 0.5), (0.22, 0.1), (0.25, 0.1)])
ssparams = [resultsUhap[1, :sweet_spots][pos...].parameters for pos in positions]
paramstring = map(ss -> map(x -> round(x, digits=2), ss), ssparams)
csdata_small = [charge_stability_scan((; params...), 1, 1, smallres) for params in ssparams]

##
colors = Makie.wong_colors()[[2, 4, 6]]
colors = cgrad(:Dark2_7, categorical=true)[[2, 3, 4]]
colors = cgrad(:rainbow, categorical=true)[[1, 2, 5, 4]]
ssdata = resultsUhap[1, :]
fixedparamstring = map(x->round(x,digits=2), ssdata[:fixedparams])
f = Figure(resolution=400 .* (1, .9), fontsize=20, backgroundcolor=:white);
f = Figure(resolution=400 .* (1, 0.9), fontsize=20, backgroundcolor=:transparent);
g = f[1, 1] = GridLayout()
gb = g[1, 1] = GridLayout()
gs = g[1, 2] = GridLayout()

datamap = x -> sign(x.gap)
# ax = Axis(gb[1, 1], xlabel=L"μ_1", ylabel=L"μ_2", aspect = 1)#,  subtitle = L"\tanh{\left(δE\right)}")
xlabel = paramstyle[ssdata[:xlabel]]
ylabel = paramstyle[ssdata[:ylabel]]
ax = Axis(gb[1, 1]; xlabel, ylabel)
hm = plot_sweet_scan!(ax, ssdata)
# ax, hm = plot_sweet_scan!(f, ssdata, (1, 1); colorbar=false)

# hm = plot_charge_stability!(ax, csdata_big; colorrange=1, datamap=x -> tanh(x.gap), colormap=:berlin)
cb = add_exp_colorbar!(gb[2, 1], hm; vertical=false, label=L"1-\mathrm{MP}", flipaxis=false, labelpadding=-5, height=12,
    alignmode=Outside())
# cb = Colorbar(gb[2, 1], hm; vertical=false, label=L"1-\mathrm{MP}", flipaxis=false,labelpadding = -20, ticks = [0,1],height = 12,
# alignmode = Outside())

# pointax = Axis(gb[1, 1])
scatter!(ax, map(params -> params[ssdata.xlabel], ssparams), map(params -> params[ssdata.ylabel], ssparams),
    color=colors, markersize=20, marker=:xcross, strokewidth=2)
# hidedecorations!(pointax)
# linkaxes!(pointax, ax)

spinecolor = [NamedTuple(map(x -> x => colors[n], [:bottomspinecolor, :leftspinecolor, :topspinecolor, :rightspinecolor])) for n in eachindex(csdata_small)]
small_axes = [Axis(gs[n, 1]; spinecolor[n]..., spinewidth=4, aspect=1) for (n, data) in enumerate(csdata_small)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
foreach((ax, data) -> plot_charge_stability!(ax, data; datamap), small_axes, csdata_small)
hidedecorations!.(small_axes)
# ax.xticks = -4:8:4
# ax.yticks = -4:8:4
# ax.xticklabelspace = 0.0#tight_xticklabel_spacing!(ax)
# ax.yticklabelspace = 0.0#tight_xticklabel_spacing!(ax)
# rowgap!(g, 5)
#labels = [Label(gs[n, 1, Bottom()], L"ϕ ≈ %$(round(ϕ,digits=2))", padding=(0, 0, -5, 2)) for (n, ϕ) in enumerate(ϕs)]

# Label(gb[1, 1, Top()], "tanh(δE)", valign = :top,
# Label(gb[1, 1, Top()], L"\tanh{(δE)}", valign = :top,
# padding = (0, 0, 2, -10))

# Label(gs[1, 1, Top()], L"\text{Parity}", valign=:bottom,
# padding=(0, 0, 3, -10))

colsize!(g, 2, Auto(0.4))
# rowsize!(gb, 2, Auto(.1))
cb.alignmode = Mixed(top=-15, bottom=0)
# rowgap!(gb, 10)
# trim!(f.layout)
f |> display

##
save(plotsdir(string("uhplot_transparent",fixedparamstring, ".png")), f, px_per_unit=4)
