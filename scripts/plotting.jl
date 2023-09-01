using DrWatson
@quickactivate "Majorana sweet spot"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))
using DataFrames
resultsUVp = collect_results(datadir("UV-scan","parallel"))
resultsUVap = collect_results(datadir("UV-scan","anti_parallel"))
resultsUhp = collect_results(datadir("Uh-scan","parallel"))
resultsUhap = collect_results(datadir("Uh-scan","anti_parallel"))


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
    fig = Figure(; resolution=400 .* (2, 1),fontsize=20, backgroundcolor=:transparent)
    xlabel1 = paramstyle[data1[:xlabel]]
    xlabel2 = paramstyle[data2[:xlabel]]
    ylabel1= paramstyle[data1[:ylabel]]
    ylabel2= paramstyle[data2[:ylabel]]
    ax1 = Axis(fig[1,1]; xlabel = xlabel1, ylabel = ylabel1)
    ax2 = Axis(fig[1,2]; xlabel = xlabel2, ylabel = ylabel2)
    g = fig[1,1:3] = GridLayout()
    hm1 = plot_sweet_scan!(ax1, data1)
    hm2 = plot_sweet_scan!(ax2, data2)
    label =  L"1-\text{MP}"
    add_exp_colorbar!(fig[1,3], hm1;)
    # add_exp_colorbar!(fig[1,3], hm)
    # ax, hm = plot_sweet_scan!(fig[1, 1], data1; colorrange)
    # ax, hm = plot_sweet_scan!(fig[1, 2], data2; colorrange)
    # add_colorbar!(fig[1, 1+end], hm; colorrange)
    paramstring1 = map(x -> round(x, digits=2), data1[:fixedparams])
    paramstring2 = map(x -> round(x, digits=2), data2[:fixedparams])

    supertitle = Label(fig[1, 3,Top()], label, fontsize=25, padding = (0,0,8,0))
    colgap!(fig.layout,1,16)
    colgap!(fig.layout,2,10)
    save(plotsdir(string("____uhvplot_transparent", paramstring1, "_", paramstring2, ".png")), fig, px_per_unit=4)
    fig
end