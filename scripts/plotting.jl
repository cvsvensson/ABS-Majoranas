using DrWatson
@quickactivate "Majorana sweet spot"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))
using DataFrames
resultsUV = collect_results(datadir("UV-scan"))
resultsUh = collect_results(datadir("Uh-scan"))
resultsUhpar = collect_results(datadir("Uh-scan","par"))

UVfig, _, _ = plot_sweet_scan(resultsUV[1, :])
display(UVfig)

Uhfig, _, _ = plot_sweet_scan(resultsUh[1, :])
display(Uhfig)
Uhfig, _, _ = plot_sweet_scan(resultsUhpar[1, :])
display(Uhfig)
fig, _, _ = plot_sweet_scan(resultsUh[1, :]; datamap=x -> x.parameters.μ1 - x.parameters.U - x.parameters.h, colorscale=identity, colorrange=(-2, 2), colorbar=false)
display(fig)
fig, _, _ = plot_sweet_scan(resultsUh[1, :]; datamap=x -> x.parameters.μ2 + x.parameters.h, colorscale=identity, colorrange=(-2, 2), colorbar=false)
display(fig)
fig, _, _ = plot_sweet_scan(resultsUh[1, :]; datamap=x -> x.parameters.ϕ, colorscale=identity, colorrange=(0, pi), colorbar=false)
display(fig)