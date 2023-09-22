using DrWatson
@quickactivate "Majorana sweet spot"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))
using DataFrames
resultsUVap = collect_results(datadir("UV-scan", "anti_parallel","final"))
resultsUhap = collect_results(datadir("Uh-scan", "anti_parallel","final"))
resultUV = combine_results(resultsUVap)
resultUh = combine_results(resultsUhap)
##
data = resultsUV
ϕs = map(x -> x.parameters.ϕ, data[:sweet_spots])
mu1s = map(x -> x.parameters.μ1, data[:sweet_spots])
mu2s = map(x -> x.parameters.μ2, data[:sweet_spots])
##

f = Figure(; resolution=300 .* (1.8, 1), fontsize=25, backgroundcolor=:white)
g = f[1, 1] = GridLayout()
gp = g[2, 1] = GridLayout()
gc = g[1, 1] = GridLayout()
flipaxis = true
height = 10
width = 135
xlabel = paramstyle[data[:xlabel]]
ylabel = paramstyle[data[:ylabel]]
aspect = 1
colormap = :batlow
x = data[:x]
y = data[:y]
tellwidth = false
# alignmode = Inside()
alignmode = Mixed(top=100, bottom=0)
alignmodeax = Mixed(top=0, bottom=-60)
ax1 = Axis(gp[1, 1]; xlabel, ylabel, aspect, tellwidth, yticks=0:1, alignmode=alignmodeax)
hm = heatmap!(ax1, x, y, ϕs; colorrange=(0, pi), colormap=:batlow, tellwidth)
cb = Colorbar(gc[1, 1], hm; alignmode, tellwidth, vertical=false, label=L"\phi", flipaxis, height, width)

ax2 = Axis(gp[1, 2]; xlabel, ylabel, aspect, tellwidth, alignmode=alignmodeax)
hm = heatmap!(ax2, x, y, mu1s; colormap)
cb = Colorbar(gc[1, 2], hm; alignmode, tellwidth, vertical=false, label=paramstyle[:μ1], flipaxis, height, width)


ax3 = Axis(gp[1, 3]; xlabel, ylabel, aspect, tellwidth, alignmode=alignmodeax)
hm = heatmap!(ax3, x, y, mu2s; colormap)
cb = Colorbar(gc[1, 3], hm; alignmode, tellwidth, vertical=false, label=paramstyle[:μ2], flipaxis, height, ticks=-2:1:2, width)
hideydecorations!(ax2)
hideydecorations!(ax3)
colgap!(gp, 25)
colgap!(gc, 25)

rowsize!(g, 1, Relative(-0.2))
f
##
save(plotsdir(string("UV_parameters.png")), f, px_per_unit=4)

########

##
data = resultsUh
ϕs = map(x -> x.parameters.ϕ, data[:sweet_spots])
mu1s = map(x -> x.parameters.μ1, data[:sweet_spots])
mu2s = map(x -> x.parameters.μ2, data[:sweet_spots])
##
f = Figure(; resolution=300 .* (1.8, 1), fontsize=25, backgroundcolor=:white)
g = f[1, 1] = GridLayout()
gp = g[2, 1] = GridLayout()
gc = g[1, 1] = GridLayout()
flipaxis = true
height = 10
width = 135
xlabel = paramstyle[data[:xlabel]]
ylabel = paramstyle[data[:ylabel]]
aspect = 1
colormap = :batlow
x = data[:x]
y = data[:y]
tellwidth = false
alignmode = Mixed(top=100, bottom=0)
alignmodeax = Mixed(top=0, bottom=-60)
ax1 = Axis(gp[1, 1]; xlabel, ylabel, aspect, tellwidth, yticks=0:2:10, alignmode=alignmodeax)
hm = heatmap!(ax1, x, y, ϕs; colorrange=(0, pi), colormap=:batlow, tellwidth)
cb = Colorbar(gc[1, 1], hm; alignmode, tellwidth, vertical=false, label=L"\phi", flipaxis, height, width)

ax2 = Axis(gp[1, 2]; xlabel, ylabel, aspect, tellwidth, alignmode=alignmodeax)
hm = heatmap!(ax2, x, y, mu1s; colormap)
cb = Colorbar(gc[1, 2], hm; alignmode, tellwidth, vertical=false, label=paramstyle[:μ1], flipaxis, height, width)


ax3 = Axis(gp[1, 3]; xlabel, ylabel, aspect, tellwidth, alignmode=alignmodeax)
hm = heatmap!(ax3, x, y, mu2s; colormap)
cb = Colorbar(gc[1, 3], hm; alignmode, tellwidth, vertical=false, label=paramstyle[:μ2], flipaxis, height, ticks=-12:3:10, width)
hideydecorations!(ax2)
hideydecorations!(ax3)
colgap!(gp, 25)
colgap!(gc, 25)

rowsize!(g, 1, Relative(-0.2))

f
##
save(plotsdir(string("Uh_parameters.png")), f, px_per_unit=4)
