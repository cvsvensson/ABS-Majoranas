using DrWatson
@quickactivate "Majorana sweet spot"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))
using DataFrames
resultsUVap = collect_results(datadir("UV-scan", "anti_parallel", "final"))
resultsUhap = collect_results(datadir("Uh-scan", "anti_parallel", "final"))
resultsUV = combine_results(resultsUVap)
resultsUh = combine_results(resultsUhap)
##
data1 = resultsUh
ϕs1 = map(x -> x.parameters.ϕ, data1[:sweet_spots])
mu1s1 = map(x -> x.parameters.μ1, data1[:sweet_spots])
mu2s1 = map(x -> x.parameters.μ2, data1[:sweet_spots])
data2 = resultsUV
ϕs2 = map(x -> x.parameters.ϕ, data2[:sweet_spots])
mu1s2 = map(x -> x.parameters.μ1, data2[:sweet_spots])
mu2s2 = map(x -> x.parameters.μ2, data2[:sweet_spots])
##
f = Figure(; resolution=500 .* (1.4, 1), fontsize=25, backgroundcolor=:white)
g = f[1, 1] = GridLayout()
gp = g[2, 1] = GridLayout()
gc = g[1, 1] = GridLayout()
flipaxis = true
xlabel1 = paramstyle[data1[:xlabel]]
ylabel1 = paramstyle[data1[:ylabel]]
xlabel2 = paramstyle[data2[:xlabel]]
ylabel2 = paramstyle[data2[:ylabel]]
aspect = 1
colormap = :batlow
x1 = data1[:x]
y1 = data1[:y]
x2 = data2[:x]
y2 = data2[:y]
tellwidth = false
cr1 = (0, pi)
cr2 = (-9, 0)
cr3 = (-1, 4)
labelspace = 20
cbarwidth = 165
cbarheight = 10
ax1_1 = Axis(gp[1, 1]; ylabel=ylabel1, aspect, yticks=0:2:10)#, alignmode=alignmodeax)
hm = heatmap!(ax1_1, x1, y1, ϕs1; colorrange=cr1, tellwidth, colormap)
cb = Colorbar(gc[1, 1,], hm; height=cbarheight, width=cbarwidth, tellwidth, vertical=false, label=L"\phi", flipaxis)
cb.ticklabelspace = labelspace
ax2_1 = Axis(gp[1, 2]; aspect)
hm = heatmap!(ax2_1, x1, y1, -mu1s1; colorrange=cr2, colormap)
cb = Colorbar(gc[1, 2,], hm; height=cbarheight, width=cbarwidth, tellwidth, vertical=false, label=paramstyle[:μ1], flipaxis)
cb.ticklabelspace = labelspace
ax3_1 = Axis(gp[1, 3]; aspect)
hm = heatmap!(ax3_1, x1, y1, -mu2s1; colorrange=cr3, colormap)
cb = Colorbar(gc[1, 3,], hm; height=cbarheight, width=cbarwidth, tellwidth, vertical=false, label=paramstyle[:μ2], flipaxis, ticks=-1:2:4)
cb.ticklabelspace = labelspace
hideydecorations!(ax2_1)
hideydecorations!(ax3_1)
hidexdecorations!.((ax1_1, ax2_1, ax3_1))

ax1_2 = Axis(gp[2, 1]; ylabel=ylabel2, aspect, tellwidth, yticks=0:1:1)#, alignmode=alignmodeax)
hm = heatmap!(ax1_2, x2, y2, ϕs2; colorrange=cr1, colormap, tellwidth)
ax2_2 = Axis(gp[2, 2]; xlabel = xlabel2, aspect, tellwidth)#, alignmode=alignmodeax)
hm = heatmap!(ax2_2, x2, y2, -mu1s2; colorrange=cr2, colormap)
ax3_2 = Axis(gp[2, 3]; aspect, tellwidth)#, alignmode=alignmodeax)
hm = heatmap!(ax3_2, x2, y2, -mu2s2; colorrange=cr3, colormap)
hideydecorations!(ax2_2)
hideydecorations!(ax3_2)
rowgap!(gp, 1, 20)

rowsize!(g, 1, Relative(-0.08))
f
##
save(plotsdir(string("parameters.png")), f, px_per_unit=4)
