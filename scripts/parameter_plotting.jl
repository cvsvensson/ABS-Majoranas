using DrWatson
@quickactivate "Majorana sweet spot"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))
using DataFrames
resultsUVap = collect_results(datadir("UV-scan", "anti_parallel", "final"))
resultsUhap = collect_results(datadir("Uh-scan", "anti_parallel", "final"))
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
alignmode = Mixed(top=100, bottom=-10)
alignmodeax = Mixed(top=0, bottom=-50)
ax1 = Axis(gp[1, 1]; ylabel, aspect, tellwidth, yticks=0:1:10, alignmode=alignmodeax)
hm = heatmap!(ax1, x, y, ϕs; colorrange=(0, pi), colormap=:batlow, tellwidth)
cb = Colorbar(gc[1, 1], hm; alignmode, tellwidth, vertical=false, label=L"\phi", flipaxis, height, width)

ax2 = Axis(gp[1, 2]; ylabel, aspect, tellwidth, alignmode=alignmodeax)
hm = heatmap!(ax2, x, y, -mu1s; colormap)
cb = Colorbar(gc[1, 2], hm; alignmode, tellwidth, vertical=false, label=paramstyle[:μ1], flipaxis, height, width)

ax3 = Axis(gp[1, 3]; ylabel, aspect, tellwidth, alignmode=alignmodeax)
hm = heatmap!(ax3, x, y, -mu2s; colormap)
cb = Colorbar(gc[1, 3], hm; alignmode, tellwidth, vertical=false, label=paramstyle[:μ2], flipaxis, height, ticks=-1:1, width)
hideydecorations!(ax2)
hideydecorations!(ax3)
linkaxes!(ax1, ax2, ax3)
Label(gp[1, 1:3, Bottom()], xlabel, tellheight=false, padding=(0, 0, 45, 0))
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
alignmode = Mixed(top=100, bottom=-10)
alignmodeax = Mixed(top=0, bottom=-50)
ax1 = Axis(gp[1, 1]; ylabel, aspect, tellwidth, yticks=0:2:10, alignmode=alignmodeax)
hm = heatmap!(ax1, x, y, ϕs; colorrange=(0, pi), colormap=:batlow, tellwidth)
cb = Colorbar(gc[1, 1], hm; alignmode, tellwidth, vertical=false, label=L"\phi", flipaxis, height, width)

ax2 = Axis(gp[1, 2]; ylabel, aspect, tellwidth, alignmode=alignmodeax)
hm = heatmap!(ax2, x, y, -mu1s; colormap)
cb = Colorbar(gc[1, 2], hm; alignmode, tellwidth, vertical=false, label=paramstyle[:μ1], flipaxis, height, width)


ax3 = Axis(gp[1, 3]; ylabel, aspect, tellwidth, alignmode=alignmodeax)
hm = heatmap!(ax3, x, y, -mu2s; colormap)
cb = Colorbar(gc[1, 3], hm; alignmode, tellwidth, vertical=false, label=paramstyle[:μ2], flipaxis, height, ticks=-12:3:10, width)
hideydecorations!(ax2)
hideydecorations!(ax3)
Label(gp[1, 1:3, Bottom()], xlabel, tellheight=false, padding=(0, 0, 45, 0))

colgap!(gp, 25)
colgap!(gc, 25)

rowsize!(g, 1, Relative(-0.2))

f
##
save(plotsdir(string("Uh_parameters.png")), f, px_per_unit=4)


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
# height = 10
# width = 135
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
# alignmode = Mixed(top=-50, bottom=0)
# alignmodeax = Mixed(top=0, bottom=-50)
cr1 = (0, pi)
cr2 = (-9, 0)
cr3 = (-1, 4)
labelspace = 20
cbarwidth = 165
cbarheight = 10
# labelpadding = 10
ax1_1 = Axis(gp[1, 1]; ylabel=ylabel1, aspect, yticks=0:2:10)#, alignmode=alignmodeax)
hm = heatmap!(ax1_1, x1, y1, ϕs1; colorrange=cr1, tellwidth, colormap)
cb = Colorbar(gc[1, 1,], hm; height=cbarheight, width=cbarwidth, tellwidth, vertical=false, label=L"\phi", flipaxis)
cb.ticklabelspace = labelspace
# cb.labelpadding = labelpadding
ax2_1 = Axis(gp[1, 2]; aspect)#, alignmode=alignmodeax)
hm = heatmap!(ax2_1, x1, y1, -mu1s1; colorrange=cr2, colormap)
cb = Colorbar(gc[1, 2,], hm; height=cbarheight, width=cbarwidth, tellwidth, vertical=false, label=paramstyle[:μ1], flipaxis)
cb.ticklabelspace = labelspace
# cb.labelpadding = labelpadding
ax3_1 = Axis(gp[1, 3]; aspect)#, alignmode=alignmodeax)
hm = heatmap!(ax3_1, x1, y1, -mu2s1; colorrange=cr3, colormap)
cb = Colorbar(gc[1, 3,], hm; height=cbarheight, width=cbarwidth, tellwidth, vertical=false, label=paramstyle[:μ2], flipaxis, ticks=-1:2:4)
cb.ticklabelspace = labelspace
# cb.labelpadding = labelpadding
hideydecorations!(ax2_1)
hideydecorations!(ax3_1)
hidexdecorations!.((ax1_1, ax2_1, ax3_1))

ax1_2 = Axis(gp[2, 1]; ylabel=ylabel2, aspect, tellwidth, yticks=0:1:1)#, alignmode=alignmodeax)
hm = heatmap!(ax1_2, x2, y2, ϕs2; colorrange=cr1, colormap, tellwidth)
ax2_2 = Axis(gp[2, 2]; xlabel, aspect, tellwidth)#, alignmode=alignmodeax)
hm = heatmap!(ax2_2, x2, y2, -mu1s2; colorrange=cr2, colormap)
ax3_2 = Axis(gp[2, 3]; aspect, tellwidth)#, alignmode=alignmodeax)
hm = heatmap!(ax3_2, x2, y2, -mu2s2; colorrange=cr3, colormap)
hideydecorations!(ax2_2)
hideydecorations!(ax3_2)

# Label(gp[2, 1:3, Bottom()], xlabel1, tellheight=true, padding=(0, 0, 10, 0))

# colgap!(gp, -100)
# colgap!(gc, -100)
# rowgap!(gp, 10)
rowgap!(gp, 1, 20)
# colgap!(gc, 25)

rowsize!(g, 1, Relative(-0.08))
f
##
save(plotsdir(string("parameters.png")), f, px_per_unit=4)
