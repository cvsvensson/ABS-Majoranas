using DrWatson
@quickactivate "ABS-Majoranas"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))

using DataFrames
# resultsUV = collect_results(datadir("UV-scan"))
resultsUh = collect_results(datadir("Uh-scan"))
##
resultsUh[1, :sweet_spots][2, 9].parameters
ssparams = resultsUh[1, :sweet_spots][2, 9].parameters #(; U=0, V=0, h=1.5, t=0.5, Δ=1.0, tratio=0.2)
csdata_big = charge_stability_scan((; ssparams..., μ1=0, μ2=0), 8, 8, 100);
ϕs = ssparams.ϕ .+ [0, pi / 2, -pi / 2]
csdata_small = [charge_stability_scan((; ssparams..., ϕ), 1, 1, 100) for ϕ in ϕs]

##
f = Figure(resolution=400 .* (1, 1.2), fontsize=20);
g = f[1, 1] = GridLayout()
gb = g[1, 1] = GridLayout()
gs = g[2, 1] = GridLayout()


datamap = x -> sign(x.gap)
ax, hm = plot_charge_stability!(gb[1, 1], csdata_big; datamap)

boxax = Axis(gb[1,1])#, bbox = BBox(1, 2, 3, 4))
hidedecorations!(boxax)
lines!(boxax, ssparams.μ1 .+ [-1,1], ssparams.μ2 .+ [1,1]; color = :black)
lines!(boxax, ssparams.μ1 .+ [-1,1], ssparams.μ2 .- [1,1]; color = :black)
lines!(boxax, ssparams.μ1 .+ [-1,-1], ssparams.μ2 .+ [-1,1]; color = :black)
lines!(boxax, ssparams.μ1 .+ [1,1], ssparams.μ2 .+ [-1,1]; color = :black)
linkaxes!(boxax,ax)
xlims!(first(csdata_big[:μ1]), last(csdata_big[:μ1]))
ylims!(first(csdata_big[:μ2]), last(csdata_big[:μ2]))
# ax1 = Axis(gs[1,1], xlabel = L"μ_1", ylabel = "after", fontsize = 1)

small_axes = [Axis(gs[1, n]; subtitle = L"ϕ=%$(round(ϕ,digits=2))" ) for (n, ϕ) in enumerate(ϕs)]
foreach((ax,data) -> plot_charge_stability!(ax, data; datamap), small_axes, csdata_small)
hidedecorations!.(small_axes)
ax.xticks = -8:4:8
ax.yticks = -8:4:8
rowgap!(g, 5)
rowsize!(g, 2, Auto(.33))


# rects = fig[1:4, 1:6] = [
#     Box(fig, color = c)
#     for c in get.(Ref(ColorSchemes.rainbow), (0:23) ./ 23)]
# Box(f[1,2], color = :red, height = 10,width = 10)
# axbox = Axis(f, bbox = BBox(140, 250, 200, 300))
# Box(axbox, color = :red)

f |> display