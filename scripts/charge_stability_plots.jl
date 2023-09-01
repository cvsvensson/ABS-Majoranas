using DrWatson
@quickactivate "ABS-Majoranas"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))

using DataFrames
# resultsUV = collect_results(datadir("UV-scan"))
resultsUh = collect_results(datadir("Uh-scan"))
##
bigres = 500 # 500 takes minutes
smallres = 500
resultsUh[1, :sweet_spots][2, 9].parameters
ssparams = resultsUh[1, :sweet_spots][2, 9].parameters #(; U=0, V=0, h=1.5, t=0.5, Δ=1.0, tratio=0.2)
csdata_big = charge_stability_scan((; ssparams..., μ1=0, μ2=0), 8, 8, bigres);
ϕs = ssparams.ϕ .+ [-pi / 2,0, pi / 2]
csdata_small = [charge_stability_scan((; ssparams..., ϕ), 1, 1, smallres) for ϕ in ϕs]

##
f = Figure(resolution=400 .* (1, 1), fontsize=20, backgroundcolor = :transparent);
f = Figure(resolution=400 .* (1, 1), fontsize=20, backgroundcolor = :white);
g = f[1, 1] = GridLayout()
gb = g[1, 1] = GridLayout()
gs = g[2, 1] = GridLayout()


datamap = x -> sign(x.gap)
ax = Axis(gb[1,1], xlabel = L"μ_1", ylabel = L"μ_2")#,  subtitle = L"\tanh{\left(δE\right)}")
hm = plot_charge_stability!(ax, csdata_big; colorrange = 1, datamap = x->tanh(x.gap), colormap = :berlin)
Colorbar(gb[1,2],hm)

boxax = Axis(gb[1,1])#, bbox = BBox(1, 2, 3, 4))
hidedecorations!(boxax)
linewidth = 2
linecolor = :black
lines!(boxax, ssparams.μ1 .+ [-1,1], ssparams.μ2 .+ [1,1]; color = linecolor, linewidth)
lines!(boxax, ssparams.μ1 .+ [-1,1], ssparams.μ2 .- [1,1]; color = linecolor, linewidth)
lines!(boxax, ssparams.μ1 .+ [-1,-1], ssparams.μ2 .+ [-1,1]; color = linecolor, linewidth)
lines!(boxax, ssparams.μ1 .+ [1,1], ssparams.μ2 .+ [-1,1]; color = linecolor, linewidth)
linkaxes!(boxax,ax)
xlims!(first(csdata_big[:μ1]), last(csdata_big[:μ1]))
ylims!(first(csdata_big[:μ2]), last(csdata_big[:μ2]))
# ax1 = Axis(gs[1,1], xlabel = L"μ_1", ylabel = "after", fontsize = 1)

spinecolor = NamedTuple(map(x-> x=> linecolor, [:bottomspinecolor, :leftspinecolor, :topspinecolor, :rightspinecolor]))
small_axes = [Axis(gs[1, n]; spinecolor..., spinewidth = linewidth) for (n, ϕ) in enumerate(ϕs)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
foreach((ax,data) -> plot_charge_stability!(ax, data; datamap), small_axes, csdata_small)
hidedecorations!.(small_axes)
ax.xticks = -8:4:8
ax.yticks = -8:4:8
rowgap!(g, 5)
rowsize!(g, 2, Auto(.35))
labels = [Label(gs[1, n,Bottom()], L"ϕ ≈ %$(round(ϕ,digits=2))", padding = (0,0,-5,2)) for (n, ϕ) in enumerate(ϕs)]

# Label(gb[1, 1, Top()], "tanh(δE)", valign = :top,
Label(gb[1, 1, Top()], L"\tanh{(δE)}", valign = :top,
    padding = (0, 0, 2, -10))
Label(gs[1, 1:3, Top()], L"\text{Ground state parity}", valign = :bottom,
    padding = (0, 0, 2, 0))


f |> display
##
paramstring = map(x->round(x,digits=2),ssparams)
##
save(plotsdir(string("charge_stability_phase", paramstring, ".pdf")), f, pt_per_unit=1)
save(plotsdir(string("charge_stability_phase", paramstring, ".png")), f, px_per_unit=4)
##
save(plotsdir(string("charge_stability_phase_transparent", paramstring, ".png")), f, px_per_unit=4)
save(plotsdir(string("charge_stability_phase_transparent", paramstring, ".pdf")), f, pt_per_unit=1)