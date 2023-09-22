using DrWatson
@quickactivate "ABS-Majoranas"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))

using DataFrames
resultsUVap = collect_results(datadir("UV-scan", "anti_parallel","final"))
resultsUhap = collect_results(datadir("Uh-scan", "anti_parallel","final"))
resultUV = combine_results(resultsUVap)
resultUh = combine_results(resultsUhap)

##
bigres = 200
smallres = 200
data = resultsUh
pos = Int.(round.(size(data[:sweet_spots]) .* (20/40, 9/40)))
ssparams = data[:sweet_spots][pos...].parameters
paramstring = map(x -> round(x, digits=2), ssparams)
μ0 = (ssparams.μ1 + ssparams.μ2)/2
csdata_big = charge_stability_scan((; ssparams..., μ1=μ0, μ2=μ0), 8, 8, bigres);
ϕs = ssparams.ϕ .+ [-pi / 2, 0, pi / 2]
csdata_small = [charge_stability_scan((; ssparams..., ϕ), 1.5, 1.5, smallres) for ϕ in ϕs]

## Horizontal fig
f = Figure(resolution=400 .* (1, 0.9), fontsize=20, backgroundcolor=:white);
f = Figure(resolution=400 .* (1, 0.9), fontsize=20, backgroundcolor=:transparent);
g = f[1, 1] = GridLayout()
gb = g[1, 1] = GridLayout()
gs = g[1, 2] = GridLayout()


datamap = x -> sign(x.gap)
ax = Axis(gb[1, 1], xlabel=paramstyle[:μ1], ylabel=paramstyle[:μ2], aspect=1)#,  subtitle = L"\tanh{\left(δE\right)}")
hm = plot_charge_stability!(ax, csdata_big; colorrange=(-1,1), datamap=x -> tanh(x.gap), colormap=:redsblues)

cb = Colorbar(gb[2, 1], hm; vertical=false, label=L"\tanh{(δE)}", flipaxis=false, labelpadding=-20, ticks=[-1, 1], height=12,
    alignmode=Outside())

boxax = Axis(gb[1, 1], aspect=1)
hidedecorations!(boxax)
linewidth = 2
linecolor = :black
lines!(boxax, ssparams.μ1 .+ [-1, 1], ssparams.μ2 .+ [1, 1]; color=linecolor, linewidth)
lines!(boxax, ssparams.μ1 .+ [-1, 1], ssparams.μ2 .- [1, 1]; color=linecolor, linewidth)
lines!(boxax, ssparams.μ1 .+ [-1, -1], ssparams.μ2 .+ [-1, 1]; color=linecolor, linewidth)
lines!(boxax, ssparams.μ1 .+ [1, 1], ssparams.μ2 .+ [-1, 1]; color=linecolor, linewidth)
linkaxes!(boxax, ax)
xlims!(first(csdata_big[:μ1]), last(csdata_big[:μ1]))
ylims!(first(csdata_big[:μ2]), last(csdata_big[:μ2]))
# ax1 = Axis(gs[1,1], xlabel = L"μ_1", ylabel = "after", fontsize = 1)

spinecolor = NamedTuple(map(x -> x => linecolor, [:bottomspinecolor, :leftspinecolor, :topspinecolor, :rightspinecolor]))
small_axes = [Axis(gs[n, 1]; spinecolor..., spinewidth=linewidth, aspect=1) for (n, ϕ) in enumerate(ϕs)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
foreach((ax, data) -> plot_charge_stability!(ax, data; datamap), small_axes, csdata_small)
hidedecorations!.(small_axes)
ax.xticks = -2:2:4
ax.yticks = -2:2:4
ax.xticklabelspace = 15.0#tight_xticklabel_spacing!(ax)
ax.yticklabelspace = 5.0#tight_xticklabel_spacing!(ax)
# rowgap!(g, 5)
# labels = [Label(gs[n, 1, Bottom()], L"ϕ ≈ %$(round(ϕ,digits=2))", padding=(0, 0, -5, 2)) for (n, ϕ) in enumerate(ϕs)]
[Label(gs[n, 1, Bottom()], L"\delta ϕ %$s \delta ϕ_\star", padding=(0, 0, -5, 2)) for (n, s) in enumerate(["<","=",">"])]

# Label(gb[1, 1, Top()], "tanh(δE)", valign = :top,
# Label(gb[1, 1, Top()], L"\tanh{(δE)}", valign = :top,
# padding = (0, 0, 2, -10))
Label(gs[1, 1, Top()], L"\text{Parity}", valign=:bottom,
    padding=(0, 0, 3, -10))

colsize!(g, 2, Auto(0.4))
# rowsize!(gb, 2, Auto(.1))
cb.alignmode = Mixed(top=-10, bottom=-10)
# rowgap!(gb, 10)
# trim!(f.layout)
Label(gb[1, 1, TopLeft()], L"a)", padding=(0, 25, 0, -5), tellheight=false, tellwidth=false)
Label(gs[1, 1, TopLeft()], L"b)", padding=(0, 10, 0, -5), tellheight=false, tellwidth=false)


f |> display
##
save(plotsdir(string("horizontal_charge_stability_phase", paramstring, ".pdf")), f, pt_per_unit=1)
save(plotsdir(string("horizontal_charge_stability_phase", paramstring, ".png")), f, px_per_unit=4)
##
save(plotsdir(string("horizontal_charge_stability_phase_transparent", paramstring, ".png")), f, px_per_unit=4)
save(plotsdir(string("horizontal_charge_stability_phase_transparent", paramstring, ".pdf")), f, pt_per_unit=1)
