using DrWatson
@quickactivate "ABS-Majoranas"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))

using DataFrames
resultsUVap = collect_results(datadir("UV-scan", "anti_parallel", "final"))
resultsUhap = collect_results(datadir("Uh-scan", "anti_parallel", "final"))
resultUV = combine_results(resultsUVap)
resultUh = combine_results(resultsUhap)

##
bigres = 200
smallres = 100
# data = resultsUh
# pos = Int.(round.(size(data[:sweet_spots]) .* (20 / 40, 9 / 40)))
# ssparams = data[:sweet_spots][pos...].parameters
# ssparams = @set ssparams.V = .2
ssparams = (Δ=1, tratio=0.2, h=1.25, U=2.5, V=0.1, t=0.5, ϕ=2.083212845624367, μ1=3.5010504669947298, μ2=-0.801050470769927)
paramstring = map(x -> round(x, digits=2), ssparams)
μ0 = (ssparams.μ1 + ssparams.μ2) / 2
csdata_big = charge_stability_scan((; ssparams..., μ1=μ0, μ2=μ0), 10, 10, bigres);
ϕs = ssparams.ϕ .+ [-pi / 2, 0, pi / 2]
csdata_small = [charge_stability_scan((; ssparams..., ϕ), 1.5, 1.5, smallres) for ϕ in ϕs]

## 
labs = collect(keys(csdata_big[:data][1].majcoeffs[1]))
for lab in labs
    heatmap(map(x -> 1 - x.vacuumnorms.odd[lab], csdata_big[:data]),
        colorrange=(0, 1), axis=(; subtitle="$lab in odd vacuum")) |> display
end
for lab in labs
    heatmap(map(x -> 1 - x.vacuumnorms.even[lab], csdata_big[:data]),
        colorrange=(0, 1), axis=(; subtitle="$lab in even vacuum")) |> display
end
## Horizontal fig
f = Figure(resolution=400 .* (1, 0.9), fontsize=20, backgroundcolor=:transparent);
f = Figure(resolution=400 .* (1, 0.9), fontsize=20, backgroundcolor=:white);
g = f[1, 1] = GridLayout()
gb = g[1, 1] = GridLayout()
gs = g[1, 2] = GridLayout()


datamap = x -> sign(x.gap)
ax = Axis(gb[1, 1], xlabel=paramstyle[:μ1], ylabel=paramstyle[:μ2], aspect=1)#,  subtitle = L"\tanh{\left(δE\right)}")
hm = plot_charge_stability!(ax, csdata_big; colorrange=(-1, 1), datamap=x -> tanh(x.gap), colormap=Reverse(:redsblues))
cb = Colorbar(gb[2, 1], hm; valign=:top, vertical=false, label=L"\tanh{(δE)}", flipaxis=true, labelpadding=-20, ticks=[-1, 1], height=12
)
rowsize!(gb, 2, Relative(0.15))


boxax = Axis(gb[1, 1], aspect=1)
hidedecorations!(boxax)
linewidth = 2
linecolor = :black
lines!(boxax, -ssparams.μ1 .+ [-1, 1], -ssparams.μ2 .+ [1, 1]; color=linecolor, linewidth)
lines!(boxax, -ssparams.μ1 .+ [-1, 1], -ssparams.μ2 .- [1, 1]; color=linecolor, linewidth)
lines!(boxax, -ssparams.μ1 .+ [-1, -1], -ssparams.μ2 .+ [-1, 1]; color=linecolor, linewidth)
lines!(boxax, -ssparams.μ1 .+ [1, 1], -ssparams.μ2 .+ [-1, 1]; color=linecolor, linewidth)
linkaxes!(boxax, ax)
xlims!(minimum(csdata_big[:ϵ1]), maximum(csdata_big[:ϵ1]))
ylims!(minimum(csdata_big[:ϵ2]), maximum(csdata_big[:ϵ2]))
# ax1 = Axis(gs[1,1], xlabel = L"μ_1", ylabel = "after", fontsize = 1)

spinecolor = NamedTuple(map(x -> x => linecolor, [:bottomspinecolor, :leftspinecolor, :topspinecolor, :rightspinecolor]))
small_axes = [Axis(gs[n, 1]; spinecolor..., spinewidth=linewidth, aspect=1) for (n, ϕ) in enumerate(ϕs)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
foreach((ax, data) -> plot_charge_stability!(ax, data; datamap), small_axes, csdata_small)
hidedecorations!.(small_axes)
ax.xticks = -6:3:10
ax.yticks = -6:3:10
ax.xticklabelspace = 15.0#tight_xticklabel_spacing!(ax)
ax.yticklabelspace = 5.0#tight_xticklabel_spacing!(ax)
# rowgap!(g, 5)
# labels = [Label(gs[n, 1, Bottom()], L"ϕ ≈ %$(round(ϕ,digits=2))", padding=(0, 0, -5, 2)) for (n, ϕ) in enumerate(ϕs)]
[Label(gs[n, 1, Bottom()], L"\delta ϕ %$s \delta ϕ_\star", padding=(0, 0, -5, 2)) for (n, s) in enumerate(["<", "=", ">"])]

# Label(gb[1, 1, Top()], "tanh(δE)", valign = :top,
# Label(gb[1, 1, Top()], L"\tanh{(δE)}", valign = :top,
# padding = (0, 0, 2, -10))
# Label(gs[1, 1, Top()], L"\text{Parity}", valign=:bottom,
#     padding=(0, 0, 3, -10))

elems = [PolyElement(color=cgrad(:berlin)[1], strokewidth=1),
    PolyElement(color=cgrad(:berlin)[end], strokewidth=1)]
Legend(gb[end, 1], elems,
    [" "^27, ""],
    # [L"\text{ Odd  \qquad \qquad}", L"\text{Even}"],
    framevisible=false, orientation=:horizontal,
    tellwidth=false, tellheight=false, padding=(5, 0, 0, 0),
    margin=(0, 0, -60, 0))

colsize!(g, 2, Relative(0.3))
# rowsize!(gb, 2, Auto(.1))
cb.alignmode = Mixed(top=-10, bottom=-10)
# rowgap!(gb, 10)
# trim!(f.layout)
Label(gb[1, 1, TopLeft()], L"a)", padding=(0, 25, 0, -5), tellheight=false, tellwidth=false)
Label(gs[1, 1, TopLeft()], L"b)", padding=(0, 15, 0, -5), tellheight=false, tellwidth=false)
Label(gb[2, 1, BottomLeft()], L"\text{Odd}", padding=(125, 0, -28, 0), tellheight=false, tellwidth=false)
Label(gb[2, 1, BottomRight()], L"\text{Even}", padding=(-95, 0, -28, 0), tellheight=false, tellwidth=false)

# ph = 0.3
# pw = 0.8
# pc = (-1, -1)
# poly!(ax, map(p -> p .+ pc, [(-pw, -ph), (-pw, ph), (pw, ph), (pw, -ph)]), color=:white, linewidth=linewidth)

# Label(gb[1, 1], "↓₁↓₂", color=:black, font = :bold, fontsize=20, padding=(0, 0, 0, 0), tellheight=false, tellwidth=false)
# Label(gb[1, 1], L"|↓_1↓_2 \rangle", color=:black, fontsize=20, padding=(0, 0, 0, 0), tellheight=false, tellwidth=false)

for (point, label) in zip(Base.product(range(-6.2, 1.3, 3), range(-6.2, 1.7, 3)),
    Base.product(["↓₁↑₁", "↓₁", ""], ["↓₂↑₂", "↓₂", ""]))
    # Base.product(["\\downarrow_1 \\uparrow_1", "\\downarrow_1", ""], ["\\downarrow_2 \\uparrow_2", "\\downarrow_2", ""]))
    # for (point, label, align) in zip(Base.product(range(-6.2, 3.3, 3), range(-4, 1.5, 3)), 
    #     Base.product(["↓₁↑₁", "↓₁", ""], ["↓₂↑₂", "↓₂", ""]), 
    #     Base.product([:left, :left, :right], [:top, :center, :bottom]))
    # for (point, label) in zip(Base.product(range(-6.2, 1.3, 3), range(-5, 2.5, 3)), Base.product(["↓₁↑₁", "↓₁", ""], ["↓₂↑₂", "↓₂", ""]))
    # println(align)
    # println(label[1]*label[2])
    text!(ax, label[1] * "\n" * label[2], font=:bold, color=:white, position=point, align=(:left, :bottom))
end
f |> display
##
save(plotsdir(string("horizontal_charge_stability_phase", paramstring, ".pdf")), f, pt_per_unit=1)
save(plotsdir(string("horizontal_charge_stability_phase", paramstring, ".png")), f, px_per_unit=4)
##
save(plotsdir(string("horizontal_charge_stability_phase_transparent", paramstring, ".png")), f, px_per_unit=4)
save(plotsdir(string("horizontal_charge_stability_phase_transparent", paramstring, ".pdf")), f, pt_per_unit=1)
