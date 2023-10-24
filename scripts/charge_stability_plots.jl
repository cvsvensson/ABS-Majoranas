#This script produces the charge stability plot for the paper.
using DrWatson
@quickactivate "ABS-Majoranas"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))

##
bigres = 300
smallres = 200
ssparams = (Δ=1, tratio=0.2, h=1.25, U=2.5, V=0.1, t=0.5, ϕ=2.083212845624367, μ1=3.5010504669947298, μ2=-0.801050470769927)
paramstring = map(x -> round(x, digits=2), ssparams)
μ0 = (ssparams.μ1 + ssparams.μ2) / 2
csdata_big = charge_stability_scan((; ssparams..., μ1=μ0, μ2=μ0), 10, 10, bigres);
ϕs = ssparams.ϕ .+ [-pi / 2, 0, pi / 2]
csdata_small = [charge_stability_scan((; ssparams..., ϕ), 1.5, 1.5, smallres) for ϕ in ϕs]


## Horizontal fig
f = Figure(resolution=400 .* (1, 0.9), fontsize=20, backgroundcolor=:transparent);
f = Figure(resolution=400 .* (1, 0.9), fontsize=20, backgroundcolor=:white);
g = f[1, 1] = GridLayout()
gb = g[1, 1] = GridLayout()
gs = g[1, 2] = GridLayout()


datamap = x -> sign(x.gap)
ax = Axis(gb[1, 1], xlabel=paramstyle[:μ1], ylabel=paramstyle[:μ2], aspect=1)
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

spinecolor = NamedTuple(map(x -> x => linecolor, [:bottomspinecolor, :leftspinecolor, :topspinecolor, :rightspinecolor]))
small_axes = [Axis(gs[n, 1]; spinecolor..., spinewidth=linewidth, aspect=1) for (n, ϕ) in enumerate(ϕs)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
foreach((ax, data) -> plot_charge_stability!(ax, data; datamap), small_axes, csdata_small)
hidedecorations!.(small_axes)
ax.xticks = -6:3:10
ax.yticks = -6:3:10
ax.xticklabelspace = 15.0
ax.yticklabelspace = 5.0
[Label(gs[n, 1, Bottom()], L"\delta ϕ %$s \delta ϕ_\star", padding=(0, 0, -5, 2)) for (n, s) in enumerate(["<", "=", ">"])]
# [Label(gs[n, 1, Bottom()], s, padding=(0, 0, -5, 2)) for (n, s) in enumerate([L"\delta\phi_\star-\pi/2", L"\delta\phi_\star", L"\delta\phi_\star+\pi/2"])]


elems = [PolyElement(color=cgrad(:berlin)[1], strokewidth=1),
    PolyElement(color=cgrad(:berlin)[end], strokewidth=1)]
Legend(gb[end, 1], elems,
    [" "^27, ""],
    framevisible=false, orientation=:horizontal,
    tellwidth=false, tellheight=false, padding=(5, 0, 0, 0),
    margin=(0, 0, -60, 0))

colsize!(g, 2, Relative(0.3))
cb.alignmode = Mixed(top=-10, bottom=-10)

Label(gb[1, 1, TopLeft()], L"a)", padding=(0, 25, 0, -5), tellheight=false, tellwidth=false)
Label(gs[1, 1, TopLeft()], L"b)", padding=(0, 15, 0, -5), tellheight=false, tellwidth=false)
Label(gb[2, 1, BottomLeft()], L"\text{Odd}", padding=(125, 0, -28, 0), tellheight=false, tellwidth=false)
Label(gb[2, 1, BottomRight()], L"\text{Even}", padding=(-95, 0, -28, 0), tellheight=false, tellwidth=false)

function spinlabel(l1, l2)
    string("(", l1, ",", l2, ")")
end

for (point, label) in zip(Base.product(range(-4.9, 2.4, 3), range(-5.8, 2.2, 3)),
    Base.product(["↓↑", "↓", "0"], ["↓↑", "↓", "0"]))
    text!(ax, spinlabel(label...), fontsize=15, font=:bold, color=:yellow, position=point, align=(:center, :bottom))
end
f |> display
##
save(plotsdir(string("charge_stability_phase", paramstring, ".png")), f, px_per_unit=4)



## Critical current
ϕs = ssparams.ϕ .+ range(-pi, pi, 100)
dsols = [dsolve(p -> abs_hamiltonian(c; p...), merge(ssparams, (; ϕ)); dp=(; ϕ=1e-6)) for ϕ in ϕs]
ssparams.ϕ
sols = [solve(abs_hamiltonian(c; ssparams..., ϕ)) for ϕ in ϕs]
plot(ϕs, map(LD, sols))
plot(ϕs, map(x -> x.gap, sols))
energies = reduce(hcat, map(x -> [x.energies[1][1:2]..., x.energies[2][1:2]...], sols))
labels = ["o1", "o2", "e1", "e2"]
fig, ax, sp = series(ϕs, energies; labels);
axislegend(ax);
fig 

fig2, ax2, sp2 = series(ϕs[1:end-1], diff(energies; dims=2); labels = "d" .* labels, axis=(; limits=(nothing, 0.01 .* (-1, 1))));
axislegend(ax2);
fig2 

fig3, ax3, sp3 = series(ϕs[1:end-2], diff(diff(energies; dims=2); dims=2); labels = "d²" .* labels, axis=(; limits=(nothing, 0.001 .* (-1, 1))));
axislegend(ax3);
fig3

plot(ϕs[1:end-1], map(x -> x.dsol.energies[1][1], dsols) |> diff)
plot(ϕs[1:end-1], map(x -> x.dsol.energies[1][2], dsols) |> diff)
