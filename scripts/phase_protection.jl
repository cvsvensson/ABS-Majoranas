# Here, we check if the energy degeneracy is protected from phase fluctuations
using DrWatson
@quickactivate "Majorana sweet spot"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))
using DataFrames
using FiniteDifferences

resultsUVap = collect_results(datadir("UV-scan", "anti_parallel", "final"))
resultsUhap = collect_results(datadir("Uh-scan", "anti_parallel", "final"))
resultsUh = combine_results(resultsUhap)
resultsUV = combine_results(resultsUVap)

function phase_slope(ss)
    get_gap(ϕ) = solve(ss.optimization.hamfunc(ϕ, ss.parameters.μ1, ss.parameters.μ2)).gap
    return central_fdm(5, 1)(get_gap, ss.parameters.ϕ)
end
##
@time phase_slope(resultsUh.sweet_spots[1])
Uhslopes = map(phase_slope, resultsUh.sweet_spots)
##
nh = 40
fig, ax, s1 = plot((abs.(Uhslopes[:, nh])))
s2 = plot!(ax, map(Base.Fix1(*, 1) ∘ sqrt ∘ LD, (resultsUh.sweet_spots[:, nh])), color=:blue)
axislegend(ax,
    [s1, s2],
    ["slope", "LD"])
fig
##
nU = 25
fig = Figure()
g = fig[1, 1] = GridLayout()
gs = g[1, 1] = GridLayout()
gb = g[1, 2] = GridLayout()
ax1 = Axis(gs[1, 1]; xlabel = paramstyle[resultsUh.ylabel], xlabelsize = 25)
s1 = scatter!(ax1, resultsUh.y , abs.(Uhslopes[nU, :]);)
s2 = scatter!(ax1, resultsUh.y ,map(Base.Fix1(*, 1) ∘ sqrt ∘ LD, (resultsUh.sweet_spots[nU, :])))
axislegend(ax1,
    [s1, s2],
    ["∂ΔE/∂Δϕ", "LD"])
Label(gs[1,1,Top()], "U = $(round(resultsUh.x[nU],digits=3))", fontsize=30)
fig
##
ssdata = resultsUh
nx, ny = size(ssdata[:sweet_spots])
smallres = 300
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.4, 0.7), (0.25, 0.5), (0.25, 0.38), (0.25, 0.33)])
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.5, 0.8), (0.3, 0.4), (0.1, 0.1)])
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.65, 0.8), (0.3, 0.4)])
sweet_spots = [ssdata[:sweet_spots][pos...] for pos in positions]
paramstring = map(ss -> map(x -> round(x, digits=2), ss.parameters), sweet_spots)
transport = missing

##
fixedparamstring = map(x -> round(x, digits=2), ssdata[:fixedparams])

xlabel = paramstyle[ssdata[:xlabel]]
ylabel = paramstyle[ssdata[:ylabel]]
ax = Axis(gb[1, 1]; xlabel, ylabel)
hm = plot_sweet_scan!(ax, ssdata; datamap=MPU, colorrange=(10.0^(-2.3), 1))
minimum(MP, ssdata[:sweet_spots]) |> log10

cb = add_exp_colorbar!(gb[2, 1], hm; vertical=false, label=L"1-\mathrm{MP}", flipaxis=false, labelpadding=-5, height=12)

Label(gb[1, 1, TopLeft()], L"a)", padding=(0, 50, 0, -5), tellheight=false, tellwidth=false)

# colsize!(g, 2, Relative(0.28))
cb.alignmode = Mixed(top=-15, bottom=0)
lines!(ax, resultsUh.x[nU] .+ 1e-5 .*[-1,1], [0,6], color = :black, linewidth = 4, linestyle = :dash)
fig |> display
##
save(plotsdir(string("phase_protection.png")), fig, px_per_unit=4)
