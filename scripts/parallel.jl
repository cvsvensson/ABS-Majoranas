using DrWatson
@quickactivate "Majorana sweet spot"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))
using DataFrames
resultsUVp = collect_results(datadir("UV-scan", "parallel"))
# resultsUVap = collect_results(datadir("UV-scan", "anti_parallel"))
resultsUhp = collect_results(datadir("Uh-scan", "parallel"))
# resultsUhap = collect_results(datadir("Uh-scan", "anti_parallel"))

##
ssdata = resultsUhp[1, :]
nx, ny = size(ssdata[:sweet_spots])
smallres = 100
# positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.55, 0.75), (0.25, 0.5), (0.22, 0.1), (0.25, 0.1)])
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.15, 0.15), (0.05, 0.2), (0.2, 0.15)])
sweet_spots = [ssdata[:sweet_spots][pos...] for pos in positions]
paramstring = map(ss -> map(x -> round(x, digits=2), ss.parameters), sweet_spots)
transport = Transport(QuantumDots.Pauli(), (; T=1 / 20, μ=(0.0, 0.0)))
csdata_small = [charge_stability_scan((; ss.parameters...), 0.8, 0.8, smallres; transport) for ss in sweet_spots]

##
# colors = Makie.wong_colors()[[2, 4, 6]]
# colors = cgrad(:Dark2_7, categorical=true)[[2, 3, 4]]
colors = cgrad(:rainbow, categorical=true)[[1, 2, 5, 4][1:length(csdata_small)]]
fixedparamstring = map(x -> round(x, digits=2), ssdata[:fixedparams])
f = Figure(resolution=400 .* (1.4, 1), fontsize=20, backgroundcolor=:transparent);
f = Figure(resolution=400 .* (1.4, 1), fontsize=20, backgroundcolor=:white);
g = f[1, 1] = GridLayout()
gb = g[1, 1] = GridLayout()
gs = g[1, 2] = GridLayout()

xlabel = paramstyle[ssdata[:xlabel]]
ylabel = paramstyle[ssdata[:ylabel]]
ax = Axis(gb[1, 1]; xlabel, ylabel)
hm = plot_sweet_scan!(ax, ssdata)

cb = add_exp_colorbar!(gb[2, 1], hm; vertical=false, label=L"1-\mathrm{MP}", flipaxis=false, labelpadding=-5, height=12)

scatter!(ax, map(ss -> ss.parameters[ssdata.xlabel], sweet_spots), map(ss -> ss.parameters[ssdata.ylabel], sweet_spots),
    color=colors, markersize=20, marker=:xcross, strokewidth=2)

datamap = x -> sign(x.gap)
datamap2 = x -> real(x.conductance)[1, 2] * real(x.conductance)[1, 1]
datamap2 = x -> real(x.conductance)[1, 2] - real(x.conductance)[2, 1]
datamap2 = x -> real(x.conductance)[1, 2]
spinecolor = [NamedTuple(map(x -> x => colors[n], [:bottomspinecolor, :leftspinecolor, :topspinecolor, :rightspinecolor])) for n in eachindex(csdata_small)]
small_axes = [Axis(gs[n, 1]; spinecolor[n]..., spinewidth=4, aspect=1) for (n, data) in enumerate(csdata_small)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
small_axes2 = [Axis(gs[n, 2]; spinecolor[n]..., spinewidth=4, aspect=1) for (n, data) in enumerate(csdata_small)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
foreach((ax, data) -> plot_charge_stability!(ax, data; datamap), small_axes, csdata_small)
foreach((ax, data) -> plot_charge_stability!(ax, data; datamap=datamap2, colormap=:vik), small_axes2, csdata_small)
hidedecorations!.(small_axes)
hidedecorations!.(small_axes2)
labels = [Label(gs[n, 1:2, Bottom()], L"MP ≈ %$(round(1-MPU(ss),digits=2))", padding=(0, 0, -15, 4)) for (n, ss) in enumerate(sweet_spots)]
Label(gs[1, 1, Top()], L"\text{Parity}", padding=(0, 0, 2, -10))
Label(gs[1, 2, Top()], L"G_{\text{nl}}", padding=(0, 0, 2, -10))

# Label(gb[1, 1, Top()], "tanh(δE)", valign = :top,
# Label(gb[1, 1, Top()], L"\tanh{(δE)}", valign = :top,
# padding = (0, 0, 2, -10))

# Label(gs[1, 1, Top()], L"\text{Parity}", valign=:bottom,
# padding=(0, 0, 3, -10))

colsize!(g, 2, Auto(0.4))
colsize!(g, 2, Auto(0.9))
# rowsize!(gb, 2, Auto(.1))
cb.alignmode = Mixed(top=-15, bottom=0)
# rowgap!(gb, 10)
# trim!(f.layout)

f |> display