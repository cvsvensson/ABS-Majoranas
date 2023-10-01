using DrWatson
@quickactivate "Majorana sweet spot"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))
using DataFrames

resultsUVap = collect_results(datadir("UV-scan", "anti_parallel", "final"))
resultsUhap = collect_results(datadir("Uh-scan", "anti_parallel", "final"))
resultsUh1 = collect_results(datadir("Uh-scan", "anti_parallel", "apvsp"))
resultsUV1 = collect_results(datadir("UV-scan", "anti_parallel", "apvsp"))
resultsUh = combine_results(resultsUhap)
resultsUV = combine_results(resultsUVap)

##
ssdata = resultsUh
nx, ny = size(ssdata[:sweet_spots])
smallres = 100
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.4, 0.7), (0.25, 0.5), (0.25, 0.38), (0.25, 0.33)])
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.5, 0.8), (0.3, 0.4), (0.1, 0.1)])
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.65, 0.8), (0.3, 0.4)])
sweet_spots = [ssdata[:sweet_spots][pos...] for pos in positions]
paramstring = map(ss -> map(x -> round(x, digits=2), ss.parameters), sweet_spots)
# transport = Transport(QuantumDots.Pauli(), (; T=1 / 20, μ=(0.0, 0.0)))
transport = missing
csdata_small = [charge_stability_scan((; ss.parameters...), 1.5, 1.5, smallres; transport) for ss in sweet_spots]

##
fixedparamstring = map(x -> round(x, digits=2), ssdata[:fixedparams])
f = Figure(resolution=400 .* (1.4, 1), fontsize=20, backgroundcolor=:transparent);
f = Figure(resolution=400 .* (1.4, 1), fontsize=20, backgroundcolor=:white);
g = f[1, 1] = GridLayout()
gb = g[1, 1] = GridLayout()
gs = g[1, 2] = GridLayout()

xlabel = paramstyle[ssdata[:xlabel]]
ylabel = paramstyle[ssdata[:ylabel]]
ax = Axis(gb[1, 1]; xlabel, ylabel)
hm = plot_sweet_scan!(ax, ssdata; datamap = MP)

cb = add_exp_colorbar!(gb[2, 1], hm; vertical=false, label=L"1-\mathrm{MP}", flipaxis=false, labelpadding=-5, height=12)

scatter!(ax, map(ss -> ss.parameters[ssdata.xlabel], sweet_spots), map(ss -> ss.parameters[ssdata.ylabel], sweet_spots),
    color=ss_colors[1:length(csdata_small)], markersize=20, marker=:xcross, strokewidth=2)

datamap = x -> sign(x.gap)
spinecolor = [NamedTuple(map(x -> x => ss_colors[n], [:bottomspinecolor, :leftspinecolor, :topspinecolor, :rightspinecolor])) for n in eachindex(csdata_small)]
small_axes = [Axis(gs[n, 1]; spinecolor[n]..., spinewidth=4, aspect=1, xlabel=paramstyle[:μ1]) for (n, data) in enumerate(csdata_small)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 

foreach((n, ax, data) ->
        begin
            plot_charge_stability!(ax, data; datamap)
            scatter!(ax, -sweet_spots[n].parameters.μ1, -sweet_spots[n].parameters.μ2,
                color=ss_colors[n], markersize=20, marker=:xcross, strokewidth=2)
        end,
    eachindex(small_axes), small_axes, csdata_small)

elems = [PolyElement(color=cgrad(:berlin)[1], strokewidth=1),
    PolyElement(color=cgrad(:berlin)[end], strokewidth=1)]
Legend(gs[end+1, 1], elems, [L"\text{Odd}", L"\text{Even}"], L"\text{Parity}", framevisible=false)

hidedecorations!.(small_axes)
labels = [Label(gs[n, 1, Top()], L"\mathrm{MP} ≈ %$(round(1-MP(ss),digits=2))", padding=(0, 0, 10, 0), tellheight=false, tellwidth=false) for (n, ss) in enumerate(sweet_spots)]
Label(gb[1, 1, TopLeft()], L"a)", padding=(0, 50, 0, -5), tellheight=false, tellwidth=false)
Label(gs[1, 1, TopLeft()], L"b)", padding=(0, 35, 0, -5), tellheight=false, tellwidth=false)
Label(gs[2, 1, Left()], L"\varepsilon_1", rotation=pi / 2, padding=(0, 30, 0, 0), tellheight=false, tellwidth=false)
Label(gs[2, 1, Bottom()], L"\varepsilon_2", padding=(0, 0, 0, 10), tellheight=false, tellwidth=false)

colsize!(g, 2, Relative(0.28))
cb.alignmode = Mixed(top=-15, bottom=0)
f |> display

##
save(plotsdir(string("uhplot", fixedparamstring, ".png")), f, px_per_unit=4)
######
##
ssdata = resultsUV
nx, ny = size(ssdata[:sweet_spots])
smallres = 300
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.55, 0.75), (0.25, 0.5), (0.22, 0.1), (0.25, 0.1)])
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.1, 0.6), (0.4, 0.2), (0.6, 0.6)])
positions = map(x -> Int.(floor.((nx, ny) .* x)), [(0.5, 0.4), (0.4, 0.2)])
sweet_spots = [ssdata[:sweet_spots][pos...] for pos in positions]
# transport = Transport(QuantumDots.Pauli(), (; T=1 / 20, μ=(0.0, 0.0)))
transport = missing
paramstring = map(ss -> map(x -> round(x, digits=2), ss.parameters), sweet_spots)

csdata_small = [charge_stability_scan((; ss.parameters..., ϕ=ss.parameters.ϕ, μ1=ss.parameters.μ1, μ2=ss.parameters.μ2), 1.5, 1.5, smallres; transport) for ss in sweet_spots]

##
fixedparamstring = map(x -> round(x, digits=2), ssdata[:fixedparams])
f = Figure(resolution=400 .* (1.4, 1), fontsize=20, backgroundcolor=:white);
f = Figure(resolution=400 .* (1.4, 1), fontsize=20, backgroundcolor=:transparent);
g = f[1, 1] = GridLayout()
gb = g[1, 1] = GridLayout()
gs = g[1, 2] = GridLayout()

xlabel = paramstyle[ssdata[:xlabel]]
ylabel = paramstyle[ssdata[:ylabel]]
ax = Axis(gb[1, 1]; xlabel, ylabel)
hm = plot_sweet_scan!(ax, ssdata; datamap = MP)

cb = add_exp_colorbar!(gb[2, 1], hm; vertical=false, label=L"1-\mathrm{MP}", flipaxis=false, labelpadding=-5, height=12)

scatter!(ax, map(ss -> ss.parameters[ssdata.xlabel], sweet_spots), map(ss -> ss.parameters[ssdata.ylabel], sweet_spots),
    color=ss_colors[1:length(csdata_small)], markersize=20, marker=:xcross, strokewidth=2)

datamap = x -> sign(x.gap)
spinecolor = [NamedTuple(map(x -> x => ss_colors[n], [:bottomspinecolor, :leftspinecolor, :topspinecolor, :rightspinecolor])) for n in eachindex(csdata_small)]
small_axes = [Axis(gs[n, 1]; spinecolor[n]..., spinewidth=4, aspect=1, xlabel=paramstyle[:μ1]) for (n, data) in enumerate(csdata_small)] #subtitle = L"ϕ=%$(round(ϕ,digits=2))" 
foreach((n, ax, data) ->
        begin
            plot_charge_stability!(ax, data; datamap)
            scatter!(ax, -sweet_spots[n].parameters.μ1, -sweet_spots[n].parameters.μ2,
                color=ss_colors[n], markersize=20, marker=:xcross, strokewidth=2)
        end,
    eachindex(small_axes), small_axes, csdata_small)

elems = [PolyElement(color=cgrad(:berlin)[1], strokewidth=1),
    PolyElement(color=cgrad(:berlin)[end], strokewidth=1)]
Legend(gs[end+1, 1], elems, [L"\text{Odd}", L"\text{Even}"], L"\text{Parity}", framevisible=false)

hidedecorations!.(small_axes)
labels = [Label(gs[n, 1, Top()], L"\mathrm{MP} ≈ %$(round(1-MP(ss),digits=2))", padding=(0, 0, 10, 0), tellheight=false, tellwidth=false) for (n, ss) in enumerate(sweet_spots)]
Label(gb[1, 1, TopLeft()], L"a)", padding=(0, 50, 0, -5), tellheight=false, tellwidth=false)
Label(gs[1, 1, TopLeft()], L"b)", padding=(0, 35, 0, -5), tellheight=false, tellwidth=false)
Label(gs[2, 1, Left()], L"\varepsilon_1", rotation=pi / 2, padding=(0, 30, 0, 0), tellheight=false, tellwidth=false)
Label(gs[2, 1, Bottom()], L"\varepsilon_2", padding=(0, 0, 0, 10), tellheight=false, tellwidth=false)

let Us = 0:1
    band!(gb[1, 1], Us, Us, maximum(Us), color=(:white, 0.5))
end

ax.yticks = 0:1

colsize!(g, 2, Relative(0.28))
cb.alignmode = Mixed(top=-15, bottom=0)

f |> display
##
save(plotsdir(string("uvplot", fixedparamstring, ".png")), f, px_per_unit=4)

##
data = resultsUV1[2,:]
fgap = Figure(; resolution=400 .* (1.4, 1), fontsize=20, backgroundcolor=:white);
grid = fgap[1, 1] = GridLayout()
ax, hm = heatmap(grid[1, 1], data[:x], data[:y], map(x -> log10(abs(x.gap)), data[:sweet_spots]),
    colorrange=(-10, -3), highclip=:red, colormap=Reverse(:viridis))
Colorbar(grid[2, 1], hm, label="log10(|δE|/Δ)", vertical=false)
ax2, hm2 = heatmap(grid[1, 2], data[:x], data[:y], map(x -> excgap(x.energies...), data[:sweet_spots]),
    colorrange=(0, 1), highclip=:red, colormap=Reverse(:viridis))
Colorbar(grid[2, 2], hm2, label="excgap", vertical=false)
fgap

##