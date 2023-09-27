using DrWatson
@quickactivate "ABS-Majoranas"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))

transport = Transport(QuantumDots.Pauli(), (; T=1 / 40, μ=(0.0, 0.0)))
## Example
fixedparams = (; U=(2.0, 3.0), V=0.01, Δ=(1.0, 2.0), tratio=0.2, t=0.5)
opt = Optimizer(
    hamfunc=(t, ϕ, μ1, μ2, h) -> abs_hamiltonian(c; μ1, μ2, t, ϕ, h, fixedparams...),
    ranges=[(0.1, 10.0) .* first(fixedparams.Δ), (0.0, 2.0π), (0.0, 20.0) .* first(fixedparams.Δ), (-20.0, 0.0) .* first(fixedparams.Δ), (1.0, 20.0) .* first(fixedparams.Δ)],
    initials=[1, π, first(fixedparams.Δ), -first(fixedparams.Δ), fixedparams.tratio^-1 * first(fixedparams.Δ)];
    MaxTime=5, minexcgap=first(fixedparams.Δ) / 4,
    exps=collect(range(0.1, 4, length=4)),
    tracemode=:compact,
    target=LD,
    ϵ=0.01, PopulationSize=50)

ss1 = get_sweet_spot(opt)
optsol1 = solve(opt.hamfunc(best_candidate(ss1)...); transport)
csdata = charge_stability_scan(merge(fixedparams, NamedTuple(zip((:t, :ϕ, :μ1, :μ2, :h), best_candidate(ss)))), 5, 5)
plot_charge_stability(csdata)[1]

##
targets = [MPU, MP, LD]
fixedparams = (; U=2.5, V=0.2, Δ=1, tratio=0.2, t=0.5, h=1.25)
optparams = (; MaxTime=5, TargetFitness=1e-6)
Us = range(0, 5, 20)
hs = range(0, 6, 20)
iter = collect(Base.product(Us, hs))
sweet_spots = [Folds.map(Uh -> anti_parallel_sweet_spot(; fixedparams..., U=Uh[1], h=Uh[2], optparams..., target), iter) for target in targets]
##
results = Dict(:sweet_spots => sweet_spots, :Us => Us, :hs => hs, :fixedparams => fixedparams, :optparams => optparams, :targets => targets)
tagsave(datadir("Uh-scan", "anti_parallel", "target_comparison", savename(results, "jld2"; allowedtypes=(Number, String, NamedTuple, Symbol))), results)

##
fcomparison = Figure();
g = fcomparison[1, 1] = GridLayout()
asymmetry(x) = abs(only(diff(x.reduced.cells)))
asymmetry(x) = abs(x.mps.left.mpu - x.mps.right.mpu)
for (n, target) in enumerate(targets)
    for (m, comp) in enumerate([targets..., asymmetry])
        ss = sweet_spots[n]
        heatmap(g[n, m], Us, hs, map(comp, ss);
            axis=(; subtitle="Target: $(string(target)), plotted: $(string(comp))"),
            colorrange=(10.0^(-1.5), 1), colorscale=log10,
            # colorrange=(0, 1),
            colormap=Reverse(:viridis))
    end
end
fcomparison
##
compmap = x -> x.mps.left.mp - x.mps.right.mp
heatmap(map(compmap, sweet_spots[1]) .-
        map(compmap, sweet_spots[2]), colorrange=(-2, 2), colormap = :redsblues)
heatmap(map(compmap, sweet_spots[1]) .-
        map(compmap, sweet_spots[3]), colorrange=(-2, 2), colormap = :redsblues)

##
fparams = Figure();
gparams = fparams[1, 1] = GridLayout()
for (n, target) in enumerate(targets)
    for (m, comp) in enumerate([:μ1, :μ2, :ϕ])
        ss = sweet_spots[n]
        heatmap(gparams[n, m], Us, hs, map(x -> x.parameters[comp], ss);
            axis=(; subtitle="Target: $(string(target)), plotted: $(string(comp))"),
            colorrange=(-3, 3), colormap=:turbo)
    end
end
fparams
##     
map((ss, target) -> heatmap(Us, hs, map(MPU, ss); axis=(; subtitle=string(target))), sweet_spots, targets)

pss = parallel_sweet_spot(; fixedparams, optparams, target=MPU)
apss = anti_parallel_sweet_spot(; fixedparams, optparams, target=MPU)
apss.parameters
LD(pss) #bad
MP(pss) #bad
LD(apss) #Small is good
MP(apss) #Small is good

##
@time csdata = charge_stability_scan((; fixedparams..., ϕ=0.6, U=0, V=0, μ1=0, μ2=0, h=1.5), 8, 8, 100);
@time csdata = charge_stability_scan((; fixedparams..., ϕ=0.6, U=0, V=0, μ1=0, μ2=0, h=1.5), 8, 8, 100; transport);
csfig, ax, hm = plot_charge_stability(csdata; colorrange=0.1)
display(csfig)

#Plot the non-local conductance
nlcondfig, _, _ = plot_charge_stability(csdata; datamap=x -> real(x.conductance[1, 2]), colormap=:vik, colorrange=2)
display(nlcondfig)


##Zoomed in charge-stability
@time csdata = charge_stability_scan((; fixedparams..., ϕ=2.9, U=4, V=0, μ1=-1.1, μ2=1.1 + 3.5, h=0.5), 3, 3, 100; transport);
csfig, ax, hm = plot_charge_stability(csdata; colorrange=0.05)
display(csfig)
nlcondfig, _, _ = plot_charge_stability(csdata; datamap=x -> real(x.conductance[1, 2]), colormap=:vik, colorrange=2)
display(nlcondfig)
parityfig, _, _ = plot_charge_stability(csdata; datamap=x -> sign(x.gap), colormap=:vik, colorrange=2)
display(parityfig)


