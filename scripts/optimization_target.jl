using DrWatson
@quickactivate "ABS-Majoranas"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))

## Compare different optimzation targets
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
