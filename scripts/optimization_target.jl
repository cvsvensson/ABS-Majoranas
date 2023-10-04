# This script contains code to compare different optimization targets for the sweet spot. MP, MPU and LD are compared.

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
fcomparison = Figure();
g = fcomparison[1, 1] = GridLayout()
for (n, target) in enumerate(targets)
    for (m, comp) in enumerate(targets)
        ss = sweet_spots[n]
        heatmap(g[n, m], Us, hs, map(comp, ss);
            axis=(; subtitle="Target: $(string(target)), plotted: $(string(comp))"),
            colorrange=(10.0^(-1.5), 1), colorscale=log10,
            colormap=Reverse(:viridis))
    end
end
fcomparison #MPU and LD gives similar sweet spots, MP stands out a little.
## 
fasymmetry = Figure();
g = fasymmetry[1, 1] = GridLayout()
asymmetry(x) = x.mps.left.mp - x.mps.right.mp
for (n, target) in enumerate(targets)
    heatmap(g[1, n], map(abs ∘ asymmetry, sweet_spots[n]), colorrange=(0, 1), colormap=:viridis)
end
fasymmetry #Optimizing for MP gives more asymmetric sweet spots

##
fparams = Figure();
gparams = fparams[1, 1] = GridLayout()
for (n, target) in enumerate(targets)
    for (m, comp) in enumerate([:μ1, :μ2, :ϕ])
        ss = sweet_spots[n]
        heatmap(gparams[n, m], Us, hs, map(x -> x.parameters[comp], ss);
            axis=(; subtitle="Target: $(string(target)), plotted: $(string(comp))"),
            colorrange=(-3, 3), colormap=:batlow)
    end
end
fparams
