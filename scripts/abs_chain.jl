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
fixedparamsap = (; U=2.5, V=0.2, Δ=1, tratio=0.2, t=0.5, h=1.25)
fixedparamsp = (; U=2.5, V=0.2, Δ=1, tratio=1/0.2, t=0.5, h=1.25)
optparams = (; MaxTime=5, TargetFitness=1e-6)
pss = parallel_sweet_spot(; fixedparamsp..., optparams..., target=MPU)
apss = anti_parallel_sweet_spot(; fixedparamsap..., optparams..., target=MPU)
pss2 = parallel_sweet_spot(; fixedparamsap..., optparams..., target=MPU)
apss2 = anti_parallel_sweet_spot(; fixedparamsp..., optparams..., target=MPU)
apss.parameters
LD(pss) #bad
LD(apss) #Small is good
LD(pss2) #bad
LD(apss2) #Small is good

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


