using DrWatson
@quickactivate "ABS-Majoranas"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))

transport = Transport(QuantumDots.Pauli(), (; T=1 / 40, μ=(0.0, 0.0)))
## Example
fixedparams = (; U=2, V=0.01, Δ=1, tratio=0.2, t=0.5)
opt = Optimizer(
    hamfunc=(t, ϕ, μ1, μ2, h) -> abs_hamiltonian(c; μ1, μ2, t, ϕ, h, fixedparams...),
    ranges=[(0.1, 10.0) .* first(fixedparams.Δ), (0.0, 2.0π), (0.0, 20.0) .* first(fixedparams.Δ), (-20.0, 0.0) .* first(fixedparams.Δ), (1.0, 20.0) .* first(fixedparams.Δ)],
    initials=[1, π, first(fixedparams.Δ), -first(fixedparams.Δ), fixedparams.tratio^-1 * first(fixedparams.Δ)];
    MaxTime=5, minexcgap=first(fixedparams.Δ) / 4,
    exps=collect(range(0.1, 4, length=4)),
    tracemode=:compact,
    target=LD,
    ϵ=0.01, PopulationSize=50)


# hamfunc=(t, ϕ, μ1, μ2, h) -> (m = Matrix(abs_hamiltonian_odd(c; μ1, μ2, t, ϕ, h, fixedparams...)); [vec(real(m)) vec(imag(m))])
ss1 = get_sweet_spot(opt)
ss2 = get_sweet_spot(opt)
ss3 = get_sweet_spot(opt)
ss4 = get_sweet_spot(opt)
ssb = get_sweet_spot_borg(opt)
optsol1 = solve(opt.hamfunc(best_candidate(ss1)...); transport)
optsol2 = solve(opt.hamfunc(best_candidate(ss2)...); transport)
optsol3 = solve(opt.hamfunc(best_candidate(ss3)...); transport)
optsol4 = solve(opt.hamfunc(best_candidate(ss4)...); transport)
optsolb = solve(opt.hamfunc(best_candidate(ssb)...); transport)
optsols = [optsol1, optsol2, optsol3, optsol4, optsolb]
map(MPU, optsols)
map(LD, optsols)
map(x -> x.gap, optsols)
optsol1.mps
optsolb.mps
optsol2.mps
optsol3.mps
optsol.reduced.cells
optsol2.reduced.cells
optsol.gap
optsol2.gap
optsol3.gap

pf = pareto_frontier(ss2)
best_obj1, idx_obj1 = findmin(map(elm -> fitness(elm)[1], pf))
best_obj2, idx_obj2 = findmin(map(elm -> fitness(elm)[2], pf))
bo1_solution = params(pf[idx_obj1]) # get the solution candidate itself...
bo2_solution = params(pf[idx_obj2]) # get the solution candidate itself...

csdata = charge_stability_scan(merge(fixedparams, NamedTuple(zip((:t, :ϕ, :μ1, :μ2, :h), best_candidate(ss)))), 5, 5)
csdata2 = charge_stability_scan(merge(fixedparams, NamedTuple(zip((:t, :ϕ, :μ1, :μ2, :h), best_candidate(ss2)))), 5, 5)
plot_charge_stability(csdata)[1]
plot_charge_stability(csdata2)[1]

##
exps = collect(range(0.5, 3, 4))
PopulationSize = 50
ϵ = 0.05
target = MPU
MaxTime = 2

fixedparams = (; Δ=1.0, tratio=0.2, h=1.5, U=3, V=0.1, t=0.5)
hamfunc(ϕ, μ1, μ2) = abs_hamiltonian(c; μ1, μ2, ϕ, fixedparams...)
extra_cost((ϕ, μ1, μ2), e) = exp(-(e * abs(μ1 - μ2) + 1)^4) + e * (μ2 > μ1) * (μ2 - μ1)
function cf(args)
    sol = solve(hamfunc(args...))
    exp = 3
    cost_function(sol.energies, MPU(sol); exp, minexcgap = 0) + extra_cost(args, exp)
end
co = compare_optimizers(cf;
        SearchRange = [(0.0,1.0pi), (0,10), (-10,0)], NumDimensions=3,MaxTime=1, TraceInterval=10.0,
        TraceMode=:silent)
co

Methods = [:dxnes, :xnes, :de_rand_1_bin_radiuslimited, :adaptive_de_rand_1_bin_radiuslimited, :generating_set_search, :adaptive_de_rand_1_bin, :separable_nes, :probabilistic_descent, :resampling_memetic_search, :resampling_inheritance_memetic_search]
ss2 = [anti_parallel_sweet_spot(; fixedparams..., MaxTime, target, PopulationSize, exps, ϵ, Method) for Method in Methods]

map(MPU, ss)
map(LD, ss)
map(x -> x.gap, ss)
barplot(eachindex(Methods), map(MPU, ss); axis)
barplot(eachindex(Methods), map(LD, ss); axis)
barplot(eachindex(Methods), map(x -> log(abs(x.gap)),ss); axis)
barplot(eachindex(Methods), map(MPU, ss2); axis)
barplot(eachindex(Methods), map(LD, ss2); axis)
barplot(eachindex(Methods), map(x -> log(abs(x.gap)),ss2); axis)
let ss = ss2
    axis = (; xticks=(eachindex(Methods), string.(Methods)), xticklabelrotation=pi / 4)
    barplot(eachindex(Methods), map(MPU, ss); axis) |> display
    barplot(eachindex(Methods), map(LD, ss); axis) |> display
    barplot(eachindex(Methods), map(x -> log(abs(x.gap)),ss); axis) |> display
end

##
pss = parallel_sweet_spot(; U=2, V=0.01, Δ=1, tratio=0.2, t=0.5, h=1.5, MaxTime=5, target=LD)
apss = anti_parallel_sweet_spot(; U=2, V=0.01, Δ=1, tratio=0.2, t=0.5, h=1.5, MaxTime=5, target=LD)
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


