#This script calculates the sweet spots as a function of local interactions U and non-local interactions V.
using DrWatson
@quickactivate "ABS-Majoranas"
include(srcdir("abs_chain_misc.jl"))

fixedparams = (; Δ=1.0, t=0.5, tratio=0.2, h=1.25)
MaxTime = 5 # Maximum time in seconds for each sweet spot optimization 
target = MPU # Choose LD, MP or MPU
res = 20 #Resolution 
us = range(0, 5, length=res)
vs = range(0, 1, length=res)
exps = collect(range(0.5, 3, length=4))
Methods = [:probabilistic_descent, :probabilistic_descent, :probabilistic_descent]

for (n, Method) in enumerate(Methods)
    results = sweet_spot_scan((us, :U), (vs, :V), anti_parallel_sweet_spot; fixedparams, MaxTime, target, Method)
    tagsave(datadir("UV-scan", "anti_parallel", "final", savename("$n", results, "jld2"; allowedtypes=(Number, String, NamedTuple, Symbol))), results)
end
