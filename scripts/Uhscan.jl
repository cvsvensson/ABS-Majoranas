using DrWatson
@quickactivate "ABS-Majoranas"
include(srcdir("abs_chain_misc.jl"))

V = parse(Float64, ARGS[1])
fixedparams = (; Δ=1.0, t=0.5, tratio=0.2, V) 
MaxTime = 200 # Maximum time in seconds for each sweet spot optimization 
target = MPU # Choose LD, MP or MPU
res = 100 #Resolution 
us = range(0, 5, length=res)
hs = range(0, 6, length=res)
exps = collect(range(0.5, 6, length=8))

results = sweet_spot_scan((us, :U), (hs, :h), anti_parallel_sweet_spot; exps, fixedparams, MaxTime, target);
tagsave(datadir("Uh-scan","anti_parallel", savename(results, "jld2"; allowedtypes=(Number, String, NamedTuple))), results)

#results = sweet_spot_scan((us, :U), (hs, :h), parallel_sweet_spot; fixedparams, MaxTime, #target);
#tagsave(datadir("Uh-scan","parallel", savename(results, "jld2"; allowedtypes=(Number, String, #NamedTuple))), results)
