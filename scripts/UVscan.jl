using DrWatson
@quickactivate "ABS-Majoranas"
include(srcdir("abs_chain_misc.jl"))

fixedparams = (; Δ=1.0, t=0.5, tratio=0.2, h=1.5) 
MaxTime = 10 # Maximum time in seconds for each sweet spot optimization 
cost = MPU # Choose LD, MP or MPU
res = 20 #Resolution 
us = range(0, 5, length=res)
vs = range(0, 1, length=res)

results = sweet_spot_scan((us, :U), (vs, :V), anti_parallel_sweet_spot; fixedparams, MaxTime, cost);
tagsave(datadir("UV-scan","anti_parallel", savename(results, "jld2"; allowedtypes=(Number, String, NamedTuple))), results)

results = sweet_spot_scan((us, :U), (vs, :V), parallel_sweet_spot; fixedparams, MaxTime, cost);
tagsave(datadir("UV-scan","parallel", savename(results, "jld2"; allowedtypes=(Number, String, NamedTuple))), results)
