using DrWatson
@quickactivate "ABS-Majoranas"
include(srcdir("abs_chain_misc.jl"))

fixedparams = (; Δ=1.0, t=0.5, tratio=0.2, V=0) 
MaxTime = 10 # Maximum time in seconds for each sweet spot optimization 
cost = MPU # Choose LD, MP or MPU
res = 20 #Resolution 
us = range(0, 5, length=res)
hs = range(0, 6, length=res)
results = sweet_spot_scan((us, :U), (hs, :h), anti_parallel_sweet_spot; fixedparams, MaxTime, cost);
tagsave(datadir("Uh-scan","anti_parallel", savename(results, "jld2"; allowedtypes=(Number, String, NamedTuple))), results)

results = sweet_spot_scan((us, :U), (hs, :h), parallel_sweet_spot; fixedparams, MaxTime, cost);
tagsave(datadir("Uh-scan","parallel", savename(results, "jld2"; allowedtypes=(Number, String, NamedTuple))), results)
