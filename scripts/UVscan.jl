using DrWatson
@quickactivate "ABS-Majoranas"
include(srcdir("abs_chain_misc.jl"))

h = parse(Float64, ARGS[1])
fixedparams = (; Δ=1.0, t=0.5, tratio=0.2, h) 
MaxTime = 60 # Maximum time in seconds for each sweet spot optimization 
target = MPU # Choose LD, MP or MPU
res = 40 #Resolution 
us = range(0, 5, length=res)
vs = range(0, 1, length=res)

results = sweet_spot_scan((us, :U), (vs, :V), anti_parallel_sweet_spot; fixedparams, MaxTime, target);
tagsave(datadir("UV-scan","anti_parallel", savename(results, "jld2"; allowedtypes=(Number, String, NamedTuple))), results)

#results = sweet_spot_scan((us, :U), (vs, :V), parallel_sweet_spot; fixedparams, MaxTime, #target);
#tagsave(datadir("UV-scan","parallel", savename(results, "jld2"; allowedtypes=(Number, String, #NamedTuple))), results)
