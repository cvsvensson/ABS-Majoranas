using DrWatson
@quickactivate "ABS-Majoranas"
include(srcdir("abs_chain_misc.jl"))

fixedparams = (; Δ=1.0, t=0.5, tratio=0.2, V=0)
MaxTime = 2 # Maximum time in seconds for each sweet spot optimization 
target = MPU # Choose LD, MP or MPU
res = 2 #Resolution 
us = range(0, 5, length=res)
hs = range(0, 6, length=res)
exps = collect(range(0.5, 3, 3))
Methods = [:generating_set_search, :probabilistic_descent, :adaptive_de_rand_1_bin_radiuslimited]
for Method in Methods
    results = sweet_spot_scan((us, :U), (hs, :h), anti_parallel_sweet_spot; fixedparams, MaxTime, target, exps, Method)
    tagsave(datadir("Uh-scan", "anti_parallel", string(Method), savename(results, "jld2"; allowedtypes=(Number, String, NamedTuple))), results)
end
#results = sweet_spot_scan((us, :U), (hs, :h), parallel_sweet_spot; fixedparams, MaxTime, target);
#tagsave(datadir("Uh-scan","parallel", savename(results, "jld2"; allowedtypes=(Number, String, NamedTuple))), results)
