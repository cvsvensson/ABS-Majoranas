using DrWatson
@quickactivate "ABS-Majoranas"
include(srcdir("abs_chain_misc.jl"))

V = parse(Float64, ARGS[1])
fixedparams = (; Δ=1.0, t=0.5, tratio=0.2, V) 
MaxTime = 30 # Maximum time in seconds for each sweet spot optimization 
target = MPU # Choose LD, MP or MPU
res = 40 #Resolution 
us = range(0, 5, length=res)
hs = range(0, 6, length=res)
exps = collect(range(0.5, 4, length=4))
Methods = [:generating_set_search, :probabilistic_descent, :adaptive_de_rand_1_bin_radiuslimited]

for Method in Methods
    results = sweet_spot_scan((us, :U), (hs, :h), anti_parallel_sweet_spot; fixedparams, MaxTime, target, exps, Method)
    tagsave(datadir("Uh-scan", "anti_parallel", string(Method), savename(results, "jld2"; allowedtypes=(Number, String, NamedTuple))), results)
end
#results = sweet_spot_scan((us, :U), (hs, :h), parallel_sweet_spot; fixedparams, MaxTime, target);
#tagsave(datadir("Uh-scan","parallel", savename(results, "jld2"; allowedtypes=(Number, String, NamedTuple))), results)

