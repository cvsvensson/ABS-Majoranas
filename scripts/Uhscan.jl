using DrWatson
@quickactivate "ABS-Majoranas"
include(srcdir("abs_chain_misc.jl"))

V = parse(Float64, ARGS[1])
fixedparams = (; Δ=1.0, t=0.5, tratio=0.2, V) 
MaxTime = 20 # Maximum time in seconds for each sweet spot optimization 
target = MPU # Choose LD, MP or MPU
res = 20 #Resolution 
us = range(0, 5, length=res)
hs = range(0, 6, length=res)
exps = collect(range(0.5, 3, length=4))
Methods = [:probabilistic_descent, :probabilistic_descent, :generating_set_search]
#Methods = [:dxnes, :xnes, :de_rand_1_bin_radiuslimited, :adaptive_de_rand_1_bin_radiuslimited, :generating_set_search, :adaptive_de_rand_1_bin, :separable_nes, :probabilistic_descent, :resampling_memetic_search, :resampling_inheritance_memetic_search]

for (n,Method) in enumerate(Methods)
    results = sweet_spot_scan((us, :U), (hs, :h), anti_parallel_sweet_spot; fixedparams, MaxTime, target, exps, Method)
    tagsave(datadir("Uh-scan", "anti_parallel", "methods", savename("$n",results, "jld2"; allowedtypes=(Number, String, NamedTuple, Symbol))), results)
end
#results = sweet_spot_scan((us, :U), (hs, :h), parallel_sweet_spot; fixedparams, MaxTime, target);
#tagsave(datadir("Uh-scan","parallel", savename(results, "jld2"; allowedtypes=(Number, String, NamedTuple))), results)

