using DrWatson
@quickactivate "Majorana sweet spot"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))
using DataFrames
resultsUVmethods = collect_results(datadir("UV-scan", "anti_parallel", "methods"))
resultsUhmethods = collect_results(datadir("Uh-scan", "anti_parallel", "methods"))

##
for n in 1:6
    scanfig = let data1 = resultsUhmethods[n, :], data2 = resultsUVmethods[n, :]
        fig = Figure(; resolution=400 .* (2, 1), fontsize=20, backgroundcolor=:white)
        xlabel1 = paramstyle[data1[:xlabel]]
        xlabel2 = paramstyle[data2[:xlabel]]
        ylabel1 = paramstyle[data1[:ylabel]]
        ylabel2 = paramstyle[data2[:ylabel]]
        ax1 = Axis(fig[1, 1]; xlabel=xlabel1, ylabel=ylabel1)
        ax2 = Axis(fig[1, 2]; xlabel=xlabel2, ylabel=ylabel2)
        g = fig[1, 1:3] = GridLayout()
        hm1 = plot_sweet_scan!(ax1, data1)#; datamap = x->abs(x.gap), colorrange = (1e-5,1))
        hm2 = plot_sweet_scan!(ax2, data2)#; datamap = x->abs(x.gap), colorrange = (1e-5,1))
        label = L"1-\text{MP}"
        add_exp_colorbar!(fig[1, 3], hm1;)
        paramstring1 = map(x -> round.(x, digits=2), data1[:fixedparams])
        paramstring2 = map(x -> round.(x, digits=2), data2[:fixedparams])

        supertitle = Label(fig[1, 3, Top()], label, fontsize=25, padding=(0, 0, 8, 0))
        supertitle2 = Label(fig[1, 1, Top()], string(data1[:sweet_spots][1].optimization.Method), fontsize=20, padding=(0, 0, 8, 0))
        supertitle2 = Label(fig[1, 2, Top()], string(data2[:sweet_spots][1].optimization.Method), fontsize=20, padding=(0, 0, 8, 0))
        colgap!(fig.layout, 1, 16)
        colgap!(fig.layout, 2, 10)
        fig |> display
    end
end

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