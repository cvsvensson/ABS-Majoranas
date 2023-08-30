using DrWatson
@quickactivate "ABS-Majoranas"
using QuantumDots, LinearAlgebra, QuantumDots.BlockDiagonals#, Base.Threads#, QuantumDots.BlockDiagonals
using BlackBoxOptim, LinearSolve
using Folds
using ForwardDiff
const N = 2
const c = FermionBasis((1, 2), (:↑, :↓); qn=QuantumDots.parity)
includet(srcdir("abs_chain_misc.jl"))

## Example
fixedparams = (; U=2, V=0.01, Δ=1, tratio=0.2, t=0.5)
opt = Optimizer(
    hamfunc=(t, ϕ, μ1, μ2, h) -> abs_hamiltonian(c; μ=[μ1, μ2], t, ϕ, h, fixedparams...),
    ranges=[(0.1, 10.0) .* fixedparams.Δ, (0.0, 2.0π), (0.0, 20.0) .* fixedparams.Δ, (-20.0, 0.0) .* fixedparams.Δ, (1.0, 20.0) .* fixedparams.Δ],
    initials=[1, π, -fixedparams.Δ, fixedparams.Δ, fixedparams.tratio^-1 * fixedparams.Δ];
    MaxTime=10, minexcgap=fixedparams.Δ / 4,
    exps=collect(range(1, 6, length=4)),
    tracemode=:compact,
    target=MPU())
ss = get_sweet_spot(opt)
optsol = solve(opt.hamfunc(ss...); transport=Transport(QuantumDots.Pauli(), (; T=1, μ=(0.1, 0.2), Γ=1)))
csdata = charge_stability_scan(merge(fixedparams, NamedTuple(zip((:t, :ϕ, :h), ss[[1,2,5]])), (;μ=ss[[3,4]])))

##

x = 0:5
y = range(0, 0.5, 5)
@time res = sweet_spot_scan((x, :U), (y, :V); fixedparams, MaxTime=2, target=LD());

data = charge_stability_scan(range(-5, 5, 40), range(-5, 5, 40); fixedparams..., ϕ=1.5, U=0, V=0);

##
function get_sweet_spot(; Δ, tratio, h, U, V, t, MaxTime, exps=collect(range(0.5, 3, length=4)), target=LD())
    fixedparams = (; Δ, tratio, h, U, V, t)
    pk = kitaev_sweet_spot2(; Δ, h, U, V, t, tratio)
    hamfunc(ϕ, μ1, μ2) = abs_hamiltonian(c; μ=[μ1, μ2], ϕ, fixedparams...)
    hamfunc(ϕ, μ) = abs_hamiltonian(c; μ, ϕ, fixedparams...)
    opt = Optimizer(;
        hamfunc,
        ranges=[(0.0, 1.0π), (0.0, 1.1 * pk.μ[1] + h + U), (-2h - U, U + V)],
        initials=Float64.([pk.ϕ, pk.μ...]),
        MaxTime, exps, target,
        refinefactor=(0.9, :hard_limit),
        tracemode=:silent,
        extra_cost=((θ, μ1, μ2), e) -> exp(-(e * abs(μ1 - μ2) + 1)^4))
    ss = get_sweet_spot(opt)
    optsol = solve(opt.hamfunc(ss...))
    parameters = merge(fixedparams, (; ϕ=ss[1], μ=(ss[2], ss[3])))
    optimization = opt
    sweet_spot = merge(optsol, (; parameters, optimization))
    return sweet_spot
end
##
MaxTime = 30
# non_int_sweet_spots = Vector{Any}(undef,N_ss)
minexcs = [0.2]#[0.05,.1,.2,.25,.3,.35,.4,.5,1,2]#range(.01,1,length=N_ss)
# N_ss = length(4
fixedparams = (; Δ=1, tratio=0.2, h=1.5, t=0.5)
ssLD = get_sweet_spot(; U=0, V=0, fixedparams..., MaxTime=10, target=LD())
dataLD = charge_stability_scan()
ssMP = get_sweet_spot(; U=0, V=0, fixedparams..., MaxTime, target=MP())
ssMPU = get_sweet_spot(; U=0, V=0, fixedparams..., MaxTime, target=MPU())

hs = range(fixedparams.Δ * 1.5, 0.6 * sqrt(fixedparams.tratio^-2 + 1) * fixedparams.Δ, length=3)
non_int_sweet_spots = stack(get_sweet_spot(; U=0, V=0, fixedparams.Δ, fixedparams.tratio, h, MaxTime, minexcs) for h in hs)
int_sweet_spots = [stack(get_sweet_spot(; U, V=0, fixedparams.Δ, fixedparams.tratio, h, MaxTime, minexcs) for h in hs) for U in range(0, 2, length=4)]
non_int_sweet_spots_mp = get_sweet_spot(; U=0, V=0, fixedparams..., MaxTime, minexcs, target=MP())
non_int_sweet_spots_mpu = get_sweet_spot(; U=0, V=0, fixedparams..., MaxTime, minexcs, target=MPU())
Uint_sweet_spots = get_sweet_spot(; U=1, V=0, fixedparams.Δ, fixedparams.tratio, h, MaxTime, minexcs)
Uint_sweet_spots_mpu = get_sweet_spot(; U=1, V=0, fixedparams..., MaxTime, minexcs, target=MPU())
Vint_sweet_spots = get_sweet_spot(; U=0, V=1 / 4, fixedparams..., MaxTime, minexcs)
UVint_sweet_spots = get_sweet_spot(; U=2, V=1 / 4, fixedparams..., MaxTime, minexcs)

Us = range(0, 5, length=6)
Us2 = range(4, 4.5, length=20)
Vs = range(0, 1, length=20)
us = range(0, 5, length=20)
vs = range(0, 1, length=20)
UVs = collect(Base.product(us, vs))
Vss = Folds.map(V -> get_sweet_spot(; U=0, V, fixedparams..., MaxTime=300, minexcs=[0.2])[1], Vs)
Uss = Folds.map(U -> get_sweet_spot(; U, V=0, fixedparams..., MaxTime=300, minexcs=[0.2])[1], Us)
Uss2 = Folds.map(U -> get_sweet_spot(; U, V=0, fixedparams..., MaxTime=60, minexcs=[0.2])[1], Us2)
Uss3 = Folds.map(U -> get_sweet_spot(; U, V=0, fixedparams..., MaxTime=60, minexcs=[0.2])[1], Us)
Uss4 = Folds.map(U -> get_sweet_spot(; U, V=0, fixedparams..., MaxTime=60, minexcs=[0.2])[1], Us)
UVss = Folds.map(UV -> get_sweet_spot(; U=UV[1], V=UV[2], fixedparams..., MaxTime=60, minexcs=[0.2])[1], UVs)

plot(Us, map(ss -> ss.mpu, Uss), marker=true, ylims=(0, 1))
plot(Us2, map(ss -> ss.mpu, Uss2), marker=true, ylims=(0, 1))
plot(Us2, map(Base.Fix1(quality_measure, LD()), Uss2), marker=true, ylims=(0, 1))
plot(Us2, map(Base.Fix1(quality_measure, MPU()), Uss2), marker=true, ylims=(0, 1))
plot(Vs, map(Base.Fix1(quality_measure, MPU()), Vss), marker=true, ylims=(0, 1))
for key in keys(first(Uss2).parameters)
    plot(Us2, map(ss -> ss.parameters[key], Uss2); marker=true, label=string(key), xlabel="U") |> display
end
for key in keys(first(Vss).parameters)
    plot(Vs, map(ss -> ss.parameters[key], Vss); marker=true, label=string(key), xlabel="V") |> display
end

heatmap(vs, us, map(log ∘ Base.Fix1(quality_measure, MP()), UVss), xlabel="V", ylabel="U", colorbar_title="log(1-MP)")
contourf(vs, us, map(log ∘ Base.Fix1(quality_measure, MPU()), UVss))

plot(Us2, map(ss -> ss.parameters.μ[1], Uss2), marker=true)
plot(Us2, stack(map(ss -> collect(ss.parameters), Uss2))', marker=true)
plot(Us, map(ss -> ss.mpu, Uss3), marker=true, ylims=(0, 1))
plot(Us, map(ss -> ss.mpu, Uss4), marker=true, ylims=(0, 1))
Us2 = range(3.2, 3.4, length=5)
Uss2 = [get_sweet_spot(; U, V=0, fixedparams..., MaxTime, minexcs=[0.2])[1] for U in Us2]
plot(Us2, map(ss -> ss.mpu, Uss2), marker=true)
Vs = range(0, 1, length=5)
Vss = [get_sweet_spot(; U=0, V, fixedparams..., MaxTime, minexcs=[0.05])[1] for V in Vs]
plot(Vs, map(ss -> ss.mpu, Vss), marker=true)
##
θs = range(0, π, length=50)
plot_mp(non_int_sweet_spots[2])
plot(θs, stack([get_mps(non_int_sweet_spots[2].majcoeffs, t) for t in θs])', ylims=(-1, 1), label=["mp" "mpu"])

##
res = 100
@time nonintssdata = map(ss -> get_param_scan_data(ss; μwidth=0.1, ΓoT=1 / 10, res), non_int_sweet_spots)
# function guess_params(Δ, tratio, h, U)
#     ϕ, μ1, μ2 = phase_μs(; tratio, h, Δ, U)
#     t = 1
#     (; Δ=Δ, tratio, h, U, ϕ, μ1, μ2, t, V=0)
# end
guessdata = map(U -> get_param_scan_data(guess_params(1, 0.2, 4, U); μwidth=0.1, res), range(0, 3, length=4))
nonintssdata_mpu = [get_param_scan_data(ss; res) for ss in non_int_sweet_spots_mpu]
Ussdata = [get_param_scan_data(ss; μwidth=0.1, res) for ss in Uss]
Ussdatalind = [get_param_scan_data(ss; μwidth=0.1, res=20, alg=:lindblad) for ss in Uss]
Ussdata2 = [get_param_scan_data(ss; μwidth=0.01, res) for ss in Uss2[[1, 25]]]
Ussdata_mpu = [get_param_scan_data(ss; res) for ss in Uint_sweet_spots_mpu]
UVssdata = [get_param_scan_data(ss; res) for ss in UVss]
Vssdata = [get_param_scan_data(ss; μwidth=0.01, res) for ss in Vss]
@profview nonintssdatalind3 = [get_param_scan_data(ss; res=10, alg=:lindblad, ΓoT=1 / 10) for ss in non_int_sweet_spots]
@time nonintssdatalind3 = [get_param_scan_data(ss; μwidth=0.1, res=20, alg=:lindblad, ΓoT=0.1) for ss in non_int_sweet_spots]
@time nonintssdatalind4 = [get_param_scan_data(ss; μwidth=0.1, res=2, alg=:lindblad, ΓoT=0.1) for ss in non_int_sweet_spots]
##
wsave(sdatadir("sims", "ABS-chain", string("data_", savename(fixedparams, "jld2"))),
    Dict(["0" => nonintssdata, "U" => Ussdata, "UV" => UVssdata, "V" => Vssdata]))
wsave(sdatadir("sims", "ABS-chain", string("sweet_spots_", savename(fixedparams, "jld2"))),
    Dict(["0" => non_int_sweet_spots, "U" => Uint_sweet_spots, "UV" => UVint_sweet_spots, "V" => Vint_sweet_spots]))
##
plotsize = 700 .* (1, 0.8)
# for (k,data) in enumerate(nonintssdata)
for (k, data) in enumerate(nonintssdatalind3)
    # for (k, data) in enumerate(guessdata)
    p = plot(param_scan_plots(data, target=MPU())..., layout=(3, 2), size=plotsize)
    p |> display
    #save(splotsdir("ABS-chain", string("nonint",k,".png")), p)
    # save(splotsdir("ABS-chain", string("nonint-high-temp",k,".png")), p)
    #save(splotsdir("ABS-chain", string("nonint-high-temp-lindblad2",k,".png")), p)
end
for (k, data) in enumerate(Ussdatalind)
    p = plot(param_scan_plots(data, target=MPU())..., layout=(3, 2), size=plotsize)
    p |> display
    # save(splotsdir("ABS-chain", string("U",k,".png")), p )
end
for (k, data) in enumerate(UVssdata)
    p = plot(param_scan_plots(data)..., layout=4, size=plotsize)
    p |> display
    # save(splotsdir("ABS-chain",string("UV",k,".png")), p)
end
for (k, data) in enumerate(Vssdata)
    p = plot(param_scan_plots(data)..., layout=(3, 2), size=plotsize)
    p |> display
    # save(splotsdir("ABS-chain", string("V",k,".png")), p)
end
##
let
    datas = [non_int_sweet_spots, Uint_sweet_spots, Vint_sweet_spots, UVint_sweet_spots]
    lds = [map(d -> norm(d.reduced.cells), data) for data in datas]
    excgaps = [map(d -> excgap(d.energies...), data) for data in datas]
    legends = ["(U,V) = ($(first(data).parameters.U),$(first(data).parameters.V))" for data in datas]
    p = plot(; xlims=(0, 0.1 + maximum(maximum.(excgaps))), ylims=(0, 0.1 + maximum(maximum.(lds))),
        xlabel="Excgap", ylabel="LD", framestyle=:box)
    foreach((excgap, ld, leg) -> plot!(p, excgap, ld; label=leg, markers=true), excgaps, lds, legends)
    save(splotsdir("ABS-chain", "scaling.png"), p)
    p
end
map(d -> excgap(d.energies...), non_int_sweet_spots)

##
Us = range(0, 5, 50)
Vs = range(0, 1, 50)
hs = [1.5, 3.0, 4.0]
results_guess = [Matrix{Any}(undef, length(Us), length(Vs)) for h in hs]
for (hn, h) in enumerate(hs)
    fp = @set fixedparams.h = h
    for (n, U) in enumerate(Us), (m, V) in enumerate(Vs)
        pk = kitaev_sweet_spot2(; U, V, fp...)
        parameters = (; U, V, μ=[pk.μ...], ϕ=pk.δϕ, fp...)
        sol = solve(abs_hamiltonian(c; parameters...))
        results_guess[hn][n, m] = merge(sol, (; parameters))
    end
end

##
using DataFrames
resultsUV = collect_results(sdatadir("sims", "ABS-chain", "UV-scan"))
resultsUh = collect_results(sdatadir("sims", "ABS-chain", "Uh-scan"))
resultsVh = collect_results(sdatadir("sims", "ABS-chain", "Vh-scan"))

for n in 1:3, (label, data) in [(:opt, results[:, :sweet_spots]), (:guess, results_guess)]
    display(heatmap(Vs, Us, map(ss -> mod(ss.parameters.ϕ, pi), data[n]), title="$n, $label, ϕ", clims=(0, pi)))
end
for n in 1:3, (label, data) in [(:opt, results[:, :sweet_spots]), (:guess, results_guess)]
    display(heatmap(Vs, Us, map(Base.Fix1(quality_measure, MP()), data[n]), title="$n, $label, (1-MP)", clims=(0, 1)))
end
for n in 1:3, (label, data) in [(:opt, results[:, :sweet_spots]), (:guess, results_guess)]
    display(heatmap(Vs, Us, map(ss -> ss.parameters.t, data[n]), title="$n, $label, t", clims=(0, 1)))
end
for n in 1:3, (label, data) in [(:opt, results[:, :sweet_spots]), (:guess, results_guess)]
    display(heatmap(Vs, Us, map(ss -> excgap(ss.energies...), data[n]), title="$n, $label, excgap", clims=(0, 1)))
end

let p = :μ
    for n in 1:3, (label, data) in [(:opt, results[:, :sweet_spots]), (:guess, results_guess)]
        display(heatmap(Vs, Us, map(ss -> ss.parameters.μs[1], data[n]), title="$n, $label, $p", clims=(0, 1)))
    end
end

us = map(ss -> ss.parameters.U, resultsUV[1, :sweet_spots][:, 1])
vs = map(ss -> ss.parameters.V, resultsUV[1, :sweet_spots][1, :])
hs = map(ss -> ss.parameters.h, resultsUh[1, :sweet_spots][1, :])
let data = resultsUh
    for n in axes(data, 1)
        heatmap(us, hs, map(ss -> mod(ss.parameters.ϕ, 2pi), data[n, :sweet_spots]), title="$(data[n,:params]) \n ϕ", clims=(0, 2pi))
        plot!(h -> sqrt(1 + data[n, :params].tratio^(-2)) - h; xlims=(first(us), last(us)), ylims=(first(hs), last(hs))) |> display
    end
end
let data = resultsUh
    for n in axes(data, 1)
        p = heatmap(hs, us, map(log ∘ Base.Fix1(quality_measure, LD()), data[n, :sweet_spots]);
            title="$(data[n,:params]) \n LD", clims=(-5, 0), xlabel="h", ylabel="U")
        # p = heatmap(hs,us,map(log ∘ Base.Fix1(quality_measure, LD()), permutedims(data[n, :sweet_spots],(2,1))), title="$(data[n,:params]) \n LD", clims=(-5, 0))
        bound(h) = sqrt(1 + data[n, :params].tratio^(-2)) - h
        plot!(p, bound; ylims=(first(us), last(us)), xlims=(first(hs), last(hs)), label="bound")
        p |> display
    end
end

let data = resultsUh
    for n in axes(data, 1)
        p = heatmap(hs, us, map(ss -> excgap(ss.energies...), data[n, :sweet_spots]);
            title="$(data[n,:params]) \n LD", xlabel="h", ylabel="U", clims=(0, 0.4))
        # p = heatmap(hs,us,map(log ∘ Base.Fix1(quality_measure, LD()), permutedims(data[n, :sweet_spots],(2,1))), title="$(data[n,:params]) \n LD", clims=(-5, 0))
        bound(h) = sqrt(1 + data[n, :params].tratio^(-2)) - h
        plot!(p, bound; ylims=(first(us), last(us)), xlims=(first(hs), last(hs)), label="bound")
        p |> display
    end
end

let data = resultsUh
    for n in axes(data, 1)
        p = heatmap(hs, us, map(ss -> ss.parameters.μ[2], data[n, :sweet_spots]);
            title="$(data[n,:params]) \n U", xlabel="h", ylabel="U")
        # p = heatmap(hs,us,map(log ∘ Base.Fix1(quality_measure, LD()), permutedims(data[n, :sweet_spots],(2,1))), title="$(data[n,:params]) \n LD", clims=(-5, 0))
        bound(U) = sqrt(1 + data[n, :params].tratio^(-2)) - U / 2
        plot!(p, bound; ylims=(first(us), last(us)), xlims=(first(hs), last(hs)))
        p |> display
    end
end

let data = resultsUV
    for n in axes(data, 1)
        p = heatmap(vs, us, map(log ∘ Base.Fix1(quality_measure, LD()), data[n, :sweet_spots]);
            title="$(data[n,:params]) \n LD", clims=(-5, 0), xlabel="V", ylabel="U")
        contour!(p, vs, us, map(ss -> KitaevParameters(; ss.parameters...) |> p -> -abs(p.V / p.t), data[n, :sweet_spots]), lw=3)

        # bound(h) = sqrt(1+data[n,:params].tratio^(-2)) - h + 1
        # plot!(p, bound; xlims = (first(us),last(us)), ylims = (first(hs),last(hs)))
        p |> display
    end
end
let data = resultsVh
    for n in axes(data, 1)
        p = heatmap(hs, vs, map(log ∘ Base.Fix1(quality_measure, LD()), data[n, :sweet_spots]);
            title="$(data[n,:params]) \n LD", clims=(-5, 0), xlabel="h", ylabel="V")
        contour!(p, hs, vs, map(ss -> KitaevParameters(; ss.parameters...) |> p -> -abs(p.V / p.t), data[n, :sweet_spots]), lw=3)
        # bound(h) = sqrt(1+data[n,:params].tratio^(-2)) - h + 1
        # plot!(p, bound; xlims = (first(us),last(us)), ylims = (first(hs),last(hs)))
        p |> display
    end
end


for n in 1:7
    data = results3
    display(heatmap(map(log ∘ Base.Fix1(quality_measure, LD()), data[n, :sweet_spots]), title="$(data[n,:params]) \n LD", clims=(-5, 0)))
end
for n in 1:7
    data = results3
    display(heatmap(map(ss -> excgap(ss.energies...), data[n, :sweet_spots]), title="$(data[n,:params]) \n excgap", clims=(0, 1)))
end


for n in 3:3
    data = results3
    display(heatmap(map(ss -> ss.parameters.μ1 - ss.parameters.U, data[n, :sweet_spots]), title="$(data[n,:params]) \n μ"))
end
##
heatmap(map(ss -> KitaevParameters(; ss.parameters...) |> p -> p.V / ss.parameters.V, resultsUV[3, :sweet_spots]),)
heatmap(map(ss -> KitaevParameters(; ss.parameters...) |> p -> (abs(p.t) - p.V / 2) / abs(p.Δ), resultsUV[4, :sweet_spots]), clims=(0.5, 1.5), c=:redsblues)
contourf(hs, vs, map(ss -> KitaevParameters(; ss.parameters...) |> p -> abs(p.V / p.Δ), resultsVh[2, :sweet_spots]))
heatmap(hs, vs, map(ss -> KitaevParameters(; ss.parameters...) |> p -> (p.μ[2])), resultsVh[2, :sweet_spots])
heatmap(hs, vs, map(ss -> KitaevParameters(; ss.parameters...) |> p -> (p.μ[1])), resultsVh[2, :sweet_spots])
##
current_resolution = 30
currents = Matrix{Any}(undef, current_resolution, current_resolution)
@time let U = 0.0, V = 0.0, Δ = 1.0, tratio = 0.2, h = 5, t = 0.3, ϕ = π
    T = Δ / 40
    Γ = T / 10
    fixedparams = (; U, V, Δ, tratio)
    @threads for (n1, V) in collect(enumerate(2h * range(-1, 1, length=current_resolution)))
        for (n2, μ) in (enumerate(2h * range(-1, 1, length=current_resolution)))
            currents[n1, n2] = current(abs_hamiltonian(c; μ=[5, -5] .+ μ, Δ, U, V=0, tratio, h, t, ϕ), leads(c, T, V, -V, Γ))
        end
    end
end
current_resolution2 = 100
currents2 = Matrix{Any}(undef, current_resolution2, current_resolution2)
@time let U = 0.0, V = 0.0, Δ = 1.0, tratio = 0.2, h = 5, t = 0.3, ϕ = π
    T = Δ / 40
    Γ = T / 10
    fixedparams = (; U, V, Δ, tratio)
    @threads for (n1, V) in collect(enumerate(2h * range(-1, 1, length=current_resolution2)))
        for (n2, μ) in (enumerate(2h * range(-1, 1, length=current_resolution2)))
            currents2[n1, n2] = current2(abs_hamiltonian(c; μ=[5, -5] .+ μ, Δ, U, V=0, tratio, h, t, ϕ), leads(c, T, V, -V, Γ))
        end
    end
end
map(sum, currents2)
heatmap(map(sum, currents2); c=:redsblues, clims=maximum(map(abs ∘ sum, currents2)) .* (-1, 1))
heatmap(map(c -> sum(c[1:2, :]), currents); c=:redsblues, clims=maximum(map(abs ∘ sum, currents2)) .* (-1, 1))
heatmap(map(c -> sum(c[1:2, :]), currents2); c=:redsblues, clims=maximum(map(abs ∘ sum, currents2)) .* (-1, 1))
heatmap(map(c -> sum(c[3:4, :]), currents2); c=:redsblues, clims=maximum(map(abs ∘ sum, currents2)) .* (-1, 1))
heatmap(map(c -> sum(c[1]), currents2); c=:redsblues, clims=maximum(map(abs ∘ sum, currents2)) .* (-1, 1))


##


##
conductance_resolution = 100
conductances = Matrix{Any}(undef, conductance_resolution, conductance_resolution)
@time let U = 0.0, V = 0.0, Δ = 1.0, tratio = 0.2, h = 5, t = 2, dV = Δ / 10, ϕ = π
    T = Δ / 40
    Γ = T / 10
    fixedparams = (; U, V, Δ, tratio)
    @threads for (n1, μ1) in collect(enumerate(-h .+ Δ * range(-1, 1, length=conductance_resolution)))
        for (n2, μ2) in (enumerate(h .+ Δ * range(-1, 1, length=conductance_resolution)))
            H = abs_hamiltonian(c; μ=[μ1, μ2], Δ, U, V, tratio, h, t, ϕ)
            c1 = current2(H, leads(c, T, dV, 0, Γ))
            c2 = current2(H, leads(c, T, -dV, 0, Γ))
            conductances[n1, n2] = (c1 .- c2) / (2dV)
        end
    end
end
map(sum, conductances)
heatmap(map(sum, conductances); c=:redsblues, clims=maximum(map(abs ∘ sum, conductances)) .* (-1, 1))
let cs = map(c -> c[1:2, :], conductances)
    heatmap(map(sum, cs); c=:redsblues, clims=maximum(map(abs ∘ sum, cs)) .* (-1, 1))
end
let cs = map(c -> c[3:4, :], conductances)
    heatmap(map(sum, cs); c=:redsblues, clims=maximum(map(abs ∘ sum, cs)) .* (-1, 1))
end
# heatmap(map(c->sum(c[3:4,:]),conductances); c = :redsblues,clims = maximum(map(abs ∘ sum,conductances)) .* (-1,1))
