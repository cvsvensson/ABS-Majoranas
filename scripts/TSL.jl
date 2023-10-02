
using DrWatson
@quickactivate "Majorana sweet spot"
includet(srcdir("abs_chain_misc.jl"))
includet(srcdir("plotting.jl"))

d = FermionBasis((:L, :C, :R), (:↑, :↓); qn=QuantumDots.parity)
##
Δ = 1.0
t = 0.5Δ
tsoc = 0.2t
U = 5Δ
h = 1.5Δ / 2
##
transport = Transport(QuantumDots.Pauli(), (; T=1 / 40, μ=(0.0, 0.0), Γ=1))
params = (; Δ, U, t, tsoc, h)
hamLCR(μL, μC, μR) = blockdiagonal(QuantumDots.TSL_hamiltonian(d; params..., μC, μL, μR), d)

## Fig3
res = 50
μ1s = range(-1, 1, res)
μ2s = range(-1, 1, res)
sols3 = [[solve(hamLCR(μ1, μC, μ2); basis=d) for μ1 in μ1s, μ2 in μ2s] for μC in (-0.313, -0.7, 0.0) .* params.Δ]
##
f = Figure(; resolution=400 .* (3, 1), fontsize=25);
g = f[1, 1] = GridLayout();
hms = [heatmap(g[1, n], μ1s, μ2s, map(ss -> ss.gap, data), colormap=:redsblues, colorrange=0.8 .* (-1, 1)) for (n, data) in enumerate(sols3)]
Colorbar(g[1, end+1], hms[1].plot)
f
##
res = 50
μ1s = range(-1, 1, res)
μ2s = range(-1, 1, res)
hamLR(μL, μR) = blockdiagonal(QuantumDots.TSL_hamiltonian(d; params..., μC=-0.558 * params.Δ, μL, μR), d)
sols = [solve(hamLR(μ1, μ2); basis=d, transport) for μ1 in μ1s, μ2 in μ2s]

##
heatmap(μ1s, μ2s, map(ss -> ss.conductance[1, 2], sols); colormap=:redsblues, colorrange=pi * transport.parameters.Γ^2 .* (-1, 1))
heatmap(μ1s, μ2s, map(ss -> ss.mps.left.mpu, sols); colormap=:viridis, colorrange=(0, 1))
heatmap(μ1s, μ2s, map(ss -> ss.gap, sols); colormap=:redsblues, colorrange=0.8 .* (-1, 1))


## Let's make Figure 2
res = 50
μs = range(-0.4, 0.4, res)
μCs = range(-1, 1, res)
hamGC(μG, μC) = blockdiagonal(QuantumDots.TSL_hamiltonian(d; params..., μC, μL=μG, μR=μG), d)
fig2data = [solve(hamGC(μ, μC); basis=d) for μC in μCs, μ in μs]
##
fig2 = Figure();
fig2g = fig2[1, 1] = GridLayout()
fig2ax1 = Axis(fig2g[1, 1]; xlabel=L"\mu_C", ylabel=L"\mu")
fig2ax2 = Axis(fig2g[2, 1]; xlabel=L"\mu_C", ylabel=L"\mu")
heatmap!(fig2ax1, μCs, μs, map(ss -> (ss.gap), fig2data); colormap=:redsblues, colorrange=0.1 .* (-1, 1))
heatmap!(fig2ax2, μCs, μs, map(ss -> ss.mps.left.mpu, fig2data); colormap=:viridis, colorrange=(0, 1))
fig2
## Add in optimized sweet spot
includet(srcdir("abs_chain_misc.jl"))
params_abs, params_tsl = let Δ = 1.0, tratio = 0.2, t, h
    h = 1.25Δ
    t = 0.5Δ
    t2 = t*1.2
    U = 1Δ
    (; Δ, U, t, tratio=0.2, h, V=0),
    (; Δ, U, t = t2, tsoc=t2 * tratio, h)
end
tsl_ham(μG, μC) = blockdiagonal(QuantumDots.TSL_hamiltonian(d; params_tsl..., μC, μL=μG, μR=μG), d)
abs_ham(μ1, μ2, δϕ) = abs_hamiltonian(c; μ1, μ2, ϕ=δϕ, params_abs...)
function cost(x; ham, basis, P)
    sol = solve(ham(x...); basis)
    P * abs2(sol.gap) + MP(sol)
end
tsl_cost(x; P) = cost(x; ham=tsl_ham, basis=d, P)
abs_cost(x; P) = cost(x; ham=abs_ham, basis=c, P)

##
optparams = (;  MaxTime=5, TargetFitness=1e-6)
tsl_opt = bboptimize(x -> tsl_cost(x; P=1e3); SearchRange = [(-5,5),(-5,5)], NumDimensions=2, optparams...)
abs_opt = bboptimize(x -> abs_cost(x; P=1e3), [4,-4, pi/2];  SearchRange = [(0,20),(-20,0),(0,pi)], NumDimensions=3, optparams...)
tsl_ss = best_candidate(tsl_opt)
abs_ss = best_candidate(abs_opt)
tsl_ss = solve(tsl_ham(tsl_ss...); basis=d)
abs_ss = solve(abs_ham(abs_ss...); basis=c)
excgap(tsl_ss.energies...)
excgap(abs_ss.energies...)
Dict(pairs(tsl_ss))
Dict(pairs(abs_ss))


##
scatter!(fig2ax1, [ss[2]], [ss[1]])
scatter!(fig2ax2, [ss[2]], [ss[1]])
fig2