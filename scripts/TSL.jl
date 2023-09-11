using QuantumDots

d = FermionBasis((:L, :C, :R), (:↑, :↓); qn=QuantumDots.parity)
##
Δ = 1.0
t = 0.5Δ
##
transport = Transport(QuantumDots.Pauli(), (; T=1 / 40, μ=(0.0, 0.0), Γ=1 / 1))
params = (; Δ, U=5Δ, t, tsoc=0.2t, h=1.5Δ / 2)
ham(μL, μC, μR) = blockdiagonal(QuantumDots.TSL_hamiltonian(d; params..., μC, μL, μR), d)

## Fig3
res = 50
μ1s = range(-1, 1, res)
μ2s = range(-1, 1, res)
sols3 = [[solve(ham(μ1, μC, μ2); basis=d) for μ1 in μ1s, μ2 in μ2s] for μC in (-.313, -.7, 0.0 ) .* params.Δ]
##
f = Figure(; resolution = 400 .* (3,1), fontsize=25);
g = f[1,1] = GridLayout();
hms = [heatmap(g[1,n], μ1s,μ2s, map(ss->ss.gap, data), colormap=:redsblues, colorrange=0.8 .* (-1, 1)) for (n, data) in enumerate(sols3)]
Colorbar(g[1,end+1], hms[1].plot)
f
##
res = 50
μ1s = range(-1, 1, res)
μ2s = range(-1, 1, res)
ham(μL, μR) = blockdiagonal(QuantumDots.TSL_hamiltonian(d; params..., μC=-0.558 * params.Δ, μL, μR), d)
sols = [solve(ham(μ1, μ2); basis=d, transport) for μ1 in μ1s, μ2 in μ2s]

##
heatmap(μ1s, μ2s, map(ss -> ss.conductance[1, 2], sols); colormap=:vik, colorrange=pi * transport.parameters.Γ^2 .* (-1, 1))
heatmap(μ1s, μ2s, map(ss -> ss.mps.left.mpu, sols); colormap=:viridis, colorrange=(0, 1))
heatmap(μ1s, μ2s, map(ss -> ss.gap, sols); colormap=:vik, colorrange= 0.8 .* (-1, 1))

##
res = 100
μs = range(-.4, .4, res)
μCs = range(-1, 1, res)
ham(μ, μC) = blockdiagonal(QuantumDots.TSL_hamiltonian(d; params..., μC, μL=μ, μR=μ), d)
sols = [solve(ham(μ, μC); basis=d) for μC in μCs, μ in μs]
heatmap(μCs, μs, map(ss -> ss.mps.left.mpu, sols); colormap=:viridis, colorrange=(0, 1))
heatmap(μCs, μs, map(ss -> (ss.gap), sols); colormap=:vik, colorrange= 0.1 .* (-1, 1))


##