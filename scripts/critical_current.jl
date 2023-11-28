
## Critical current
ϕs = ssparams.ϕ .+ range(-pi, pi, 100)
dsols = [dsolve(p -> abs_hamiltonian(c; p...), merge(ssparams, (; ϕ)); dp=(; ϕ=1e-6)) for ϕ in ϕs]
ssparams.ϕ
sols = [solve(abs_hamiltonian(c; ssparams..., ϕ)) for ϕ in ϕs]
plot(ϕs, map(LD, sols))
plot(ϕs, map(x -> x.gap, sols))
energies = reduce(hcat, map(x -> [x.energies[1][1:2]..., x.energies[2][1:2]...], sols))
labels = ["o1", "o2", "e1", "e2"]
fig, ax, sp = series(ϕs, energies; labels);
axislegend(ax);
fig 

fig2, ax2, sp2 = series(ϕs[1:end-1], diff(energies; dims=2); labels = "d" .* labels, axis=(; limits=(nothing, 0.01 .* (-1, 1))));
axislegend(ax2);
fig2 

fig3, ax3, sp3 = series(ϕs[1:end-2], diff(diff(energies; dims=2); dims=2); labels = "d²" .* labels, axis=(; limits=(nothing, 0.001 .* (-1, 1))));
axislegend(ax3);
fig3

plot(ϕs[1:end-1], map(x -> x.dsol.energies[1][1], dsols) |> diff)
plot(ϕs[1:end-1], map(x -> x.dsol.energies[1][2], dsols) |> diff)
