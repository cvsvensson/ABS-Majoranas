abs_hamiltonian(c; μ, Δ, t, tratio, h, ϕ, U, V) = blockdiagonal((QuantumDots.BD1_hamiltonian(c; h, t, μ, Δ=Δ * [exp(1im * ϕ / 2), exp(-1im * ϕ / 2)], Δ1=0, θ=parameter(2atan(tratio), :diff), ϕ=0, U, V)), c)

cost_function(energies, reduced::Number; exp=12.0, minexcgap=0) = cost_reduced(reduced) + cost_energy(energies; exp, minexcgap)
cost_energy(energies; minexcgap=0, exp) = cost_gapratio(gapratio(energies...); exp) + ((excgap(energies...) - minexcgap) < 0 ? 1 + abs(excgap(energies...) - minexcgap) : 0)
cost_gapratio(gr; exp) = abs(gr) > 2 * 10.0^(-exp) ? 1.0 + 10^(exp) * abs2(gr) : abs2(gr)
cost_reduced(reduced) = reduced^2# + 10^3*sqrt(maxreduced)*reduced*(1 + sign(reduced-maxreduced))
refine_interval((a, b), newmid, α::Number) = refine_interval((a, b), newmid, (α, :hard_limit))
function refine_interval((a, b), newmid, (α, s)::Tuple{Number,Symbol})
    if s == :hard_limit
        r = abs(α * (a - b) / 2)
        return (max(newmid - r, a), min(newmid + r, b))
    elseif s == :soft_limit
        r = abs(α * (a - b) / 2)
        return (newmid - r, newmid + r)
    else
        error("Option $s not supported. Supported options are :hard_limit and :soft_limit.")
    end
end

cell_labels(n, basis) = Tuple(keys(QuantumDots.cell(n, basis)))
cell_labels(basis) = Base.Fix2(cell_labels, basis)
function reduced_similarity(basis, oddvec::AbstractVector{T}, evenvec) where {T}
    o = oddvec * oddvec'
    e = evenvec * evenvec'
    fermions = map(label -> norm(QuantumDots.reduced_density_matrix(o, (label,), basis) - QuantumDots.reduced_density_matrix(e, (label,), basis)), keys(basis.dict))
    labels = cell_labels(basis)
    cells = map(n -> norm(QuantumDots.reduced_density_matrix(o, labels(n), basis) -
                          QuantumDots.reduced_density_matrix(e, labels(n), basis)),
        collect(unique(first.(keys(basis)))))
    return (; fermions, cells)
end

function half_majorana_polarizations(majcoeffs, basis)
    keys1 = sort(unique(first.(collect(keys(basis)))))
    N = length(keys1)
    n = div(N + 1, 2)
    keys1L = eachindex(keys1)[1:n]
    keys1L = eachindex(keys1)[1:n]
    keys1R = eachindex(keys1)[end:-1:end-n+1]
    keysL = filter(k -> first(k) in keys1L, keys(basis))
    keysR = filter(k -> first(k) in keys1R, keys(basis))
    left = QuantumDots.majorana_polarization(majcoeffs..., keysL)
    right = QuantumDots.majorana_polarization(majcoeffs..., keysR)
    return (; left, right)
end
function solve(H; basis=c, reduced=true, transport=Transport(missing, missing))
    eig = QuantumDots.diagonalize(H)
    sectors = blocks(eig)
    fullsectors = blocks(eig; full=true)
    oddvals = sectors[1].values
    evenvals = sectors[2].values
    oddvecs = fullsectors[1].vectors
    evenvecs = fullsectors[2].vectors
    majcoeffs = QuantumDots.majorana_coefficients(oddvecs[:, 1], evenvecs[:, 1], basis)
    mps = half_majorana_polarizations(majcoeffs, basis)
    reduced = reduced ? reduced_similarity(basis, oddvecs[:, 1], evenvecs[:, 1]) : missing
    conductance = conductance_matrix(transport, eig)
    return (; gap=first(oddvals) - first(evenvals), gapratio=gapratio(oddvals, evenvals), reduced, mps, majcoeffs, energies=(oddvals, evenvals), conductance)
end


function gapratio(oddvals, evenvals)
    δE = first(oddvals) - first(evenvals)
    Δ = min(oddvals[2], evenvals[2]) - min(first(oddvals), first(evenvals))
    return δE / Δ
end
function expand_searchrange(range, init::Number)
    if first(range) < init < last(range)
        return range
    else
        return range .+ init .- sum(range) / 2
    end
end
excgap(odd, even) = min(odd[2] - odd[1], even[2] - even[1])

abstract type MajoranaQuality end
struct LD <: MajoranaQuality end
struct MP <: MajoranaQuality end
struct MPU <: MajoranaQuality end

Base.@kwdef struct Optimizer2{f,r,i,rf}
    hamfunc::f
    ranges::Vector{r}
    initials::Vector{i}
    MaxTime::Int = 10
    minexcgap::Float64 = 0.0
    exps::Vector{Float64} = Float64.(collect(range(1, 9; length=4)))
    refinefactor::rf = 0.5
    target::MajoranaQuality = LD()
    tracemode::Symbol = :silent
    extra_cost = (x...) -> 0
end
Optimizer = Optimizer2

(::LD)(sol) = norm(sol.reduced.cells)^2
(::MP)(sol) = norm((1 - abs(sol.mps.left.mp)), (1 - abs(sol.mps.right.mp)))
(::MPU)(sol) = norm((1 - abs(sol.mps.left.mpu)), (1 - abs(sol.mps.right.mpu)))

tracemode(opt::Optimizer) = opt.tracemode
function cost(exp, opt::Optimizer)
    function _cost(args)
        sol = solve(opt.hamfunc(args...))
        cost_function(sol.energies, opt.target(sol); exp, opt.minexcgap) + opt.extra_cost(args, exp)
    end
end

function get_sweet_spot(opt::Optimizer)
    refinements = length(opt.exps)
    SearchRange = map(expand_searchrange, opt.ranges, opt.initials)
    NumDimensions = length(SearchRange)
    MaxTime = opt.MaxTime / refinements
    println("Initial point: ", opt.initials)
    println("SearchRange: ", SearchRange)
    res = bboptimize(cost(first(opt.exps), opt), opt.initials; SearchRange, NumDimensions, MaxTime, TraceInterval=10.0, TraceMode=tracemode(opt))
    for exp in Iterators.drop(opt.exps, 1)
        ss::typeof(opt.initials) = best_candidate(res)
        SearchRange = [refine_interval(sr, ss, opt.refinefactor) for (sr, ss) in zip(SearchRange, ss)]
        println("Sweet spot:", ss)
        println("SearchRange:", SearchRange)
        res = bboptimize(cost(exp, opt), ss; SearchRange, NumDimensions, MaxTime, TraceInterval=10.0, TraceMode=tracemode(opt))
    end
    bc = best_candidate(res)
    return bc
end
##
function sweet_spot_scan((xs, xlabel), (ys, ylabel); fixedparams, MaxTime, target)
    iter = collect(Base.product(xs, ys))
    ss = Folds.map((xy) -> get_sweet_spot(; fixedparams..., Dict(xlabel => xy[1], ylabel => xy[2])..., MaxTime, target), iter)
    return Dict(:sweet_spots => ss, xlabel => xs, ylabel => ys, :fixedparams => fixedparams, :MaxTime => MaxTime, :target => target)
end
function charge_stability_scan(parameters, dx=1, dy=1, res = 100)
    μ1s = range(parameters.μ[1] .- dx/2, parameters.μ[1] .+ dx/2; length=res)
    μ2s = range(parameters.μ[1] .- dy/2, parameters.μ[2] .+ dy/2; length=res)
    iter = collect(Base.product(μ1s, μ2s))
    data = Folds.map((xy) -> solve(abs_hamiltonian(c; parameters..., :μ => xy)), iter)
    return Dict(:data => data, :μs => iter, :μ1 => μ1s, :μ2 => μ2s, :fixedparams => fixedparams)
end

## Transport
struct Transport{T,NT}
    type::T
    parameters::NT
end
QuantumDots.conductance_matrix(::Transport{Missing}, eig) = missing
function QuantumDots.conductance_matrix(t::Transport, eig)
    leads = get_leads(c, t.parameters...)
    system = t.type(QuantumDots.diagonalize(QuantumDots.OpenSystem(eig, leads)))
    QuantumDots.conductance_matrix(system)
end
function get_leads(c, T, μ, Γ=1)
    N = div(QuantumDots.nbr_of_fermions(c), 2)
    left = QuantumDots.CombinedLead((Γ * c[1, :↑]', Γ * c[1, :↓]'); T, μ=μ[1])
    right = QuantumDots.CombinedLead((Γ * c[N, :↑]', Γ * c[N, :↓]'); T, μ=μ[2])
    return (; left, right)
end

## Kitaev parameters
function _Δk(t, tratio, s, c, sqp, sqm, pf)
    t * tratio * pf * (-1im * s * (sqp[1] * sqm[2] - sqm[1] * sqp[2]) +
                       c * (sqp[1] * sqm[2] + sqm[1] * sqp[2]))
end
function _tk(t, tratio, s, c, sqp, sqm, pf)
    t * pf * (c * (prod(sqm) - prod(sqp)) - 1im * s * (prod(sqm) + prod(sqp)))
end
function _μk(h, U, μ, β, Δ, sqm, sqp, V)
    (-h - U * (μ[1] + β[1]) / (2β[1]) + (μ[1]^2 + Δ * sqm[1] * sqp[1]) / β[1] - V * μ[1] * (β[2] + μ[2]) / prod(β),
        -h - U * (μ[2] + β[2]) / (2β[2]) + (μ[2]^2 + Δ * sqm[2] * sqp[2]) / β[2] - V * μ[2] * (β[1] + μ[1]) / prod(β))
end
function KitaevParameters(; t, tratio, Δ, ϕ, μ, U, V, h)
    β, sqp, sqm = get_β_sq(μ, Δ)
    pf = 1 / (2 * sqrt(prod(β)))
    s, c = sincos(ϕ / 2)
    Δk = _Δk(t, tratio, s, c, sqp, sqm, pf)
    tk = _tk(t, tratio, s, c, sqp, sqm, pf)
    Vk = V * prod(μ) / prod(β)
    μk = _μk(h, U, μ, β, Δ, sqm, sqp, V)
    (; μ=μk, t=tk, V=Vk, Δ=Δk)
end
function kitaev_μ(; h, U, μ, Δ, V)
    β, sqp, sqm = get_β_sq(μ, Δ)
    _μk(h, U, μ, β, Δ, sqm, sqp, V)
end
function get_β_sq(μ, Δ)
    β = (sqrt(μ[1]^2 + Δ^2), sqrt(μ[2]^2 + Δ^2))# map(μ -> sqrt(μ^2 + Δ^2), μ)
    sqp = (sqrt(β[1] + μ[1]), sqrt(β[2] + μ[2]))# map((β, μ) -> sqrt(β + μ), β, μ)
    sqm = (sqrt(β[1] - μ[1]), sqrt(β[2] - μ[2]))# map((β, μ) -> sqrt(β - μ), β, μ)
    return β, sqp, sqm
end
function kitaev_tΔ(; t, tratio, μ, Δ)
    β, sqp, sqm = get_β_sq(μ, Δ)
    pf = 1 / (2 * sqrt(prod(β)))
    δϕ -> begin
        s, c = sincos(δϕ / 2)
        Δk = _Δk(t, tratio, s, c, sqp, sqm, pf)
        tk = _tk(t, tratio, s, c, sqp, sqm, pf)
        (tk, Δk)
    end
end

function kitaev_sweet_spot2(; t, tratio, Δ, U, V, h)
    guess = [sqrt((h + U / 2)^2 - Δ^2) + U / 2, -sqrt((h + U / 2)^2 - Δ^2) + U / 2, pi / 2.0]
    function cost((μ1, μ2, ϕ))
        pK = KitaevParameters(; t, tratio, Δ, U, V, h, μ=(μ1, μ2), ϕ)
        abs2(abs(pK.t) + pK.V / 2 - abs(pK.Δ)) + sum(abs2, pK.μ .+ pK.V / 2)
    end
    sol = bboptimize(cost, guess;
        SearchRange=[map(μ -> (μ - 1, μ + 1), guess[1:2])..., (0, pi)], TargetFitness=1e-6, TraceMode=:silent)
    bc = best_candidate(sol)
    (; μ=(bc[1], bc[2]), ϕ=bc[3])
end

kitaev_μ_zero(Δ, h, U) = U / 2 .+ (1, -1) .* (sqrt(-4Δ^2 + U^2 + 4h * U + 4h^2) / 2)

function sweet_spot_ϕ(; t, tratio, μ, Δ, V)
    μ1, μ2 = μ
    acos((t^2 * (-1 + tratio^2) * V^2 * sqrt(Δ^4) * μ1^2 * μ2^2 -
          2 * sqrt(Δ^2 + μ1^2) * sqrt(Δ^2 + μ2^2) * sqrt(
              t^4 * tratio^2 * V^2 * Δ^4 * μ1^2 * μ2^2 * (4 * t^2 * (1 + tratio^2) -
                                                          (V^2 * μ1^2 * μ2^2) / ((Δ^2 + μ1^2) * (Δ^2 + μ2^2)))) -
          2 * t^4 * (1 + tratio^2) * sqrt(Δ^4) * ((-1 + tratio^2) * Δ^4 +
                                                  (-1 + tratio^2) * Δ^2 * (μ1^2 + μ2^2) + μ1 * μ2 *
                                                                                          ((-1 + tratio^2) * μ1 * μ2 - (1 + tratio^2) *
                                                                                                                       sqrt(Δ^2 + μ1^2) * sqrt(Δ^2 + μ2^2)))) /
         (2 * t^4 * (1 + tratio^2)^2 * Δ^4 * sqrt(Δ^2 + μ1^2) * sqrt(Δ^2 + μ2^2)))
end
