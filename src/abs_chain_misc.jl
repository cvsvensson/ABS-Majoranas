using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using BlackBoxOptim # For optimizing 
using Folds # For multithreading
using ForwardDiff, LinearSolve # For transport

const N = 2
const c = FermionBasis((1, 2), (:↑, :↓); qn=QuantumDots.parity)
abs_hamiltonian(c; μ1, μ2, Δ, t, tratio, h, ϕ, U, V) = blockdiagonal((QuantumDots.BD1_hamiltonian(c; h, t, μ=(μ1, μ2), Δ=Δ * [exp(1im * ϕ / 2), exp(-1im * ϕ / 2)], Δ1=0, θ=parameter(2atan(tratio), :diff), ϕ=0, U, V)), c)
abs_hamiltonian_odd(c; μ1, μ2, Δ, t, tratio, h, ϕ, U, V) = (n = div(QuantumDots.nbr_of_fermions(c), 2);
QuantumDots.BD1_hamiltonian(c; h, t, μ=(μ1, μ2), Δ=Δ * [exp(1im * ϕ / 2), exp(-1im * ϕ / 2)], Δ1=0, θ=parameter(2atan(tratio), :diff), ϕ=0, U, V)[1:2^n, 1:2^n])
cost_function(energies, reduced::Number; exp=12.0, minexcgap=0) = cost_reduced(reduced) + cost_energy(energies; exp, minexcgap)
cost_energy(energies; minexcgap=0, exp) = cost_gapratio(gapratio(energies...); exp) + ((excgap(energies...) - minexcgap) < 0 ? 1 + abs(excgap(energies...) - minexcgap) : 0)
cost_gapratio(gr; exp) = abs(gr) > 2 * 10.0^(-exp) ? 1.0 + 10^(exp) * abs2(gr) : abs2(gr)
cost_reduced(reduced) = reduced^2

cost_function_borg(energies, reduced::Number; exp=12.0, minexcgap=0) = (cost_reduced(reduced), 10.0^(2exp) * abs2(gapratio(energies...)))


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
function expand_searchrange(range, init::Number)
    if first(range) < init < last(range)
        return range
    else
        return range .+ init .- sum(range) / 2
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
        spatial_labels(basis))
    return (; fermions, cells)
end

function half_majorana_polarizations(majcoeffs, basis)
    keys1 = spatial_labels(basis)
    N = length(keys1)
    n = div(N + 1, 2)
    keys1L = (keys1)[1:n]
    keys1R = (keys1)[end:-1:end-n+1]
    keysL = filter(k -> first(k) in keys1L, keys(basis))
    keysR = filter(k -> first(k) in keys1R, keys(basis))
    left = QuantumDots.majorana_polarization(majcoeffs..., keysL)
    right = QuantumDots.majorana_polarization(majcoeffs..., keysR)
    return (; left, right)
end
function solve(H; basis=c, reduced=true, transport=missing)
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
    conductance = conductance_matrix(transport, eig; basis)
    return (; gap=first(oddvals) - first(evenvals), gapratio=gapratio(oddvals, evenvals), reduced, mps, majcoeffs, energies=(oddvals, evenvals), conductance)
end


function gapratio(oddvals, evenvals)
    δE = first(oddvals) - first(evenvals)
    Δ = min(oddvals[2], evenvals[2]) - min(first(oddvals), first(evenvals))
    return δE / Δ
end
excgap(odd, even) = min(odd[2] - odd[1], even[2] - even[1])
Base.@kwdef struct Optimizer4{f,r,i,rf}
    hamfunc::f
    ranges::Vector{r}
    initials::Vector{i}
    MaxTime::Int = 10
    minexcgap::Float64 = 0.0
    exps::Vector{Float64} = Float64.(collect(range(1, 9; length=4)))
    refinefactor::rf = 0.5
    target = LD
    tracemode::Symbol = :silent
    extra_cost = (x...) -> 0
    PopulationSize = 100
    ϵ = 0.001
end
Optimizer = Optimizer4

LD(sol) = norm(sol.reduced.cells)^2
MP(sol) = norm((1 - abs(sol.mps.left.mp)), (1 - abs(sol.mps.right.mp)))
MPU(sol) = norm((1 - abs(sol.mps.left.mpu)), (1 - abs(sol.mps.right.mpu)))

tracemode(opt::Optimizer) = opt.tracemode
function cost(exp, opt::Optimizer)
    function _cost(args)
        sol = solve(opt.hamfunc(args...))
        cost_function(sol.energies, opt.target(sol); exp, opt.minexcgap) + opt.extra_cost(args, exp)
    end
end
function cost_borg(exp, opt::Optimizer)
    function _cost(args)
        sol = solve(opt.hamfunc(args...))
        cost_function_borg(sol.energies, opt.target(sol); exp, opt.minexcgap)
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
    # bc = best_candidate(res)
    return res
end

function get_sweet_spot_borg(opt::Optimizer)
    refinements = length(opt.exps)
    SearchRange = map(expand_searchrange, opt.ranges, opt.initials)
    NumDimensions = length(SearchRange)
    MaxTime = opt.MaxTime / refinements
    println("Initial point: ", opt.initials)
    println("SearchRange: ", SearchRange)
    res = bboptimize(cost_borg(first(opt.exps), opt), opt.initials;
        Method=:borg_moea, PopulationSize=opt.PopulationSize,
        FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true), ϵ=opt.ϵ,
        SearchRange, NumDimensions, MaxTime, TraceInterval=10.0, TraceMode=tracemode(opt))
    for exp in Iterators.drop(opt.exps, 1)
        ss::typeof(opt.initials) = best_candidate(res)
        SearchRange = [refine_interval(sr, ss, opt.refinefactor) for (sr, ss) in zip(SearchRange, ss)]
        println("Sweet spot:", ss)
        println("SearchRange:", SearchRange)
        res = bboptimize(cost_borg(exp, opt), ss; Method=:borg_moea,
            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true), ϵ=opt.ϵ, PopulationSize=opt.PopulationSize,
            SearchRange, NumDimensions, MaxTime, TraceInterval=10.0, TraceMode=tracemode(opt))
    end
    # bc = best_candidate(res)
    return res
end
##
function anti_parallel_sweet_spot(; Δ, tratio, h, U, V, t, MaxTime, exps=collect(range(0.5, 3, length=4)), target, PopulationSize=100, ϵ=0.001)
    fixedparams = (; Δ, tratio, h, U, V, t)
    pk = kitaev_sweet_spot_guess(; Δ, h, U, V, t, tratio)
    hamfunc(ϕ, μ1, μ2) = abs_hamiltonian(c; μ1, μ2, ϕ, fixedparams...)
    opt = Optimizer(;
        hamfunc,
        ranges=[(0.0, 1.0π), (0.0, 1.1 * pk.μ1 + abs(h) + U), (-abs(2h) - U, U + V)],
        initials=Float64.([pk.ϕ, pk.μ1, pk.μ2]),
        MaxTime, exps, target,
        refinefactor=(1, :hard_limit),
        tracemode=:silent,
        extra_cost=((ϕ, μ1, μ2), e) -> exp(-(e * abs(μ1 - μ2) + 1)^4),
        PopulationSize, ϵ)
    ss = best_candidate(get_sweet_spot_borg(opt))
    optsol = solve(opt.hamfunc(ss...))
    parameters = merge(fixedparams, NamedTuple(zip((:ϕ, :μ1, :μ2), ss)))
    optimization = opt
    sweet_spot = merge(optsol, (; parameters, optimization))
    return sweet_spot
end
function parallel_sweet_spot(; Δ, tratio, h, U, V, t, MaxTime, exps=collect(range(0.5, 3, length=4)), target, lower=true)
    fixedparams = (; Δ, tratio, h, U, V, t)
    pk = kitaev_sweet_spot_guess(; Δ, h, U, V, t, tratio)
    μinitial = lower ? pk.μ2 : pk.μ1
    μrange = lower ? (-2abs(h) - U, U + V) : (0.0, 1.1 * abs(pk.μ1) + abs(h) + U)
    hamfunc(ϕ, μ) = abs_hamiltonian(c; μ1=μ, μ2=μ, ϕ, fixedparams...)
    opt = Optimizer(;
        hamfunc,
        ranges=[(0.0, 1.0π), μrange],
        initials=Float64.([pk.ϕ, μinitial]),
        MaxTime, exps, target,
        refinefactor=(1, :hard_limit),
        tracemode=:silent)
    ss = get_sweet_spot(opt)
    optsol = solve(opt.hamfunc(ss...))
    parameters = merge(fixedparams, NamedTuple(zip((:ϕ, :μ1, :μ2), ss[[1, 2, 2]])))
    optimization = opt
    sweet_spot = merge(optsol, (; parameters, optimization))
    return sweet_spot
end

function sweet_spot_scan((xs, xlabel), (ys, ylabel), get_ss=anti_parallel_sweet_spot; fixedparams, MaxTime, target)
    iter = collect(Base.product(xs, ys))
    ss = Folds.map((xy) -> get_ss(; fixedparams..., Dict(xlabel => xy[1], ylabel => xy[2])..., MaxTime, target), iter)
    return Dict(:sweet_spots => ss, :x => xs, :y => ys, :xlabel => xlabel, :ylabel => ylabel,
        :fixedparams => fixedparams, :MaxTime => MaxTime, :target => target)
end
function charge_stability_scan(parameters, dx=1, dy=1, res=100; transport=missing)
    μ1s = range(parameters.μ1 .- dx / 2, parameters.μ1 .+ dx / 2; length=res)
    μ2s = range(parameters.μ2 .- dy / 2, parameters.μ2 .+ dy / 2; length=res)
    iter = collect(Base.product(μ1s, μ2s))
    data = Folds.map((xy) -> solve(abs_hamiltonian(c; parameters..., :μ1 => xy[1], :μ2 => xy[2]); transport), iter)
    return Dict(:data => data, :μs => iter, :μ1 => μ1s, :μ2 => μ2s, :parameters => parameters)
end

## Transport
struct Transport{T,NT}
    type::T
    parameters::NT
end
QuantumDots.conductance_matrix(::Missing, eig; kwargs...) = missing
function QuantumDots.conductance_matrix(t::Transport, eig; basis=c)
    leads = get_leads(basis, t.parameters...)
    system = t.type(QuantumDots.diagonalize(QuantumDots.OpenSystem(eig, leads)))
    QuantumDots.conductance_matrix(system)
end
spatial_labels(basis) = collect(unique(first.(keys(basis))))
function get_leads(c, T, μ, Γ=1)
    # N = div(QuantumDots.nbr_of_fermions(c), 2)
    l = spatial_labels(c)
    left = QuantumDots.CombinedLead((Γ * c[l[1], :↑]', Γ * c[l[1], :↓]'); T, μ=μ[1])
    right = QuantumDots.CombinedLead((Γ * c[l[end], :↑]', Γ * c[l[end], :↓]'); T, μ=μ[2])
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

function kitaev_sweet_spot_guess(; t, tratio, Δ, U, V, h)
    guess = [sqrt(abs((h + U / 2)^2 - Δ^2)) + U / 2, -sqrt(abs((h + U / 2)^2 - Δ^2)) + U / 2, pi / 2.0]
    function cost((μ1, μ2, ϕ))
        pK = KitaevParameters(; t, tratio, Δ, U, V, h, μ=(μ1, μ2), ϕ)
        abs2(abs(pK.t) + pK.V / 2 - abs(pK.Δ)) + sum(abs2, pK.μ .+ pK.V / 2)
    end
    sol = bboptimize(cost, guess;
        SearchRange=[map(μ -> (μ - 1, μ + 1), guess[1:2])..., (0, pi)], TargetFitness=1e-6, TraceMode=:silent)
    bc = best_candidate(sol)
    (; μ1=bc[1], μ2=bc[2], ϕ=bc[3])
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
