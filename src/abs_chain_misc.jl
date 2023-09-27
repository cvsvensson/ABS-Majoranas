using QuantumDots, QuantumDots.BlockDiagonals, LinearAlgebra
using BlackBoxOptim # For optimizing 
using Folds # For multithreading
using ForwardDiff, LinearSolve # For transport

const N = 2
const c = FermionBasis((1, 2), (:↑, :↓); qn=QuantumDots.parity)
abs_hamiltonian(c; μ1, μ2, Δ, t, tratio, h, ϕ, U, V) = blockdiagonal((QuantumDots.BD1_hamiltonian(c; h, t, μ=(μ1, μ2), Δ=Δ .* [exp(1im * ϕ / 2), exp(-1im * ϕ / 2)], Δ1=0, θ=parameter(2atan(tratio), :diff), ϕ=0, U, V)), c)
cost_function(energies, reduced::Number; exp=12.0, minexcgap=0) = cost_reduced(reduced) + cost_energy(energies; exp, minexcgap)
cost_energy(energies; minexcgap=0, exp) = cost_gapratio(gapratio(energies...); exp) + ((excgap(energies...) - minexcgap) < 0 ? 1 + abs(excgap(energies...) - minexcgap) : 0)
cost_gapratio(gr; exp) = abs(gr) > 2 * 10.0^(-exp) ? 1.0 + 10^(exp) * abs2(gr) : abs2(gr)
cost_reduced(reduced) = reduced^2


cell_labels(n, basis) = Tuple(keys(QuantumDots.cell(n, basis)))
cell_labels(basis) = Base.Fix2(cell_labels, basis)
function reduced_similarity(basis, oddvec::AbstractVector{T}, evenvec) where {T}
    o = oddvec * oddvec'
    e = evenvec * evenvec'
    fermions = map(label -> norm(QuantumDots.reduced_density_matrix(o - e, (label,), basis)), keys(basis.dict))
    labels = cell_labels(basis)
    cells = map(n -> norm(QuantumDots.reduced_density_matrix(o - e, labels(n), basis)),
        spatial_labels(basis))
    return (; fermions, cells)
end

function half_majorana_polarizations(majcoeffs, basis)
    keys1 = spatial_labels(basis)
    N = length(keys1)
    n = div(N, 2)
    keys1L = keys1[1:n]
    keys1R = keys1[end:-1:end-n+1]
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
    vacuumnorms = (;odd = map(f -> norm(f*oddvecs[:, 1]), basis.dict), even = map(f -> norm(f*evenvecs[:, 1]), basis.dict))
    return (; gap=first(oddvals) - first(evenvals), gapratio=gapratio(oddvals, evenvals), reduced, mps, majcoeffs, energies=(oddvals, evenvals), conductance, vacuumnorms)
end


function gapratio(oddvals, evenvals)
    δE = first(oddvals) - first(evenvals)
    Δ = min(oddvals[2], evenvals[2]) - min(first(oddvals), first(evenvals))
    return δE / Δ
end
excgap(odd, even) = min(odd[2] - odd[1], even[2] - even[1])
Base.@kwdef struct Optimizer{f,r,i,t,ec}
    hamfunc::f
    ranges::Vector{r}
    initials::Vector{i}
    MaxTime::Int = 10
    minexcgap::Float64 = 0.0
    exps::Vector{Float64} = Float64.(collect(range(0.5, 3; length=4)))
    target::t = LD
    tracemode::Symbol = :silent
    extra_cost::ec = (x...) -> 0
    Method::Symbol = :probabilistic_descent
    PopulationSize::Int = 100
    TargetFitness::Float64 = 0.0
end

LD(sol) = norm(sol.reduced.cells)^2
MP(sol) = 1 - (abs(sol.mps.left.mp) + abs(sol.mps.right.mp)) / 2 #norm((1 - abs(sol.mps.left.mp)), (1 - abs(sol.mps.right.mp)))
MPU(sol) = 1 - (abs(sol.mps.left.mpu) + abs(sol.mps.right.mpu)) / 2 #norm((1 - abs(sol.mps.left.mpu)), (1 - abs(sol.mps.right.mpu)))

tracemode(opt::Optimizer) = opt.tracemode
function cost(exp, opt::Optimizer)
    function _cost(args)
        sol = solve(opt.hamfunc(args...))
        cost_function(sol.energies, opt.target(sol); exp, opt.minexcgap) + opt.extra_cost(args, exp)
    end
end

function get_sweet_spot(opt::Optimizer)
    refinements = length(opt.exps)
    SearchRange = opt.ranges #map(expand_searchrange, opt.ranges, opt.initials)
    NumDimensions = length(SearchRange)
    MaxTime = opt.MaxTime / refinements
    TargetFitness = opt.TargetFitness
    println("Initial point: ", opt.initials)
    println("SearchRange: ", SearchRange)
    Method = opt.Method
    res = bboptimize(cost(first(opt.exps), opt), opt.initials;
        SearchRange, NumDimensions, MaxTime, TraceInterval=10.0,
        TraceMode=tracemode(opt),
        Method, PopulationSize=opt.PopulationSize,
        TargetFitness)
    for exp in Iterators.drop(opt.exps, 1)
        ss::typeof(opt.initials) = best_candidate(res)
        #SearchRange = [refine_interval(sr, ss, opt.refinefactor) for (sr, ss) in zip(SearchRange, ss)]
        println("Sweet spot:", ss)
        println("SearchRange:", SearchRange)
        res = bboptimize(cost(exp, opt), ss; SearchRange, NumDimensions,
            MaxTime, TraceInterval=10.0, TraceMode=tracemode(opt),
            Method, TargetFitness,
            PopulationSize=opt.PopulationSize)
    end
    return res
end

##
function anti_parallel_sweet_spot(; Δ, tratio, h, U, V, t, MaxTime, exps=collect(range(0.5, 3, length=4)), target, kwargs...)
    fixedparams = (; Δ, tratio, h, U, V, t)
    μ1, μ2 = kitaev_μ_zero(Δ, h, U)
    ϕ = 0.5 * pi
    hamfunc(ϕ, μ1, μ2) = abs_hamiltonian(c; μ1, μ2, ϕ, fixedparams...)
    maxh = max(abs.(h)...)
    maxU = max(abs.(U)...)
    opt = Optimizer(;
        hamfunc,
        ranges=[(0.0, 1.0π), (0.0, 1.1 * μ1 + maxh + maxU), (-maxh - maxU + μ2, maxU + V)],
        initials=Float64.([ϕ, μ1, μ2]),
        MaxTime, exps, target,
        tracemode=:silent,
        extra_cost=((ϕ, μ1, μ2), e) -> exp(-(e * abs(μ1 - μ2) + 1)^4) + e * (μ2 > μ1) * (μ2 - μ1),
        kwargs...)
    ss = best_candidate(get_sweet_spot(opt))
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

function sweet_spot_scan((xs, xlabel), (ys, ylabel), get_ss=anti_parallel_sweet_spot; fixedparams, MaxTime, target, Method = :adaptive_de_rand_1_bin_radiuslimited, kwargs...)
    iter = collect(Base.product(xs, ys))
    ss = Folds.map((xy) -> get_ss(; fixedparams..., Dict(xlabel => xy[1], ylabel => xy[2])..., MaxTime, target, Method, kwargs...), iter)
    return Dict(:sweet_spots => ss, :x => xs, :y => ys, :xlabel => xlabel, :ylabel => ylabel,
        :fixedparams => fixedparams, :MaxTime => MaxTime, :target => target, :Method=>Method)
end
function charge_stability_scan(parameters, dx=1, dy=1, res=100; transport=missing)
    ϵ1s = -range(parameters.μ1 .- dx / 2, parameters.μ1 .+ dx / 2; length=res)
    ϵ2s = -range(parameters.μ2 .- dy / 2, parameters.μ2 .+ dy / 2; length=res)
    iter = collect(Base.product(ϵ1s, ϵ2s))
    data = Folds.map((xy) -> solve(abs_hamiltonian(c; parameters..., :μ1 => -xy[1], :μ2 => -xy[2]); transport), iter)
    return Dict(:data => data, :ϵs => iter, :ϵ1 => ϵ1s, :ϵ2 => ϵ2s, :parameters => parameters)
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
    l = spatial_labels(c)
    left = QuantumDots.CombinedLead((Γ * c[l[1], :↑]', Γ * c[l[1], :↓]'); T, μ=μ[1])
    right = QuantumDots.CombinedLead((Γ * c[l[end], :↑]', Γ * c[l[end], :↓]'); T, μ=μ[2])
    return (; left, right)
end

## Kitaev parameters
kitaev_μ_zero(Δ, h, U) = @. U / 2 + (1, -1) * (sqrt(abs(-4Δ^2 + U^2 + 4h * U + 4h^2)) / 2)
