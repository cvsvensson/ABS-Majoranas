allparams = (:U, :V, :Δ, :tratio, :t, :ϕ, :μ, :h)

struct HamGen{F,D,V,f,fA}
    fixedparams::F
    dynamicparams::D
    vecparams::V
    func::f
    args::fA
end
(H::HamGen)(args...) = H.func(args...)
function ham(fixedparams, vecparams=())
    dynamicparams = Tuple(setdiff(allparams, keys(fixedparams)))
    _process_args(args...) = process_args(dynamicparams, vecparams, args...)
    f(args...) = abs_hamiltonian(c; fixedparams..., _process_args(args...)...)
    args = Symbol[]
    for arg in dynamicparams
        if arg in vecparams
            foreach(j->push!(args, Symbol(arg, j)), 1:N)
        else
            push!(args, arg)
        end
    end
    HamGen(fixedparams, dynamicparams, vecparams, f, Tuple(args))
end
process_args(dynamicparams, vecparams::Nothing, args...) = NamedTuple(dynamicparams .=> args)
function process_args(dynamicparams, vecparams, args...)
    count = 1
    countdp = 1
    inds = eachindex(args)
    pairs = []
    while count <= length(inds)
        sym = dynamicparams[countdp]
        countdp += 1
        if sym in vecparams
            push!(pairs, sym => args[inds[count:count+1]])
            count += 2
        else
            push!(pairs, sym => args[inds[count]])
            count += 1
        end
    end
    NamedTuple(pairs)
end