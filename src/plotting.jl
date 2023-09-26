using CairoMakie
using LaTeXStrings

paramstyle = Dict(:U => L"U_l/\Delta", :V => L"U_{nl}/\Delta", :h => L"V_Z/\Delta",
 :tratio => L"t_{so}/\Delta", :μ1 => L"\epsilon_1/\Delta", :μ2 => L"\epsilon_2/\Delta")
ss_colors = cgrad(:rainbow, categorical=true)[[1, 2, 5, 4]]


 select_best_sweet_spot(spots...) = spots[findmin(MPU, spots)[2]]
 get_property(results, label) = map(x -> x[label], results)
 function combine_results(results)
     testlabels = (:fixedparams, :xlabel, :ylabel, :target, :x, :y)
     props = map(label -> unique(getindex.(eachrow(results), label)), testlabels)
     @assert all(length.(props) .== 1)
     uniqueprops = (; zip(testlabels, map(only, props))...)
     sweet_spots = map(select_best_sweet_spot, results[:, :sweet_spots]...)
     merge(uniqueprops, (; sweet_spots))
 end

struct IntegerTicks end
struct Log10IntegerTicks end
struct LLog10IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin):floor(Int, vmax)
Makie.get_tickvalues(::Log10IntegerTicks, vmin, vmax) = exp10.(ceil(Int, log10(vmin)):floor(Int, log10(vmax)))
function expticks(vmin, vmax)
    vals = exp10.(ceil(Int, log10(vmin)):floor(Int, log10(vmax)))
    (vals, map(exp2text ∘ string ∘ Int ∘ log10, vals))
end
function exp2text(s::AbstractString, base="10")
    codes = Dict(collect("-1234567890") .=> collect("⁻¹²³⁴⁵⁶⁷⁸⁹⁰"))
    return base * map(c -> codes[c], s)
end
function add_exp_colorbar!(pos, hm; kwargs...)
    colorrange = hm.colorrange.val
    ticks = expticks(colorrange...)
    return Colorbar(pos, hm; ticks, kwargs...)
end

## Charge stability
function plot_charge_stability(data; backgroundcolor=:transparent, kwargs...)
    fig = Figure(; resolution=400 .* (1.2, 1), fontsize=20, backgroundcolor)
    ax, hm = plot_charge_stability!(fig, data; kwargs...)
    return fig, ax, hm
end
function plot_charge_stability!(f, data; kwargs...)
    xlabel = paramstyle[:μ1]
    ylabel = paramstyle[:μ2]
    ax = Axis(f; xlabel, ylabel)
    hm = plot_charge_stability!(ax, data; kwargs...)
    return ax, hm
end
function plot_charge_stability!(ax::Axis, data; datamap=x -> x.gap, colormap=:berlin, kwargs...)
    x = data[:ϵ1]
    y = data[:ϵ2]
    gaps = map(datamap, data[:data])
    hm = heatmap!(ax, x, y, gaps; colormap, kwargs...)
    return hm
end

## Plot sweet spot scan
function plot_sweet_scan(data; backgroundcolor=:transparent, kwargs...)
    fig = Figure(; resolution=400 .* (1.2, 1), fontsize=20, backgroundcolor)
    ax, hm = plot_sweet_scan!(fig, data; kwargs...)
    return fig, ax, hm
end
function plot_sweet_scan!(fig::Figure, data, pos=(1, 1); colorbar=true, kwargs...)
    xlabel = paramstyle[data[:xlabel]]
    ylabel = paramstyle[data[:ylabel]]
    ax = Axis(fig[pos...]; xlabel, ylabel)
    hm = plot_sweet_scan!(ax, data; kwargs...)
    if colorbar
        add_exp_colorbar!(fig[pos[1], pos[2]+1], hm)
    end
    return ax, hm
end
function plot_sweet_scan!(ax::Axis, data; datamap=MPU, colormap=Reverse(:viridis), colorrange=(10.0^(-1.5), 1), colorscale=log10, plotbound=(data[:xlabel] == :U && data[:ylabel] == :h))
    x = data[:x]
    y = data[:y]
    z = map(datamap, data[:sweet_spots])
    hm = heatmap!(ax, x, y, z; colorrange, colormap, colorscale)
    if plotbound
        tratio = data[:fixedparams].tratio
        Δ = data[:fixedparams].Δ
        _upperbound(U) = upperbound(U, tratio, Δ)
        _lowerbound(U) = lowerbound(U, tratio, Δ)
        lines!(ax, x, _upperbound, color=:red, linestyle=:dash, linewidth=2)
        lines!(ax, x, _lowerbound, color=:red, linestyle=:dash, linewidth=2)
        ylims!(ax, first(y), last(y))
    end
    return hm
end


lowerbound(U, tratio, Δ) = first(@. Δ - U / 2)
upperbound(U, tratio, Δ) = first(@. Δ * sqrt(1 + tratio^(-2)) - U / 2)