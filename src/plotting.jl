using CairoMakie
using LaTeXStrings

paramstyle = Dict(:U => L"U", :V => L"V", :h => L"V_Z", :tratio => L"t_{so}", :μ1 => L"μ_1", :μ2 => L"μ_2")

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
    cbar = Colorbar(pos, hm; ticks, kwargs...)
    return nothing
end

## Charge stability
function plot_charge_stability(data; kwargs...)
    fig = Figure(resolution=400 .* (1.2, 1), fontsize=20, backgroundcolor=:transparent)
    ax, hm = plot_charge_stability!(fig, data; kwargs...)
    return fig, ax, hm
end
function plot_charge_stability!(fig::Figure, data, pos=(1, 1); colorbar=true, kwargs...)
    xlabel = paramstyle[:μ1]
    ylabel = paramstyle[:μ2]
    ax = Axis(fig[pos...]; xlabel, ylabel)
    hm = plot_charge_stability!(ax, data; kwargs...)
    if colorbar
        Colorbar(fig[pos[1], pos[2]+1], hm)
    end
    return ax, hm
end
function plot_charge_stability!(ax::Axis, data; datamap=x -> x.gap, colormap=:berlin, colorrange=0.5)
    x = data[:μ1]
    y = data[:μ2]
    gaps = map(datamap, data[:data])
    hm = heatmap!(ax, x, y, gaps; colormap, colorrange=(-1, 1) .* colorrange)
    return hm
end

## Plot sweet spot scan
function plot_sweet_scan(data; kwargs...)
    fig = Figure(resolution=400 .* (1.2, 1), fontsize=20, backgroundcolor=:transparent)
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
function plot_sweet_scan!(ax::Axis, data; datamap=MPU(), colormap=Reverse(:viridis), colorrange=(10.0^(-1.5), 1), colorscale=log10, plotbound=(data[:xlabel] == :U && data[:ylabel] == :h))
    x = data[:x]
    y = data[:y]
    z = map(datamap, data[:sweet_spots])
    hm = heatmap!(ax, x, y, z; colorrange, colormap, colorscale)
    if plotbound
        tratio = data[:fixedparams].tratio
        Δ = data[:fixedparams].Δ
        _bound(U) = bound(U, tratio, Δ)
        lines!(ax, x, _bound, color=:black, linestyle=:dash, linewidth=2)
        ylims!(ax, first(y), last(y))
    end
    return hm
end


# function mp_plot!(pos, x, y, data; colorrange=(-1.5, 0), cmap=default_colormap, colorscale=log10, axis=(), kwargs...)
#     mps = map(ss -> 1 - abs(ss.mpu), data)
#     ax = Axis(pos; axis...)
#     hm = heatmap!(ax, x, y, mps; colorrange, colormap=cmap, colorscale, kwargs...)
#     return ax, hm
# end
# function add_colorbar!(pos, hm; colorrange, kwargs...)
#     ticks = expticks(colorrange...)
#     cbar = Colorbar(pos, hm; ticks, kwargs...)
#     return nothing
# end
# function UVplot!(pos, data; colorrange, cmap=Reverse(:viridis), kwargs...)
#     x = data.Us
#     y = data.Vs
#     xlabel = L"U"
#     ylabel = L"V"
#     ax, hm = mp_plot!(pos, x, y, data[:sweet_spots]; cmap, colorrange, axis=(; xlabel, ylabel), kwargs...)
#     return ax, hm
# end

# function Uhplot!(pos, data; colorrange, cmap=Reverse(:viridis), kwargs...)
#     x = data.Us
#     y = sort!(unique(map(ss -> ss.parameters.h, data[:sweet_spots])))
#     xlabel = L"U"
#     ylabel = L"V_Z"
#     # bound(U) = sqrt(1 + tratio^(-2)) - U / 2# - U^2/20
#     # bound(U) = 1 / 2 * (-U + sqrt(U^2 + 4 * (-(U / tratio) + 1 + 1 / tratio^2)))
#     # bound(U) = ((-tratio * U + sqrt(4 - 4 * tratio * U + tratio^2 * (4 + U^2))) / (2 * tratio))
#     # bound(U) = sqrt(1 + tratio^2) / tratio + 1 / 2 * (-1 -  1 / sqrt(1 + tratio^2)) * U + (tratio^3 * U^2) / (8 * (1 + tratio^2)^(3 / 2))
#     # Vzmax = sqrt(1 + tratio^(-2))
#     # bound(U) = Vzmax - U/2*(1 + 1/sqrt(1+tratio^2))
#     ax, hm = mp_plot!(pos, x, y, data[:sweet_spots]; cmap, colorrange, axis=(; xlabel, ylabel), kwargs...)
#     tratio = data[:params].tratio
#     Δ = data[:params].Δ
#     _bound(U) = bound(U,tratio,Δ)
#     lines!(pos, x, _bound, color=:black, linestyle=:dash, linewidth=2)
#     ylims!(ax, first(y), last(y))
#     return ax, hm
# end

bound(U, tratio, Δ) = 1 / 2 * (-U + sqrt(
    U^2 + (4 * (-tratio * U + Δ * (1 +
                                   tratio^2))) / tratio^2))