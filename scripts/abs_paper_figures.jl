using DrWatson
@quickactivate "Majorana sweet spot"
using QuantumDots, DataFrames, CairoMakie, LaTeXStrings, Printf, QuantumDots.BlockDiagonals
using CairoMakie.Colors
using Base.Threads, LinearAlgebra
include(srcdir("misc.jl"))
includet(scriptsdir("abs_chain_misc.jl"))
resultsUV = collect_results(sdatadir("sims", "ABS-chain", "UV-scan"))
resultsUh = collect_results(sdatadir("sims", "ABS-chain", "Uh-scan"))
resultsVh = collect_results(sdatadir("sims", "ABS-chain", "Vh-scan"))
##
svg() = CairoMakie.activate!(type="svg") #For vector plots
png() = CairoMakie.activate!(type="png") #For faster plots
theme = Theme(
    c=:magma
)
##
default_colormap = :magma
default_colormap = Reverse(:tempo)
default_colormap = Reverse(:deep)
default_colormap = Reverse(:GnBu)
default_colormap = Reverse(:viridis)
distlabel = L"MP"#raw"\delta ρ_{\mu}"
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

# Makie.get_tickvalues(LLog10IntegerTicks(), 1e-1, 1)
#This function plots a heatmap of the mp
function mp_plot!(pos, x, y, data; colorrange=(-1.5, 0), cmap=default_colormap, colorscale=log10, axis=(), kwargs...)
    mps = map(ss -> 1 - abs(ss.mpu), data)
    ax = Axis(pos; axis...)
    hm = heatmap!(ax, x, y, mps; colorrange, colormap=cmap, colorscale, kwargs...)
    return ax, hm
end
function add_colorbar!(pos, hm; colorrange, kwargs...)
    ticks = expticks(colorrange...)
    cbar = Colorbar(pos, hm; ticks, kwargs...)
    return nothing
end
function UVplot!(pos, data; colorrange, cmap=Reverse(:viridis), kwargs...)
    x = data.Us
    y = data.Vs
    xlabel = L"U"
    ylabel = L"V"
    ax, hm = mp_plot!(pos, x, y, data[:sweet_spots]; cmap, colorrange, axis=(; xlabel, ylabel), kwargs...)
    return ax, hm
end
function Uhplot!(pos, data; colorrange, cmap=Reverse(:viridis), kwargs...)
    x = data.Us
    y = sort!(unique(map(ss -> ss.parameters.h, data[:sweet_spots])))
    xlabel = L"U"
    ylabel = L"V_Z"
    tratio = data[:params].tratio
    # bound(U) = sqrt(1 + tratio^(-2)) - U / 2# - U^2/20
    # bound(U) = 1 / 2 * (-U + sqrt(U^2 + 4 * (-(U / tratio) + 1 + 1 / tratio^2)))
    # bound(U) = ((-tratio * U + sqrt(4 - 4 * tratio * U + tratio^2 * (4 + U^2))) / (2 * tratio))
    # bound(U) = sqrt(1 + tratio^2) / tratio + 1 / 2 * (-1 -  1 / sqrt(1 + tratio^2)) * U + (tratio^3 * U^2) / (8 * (1 + tratio^2)^(3 / 2))
    # Vzmax = sqrt(1 + tratio^(-2))
    # bound(U) = Vzmax - U/2*(1 + 1/sqrt(1+tratio^2))
    bound(U) = 1 / 2 * (-U + sqrt(
        U^2 + (4 * (-tratio * U + 1 +
                    tratio^2)) / tratio^2))
    ax, hm = mp_plot!(pos, x, y, data[:sweet_spots]; cmap, colorrange, axis=(; xlabel, ylabel), kwargs...)
    lines!(pos, x, bound, color=:black, linestyle=:dash, linewidth=2)
    ylims!(ax, first(y), last(y))
    return ax, hm
end

##
using ColorSchemes
subtitle = L"\log{\left(1-|MP|\right)}"
fileext = ".pdf"
colorrange = (10.0^(-2), 1)
default_figure(; fontsize=20, resolution=400 .* (1.2, 1), kwargs...) = Figure(; fontsize, resolution, kwargs...)#, )
uvfig = let data = resultsUV[1, :]
    fig = default_figure()
    fig2 = default_figure(; backgroundcolor=:transparent)
    ax, hm = UVplot!(fig[1, 1], data; colorrange)
    ax, hm = UVplot!(fig2[1, 1], data; colorrange)
    add_colorbar!(fig[1, 1+end], hm; colorrange)
    add_colorbar!(fig2[1, 1+end], hm; colorrange)
    #save(splotsdir("abs-chain", "paper", string("uvplot", data[:params], fileext)), fig, pt_per_unit=1)
    #save(splotsdir("abs-chain", "paper", string("uvplot_transparent", data[:params], fileext)), fig, pt_per_unit=1)
    fig
end
uhfig = let data = resultsUh[1, :]
    fig = default_figure()
    fig2 = default_figure(; backgroundcolor=:transparent)
    ax, hm = Uhplot!(fig[1, 1], data; cmap = Reverse(:imola), colorrange)
    ax, hm = Uhplot!(fig2[1, 1], data; colorrange)
    add_colorbar!(fig[1, 1+end], hm; colorrange)
    add_colorbar!(fig2[1, 1+end], hm; colorrange)
    #save(splotsdir("abs-chain", "paper", string("uhplot", data[:params], fileext)), fig, pt_per_unit=1)
    #save(splotsdir("abs-chain", "paper", string("uhplot_transparent", data[:params], fileext)), fig2, pt_per_unit=1)
    fig
end

uhfig = Figure(fontsize=20, resolution=400 .* (1.2, 1), backgroundcolor=:transparent);
Uhplot!(uhfig, resultsUh[2, :]; cmap=:viridis)
uhfig


scanfig = let data1 = resultsUh[1, :], data2 = resultsUV[1, :]
    fig = default_figure(; resolution=400 .* (1.8, 1), backgroundcolor=:transparent, supertitle=subtitle)
    ax, hm = Uhplot!(fig[1, 1], data1; colorrange)
    ax, hm = UVplot!(fig[1, 2], data2; colorrange)
    add_colorbar!(fig[1, 1+end], hm; colorrange)
    supertitle = Label(fig[0, :], subtitle, fontsize=30)
    save(splotsdir("abs-chain", "paper", string("uhvplot_transparent", data1[:params], "_", data2[:params], fileext)), fig, pt_per_unit=1)
    fig
end


## for charge-stability, use e.g :tofino, :vanimo, :managua, :berlin, from https://s-ink.org/scientific-colour-maps 



## Check if the Kitaev V is too large. 
## It's not. So failure is due to inability to tune kitaev parameters
let data = resultsUV
    for n in axes(data, 1)[1]
        x = data[n, :Us]
        y = data[n, :Vs]
        z = map(ss -> abs(1 - ss.mp), data[n, :sweet_spots])
        fig = default_figure()
        ax = Axis(fig[1, 1])
        p = heatmap!(ax, x, y, z)#; title="$(data[n,:params]) \n LD", clims=(-5, 0), xlabel = "V", ylabel = "U")
        Colorbar(fig[1, 2], p)
        zc = map(ss -> KitaevParameters(; ss.parameters...) |> p -> -abs(p.V / 1 + 0p.t), data[n, :sweet_spots])
        contour!(ax, x, y, zc; lw=3)
        # bound(h) = sqrt(1+data[n,:params].tratio^(-2)) - h + 1
        # plot!(p, bound; xlims = (first(us),last(us)), ylims = (first(hs),last(hs)))
        fig |> display
    end
end
