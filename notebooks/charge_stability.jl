### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 82b99d8b-38ac-4dc8-8234-c96d00679718
using DrWatson

# ╔═╡ 90469b5e-c800-41fc-961f-2c314e5fdd3f
begin
	do_first = true
	@quickactivate "ABS-Majoranas"
end

# ╔═╡ b20762b2-9220-4441-94d2-7c07c3f79e79
begin
	do_first
	using CairoMakie, QuantumDots, PlutoUI
end

# ╔═╡ bfef5612-f803-443a-a694-74a6959e75e4
begin
	include(srcdir("abs_chain_misc.jl"))
	include(srcdir("plotting.jl"))
end

# ╔═╡ d4634301-c96d-4ac4-ae6a-cd7b98eb4930
html"""
<script>
	const button = document.createElement("button")

	button.addEventListener("click", () => {
		editor_state_set(old_state => ({
			notebook: {
				...old_state.notebook,
				process_status: "no_process",
			},
		})).then(() => {
			window.requestAnimationFrame(() => {
				document.querySelector("#process_status a").click()
			})
		})
	})
	button.innerText = "Restart process"

	return button
</script>
"""

# ╔═╡ d9014a00-5e58-4121-a11d-194f0b1499e4
function get_colorrange(data, datamap)
	maximum(abs ∘ datamap, data)
end

# ╔═╡ 571fa1cd-4109-4cc9-8ef5-cb6be8045feb
md"""
Transport $(@bind do_transport CheckBox())
"""

# ╔═╡ 2b86a979-7812-4f38-ac65-193339b3b69c
md"""
dμ $(@bind dμ PlutoUI.Slider(.1:0.1:10, default=5, show_value = true))
res $(@bind res PlutoUI.Slider(5:100, default=20, show_value = true))

μ1 $(@bind μ1 PlutoUI.Slider(-10:0.1:10, default=0, show_value = true)) 
μ2 $(@bind μ2 PlutoUI.Slider(-10:0.1:10, default=0, show_value = true)) 
"""

# ╔═╡ ced9eac8-9d15-49c4-be37-a1b96da5ea40
md"""
Δ $(@bind Δ PlutoUI.Slider(.1:0.1:2, default=1, show_value = true))
t $(@bind t PlutoUI.Slider(.1:0.1:2, default=.5, show_value = true)) 

tratio $(@bind tratio PlutoUI.Slider(.1:0.1:2, default=.2, show_value = true)) 
U $(@bind U PlutoUI.Slider(0:0.1:10, default=0, show_value = true)) 

Unl $(@bind V PlutoUI.Slider(0:0.1:2, default=0, show_value = true)) 
Vz $(@bind h PlutoUI.Slider(0:0.1:10, default=1.5, show_value = true)) 

ϕ $(@bind ϕ PlutoUI.Slider(0:0.1:pi, default=pi/2, show_value = true)) 
"""

# ╔═╡ e3f535f9-7526-453c-b86e-634cad32a91b
if do_transport md"""
T^(-1) $(@bind Tinv PlutoUI.Slider(1:100, default=40, show_value = true))
μL $(@bind μL PlutoUI.Slider(-5:.1:5, default=0, show_value = true))
μR $(@bind μR PlutoUI.Slider(-5:.1:5, default=0, show_value = true))
"""
end

# ╔═╡ 076d480d-dc86-49f5-a71b-82f7ccf6ec81
transport = do_transport ? Transport(QuantumDots.Pauli(), (; T=Tinv^-1, μ=(μL, μR))) : missing
 

# ╔═╡ a37b08ef-eb29-45cc-be00-0feb008d69db
begin
	csdata = charge_stability_scan((; tratio, t,Δ, ϕ, U, V, μ1, μ2, h), dμ, dμ, res; transport);
	"Calc data"
end

# ╔═╡ 6f0650d3-d6b0-4f74-976b-f17de02b45b7
begin
	csdata;
	fig =  do_transport ? Figure(resolution = 600 .* (1,1.2)) : Figure()
	"Create fig"
end

# ╔═╡ 8752401c-5de3-452a-b6e8-80a3a2ba1a84
begin
	ax1 = Axis(fig[1,1]; subtitle = "Gap")
	ax2 = Axis(fig[1,2]; subtitle = "Parity")
	ax3 = Axis(fig[2,1]; subtitle = "LD")
	ax4 = Axis(fig[2,2]; subtitle = "MPU")
	plot_charge_stability!(ax1, csdata)
	plot_charge_stability!(ax2, csdata, datamap = x->sign(x.gap))
	plot_charge_stability!(ax3, csdata; colormap = Reverse(:viridis), colorrange =(0,1), datamap = LD)
	plot_charge_stability!(ax4, csdata; colormap = Reverse(:viridis), colorrange =(0,1), datamap = MPU)
	if do_transport
		ax5 = Axis(fig[3,1]; subtitle = "left conductance")
		ax6 = Axis(fig[3,2]; subtitle = "left Non-local conductance")
		ax7 = Axis(fig[4,1]; subtitle = "right conductance")
		ax8 = Axis(fig[4,2]; subtitle = "right Non-local conductance")
		data5 = x->real(x.conductance[1,1])
		data6 = x->real(x.conductance[1,2])
		data7 = x->real(x.conductance[2,2])
		data8 = x->real(x.conductance[2,1])
		cr5 = get_colorrange(csdata[:data],data5)
		cr6 = get_colorrange(csdata[:data],data6)
		cr7 = get_colorrange(csdata[:data],data7)
		cr8 = get_colorrange(csdata[:data],data8)
		colormap = :vik
		plot_charge_stability!(ax5, csdata; colormap, colorrange = cr5, datamap = data5)
		plot_charge_stability!(ax6, csdata; colormap, colorrange = cr6, datamap = data6)
		plot_charge_stability!(ax7, csdata; colormap, colorrange = cr7, datamap = data7)
		plot_charge_stability!(ax8, csdata; colormap, colorrange = cr8, datamap = data8)
	end
	"Plotting"
end

# ╔═╡ c9aa968d-838c-40cb-92c9-44ca5b95f3bf
begin
	csdata;
	transport;
	ax1;
	fig
end

# ╔═╡ Cell order:
# ╠═82b99d8b-38ac-4dc8-8234-c96d00679718
# ╠═90469b5e-c800-41fc-961f-2c314e5fdd3f
# ╠═b20762b2-9220-4441-94d2-7c07c3f79e79
# ╠═bfef5612-f803-443a-a694-74a6959e75e4
# ╟─d4634301-c96d-4ac4-ae6a-cd7b98eb4930
# ╠═a37b08ef-eb29-45cc-be00-0feb008d69db
# ╠═6f0650d3-d6b0-4f74-976b-f17de02b45b7
# ╟─d9014a00-5e58-4121-a11d-194f0b1499e4
# ╠═8752401c-5de3-452a-b6e8-80a3a2ba1a84
# ╟─076d480d-dc86-49f5-a71b-82f7ccf6ec81
# ╟─571fa1cd-4109-4cc9-8ef5-cb6be8045feb
# ╟─2b86a979-7812-4f38-ac65-193339b3b69c
# ╟─ced9eac8-9d15-49c4-be37-a1b96da5ea40
# ╟─c9aa968d-838c-40cb-92c9-44ca5b95f3bf
# ╟─e3f535f9-7526-453c-b86e-634cad32a91b
