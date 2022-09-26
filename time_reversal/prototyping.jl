### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ e2b065a0-1fb2-11ed-300e-d79a41c0557a
import Pkg

# ╔═╡ 65b9ceeb-38b3-496a-a10a-58d322f7a5bb
begin
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 02725372-99f4-4af9-a9e8-d81444b5af98
using SparseArrays, LinearAlgebra, Plots, KrylovKit, SpinSymmetry

# ╔═╡ 7e0de451-bb15-49da-84a0-38d575da658b
begin
	N = 8
	H = xxz(N,6)
	ψ0 = random_state(N)
	δt = 0.1
	trange = 0:δt:5
end

# ╔═╡ 570cb354-ad5f-4234-9470-1968c13fad26
begin
	struct Pulse
		n::Vector{Float64}
		α::Float64
		Ω::Float64
	end
	Pulse(α::Float64,Ω::Float64) = Pulse(Vector{Float64}([1,0,0]),α,Ω)
	struct Sequence
		pulses::Vector{Pulse}
		n_pulses::Int
		τs::Vector{Float64}
		n_τs::Int
	end
	Sequence(pulses::Vector{Pulse},τs::Vector{Float64}) = Sequence(pulses,length(pulses),τs,length(τs))
	is_valid(seq::Sequence) = (seq.n_τs == 0 || seq.n_τs == seq.n_pulses + 1) ? true : false
end

# ╔═╡ 43cd2361-c747-4309-afdc-abbbb5bc0894
A = Vector{Float64}([1,2,3])

# ╔═╡ ceedba66-3b9d-4350-bf85-606cc384ea41
Pulse(A,2.,3.)

# ╔═╡ 650eea04-2ac9-4662-b2fe-19d350d9d5bd
begin
	t=Vector{Pulse}(undef,1)
	t[1] = Pulse(A,2.,3.)
end

# ╔═╡ d989a23c-8e15-4364-9217-3a78b719bec8
t

# ╔═╡ Cell order:
# ╠═e2b065a0-1fb2-11ed-300e-d79a41c0557a
# ╠═65b9ceeb-38b3-496a-a10a-58d322f7a5bb
# ╠═02725372-99f4-4af9-a9e8-d81444b5af98
# ╠═7e0de451-bb15-49da-84a0-38d575da658b
# ╠═570cb354-ad5f-4234-9470-1968c13fad26
# ╠═43cd2361-c747-4309-afdc-abbbb5bc0894
# ╠═ceedba66-3b9d-4350-bf85-606cc384ea41
# ╠═650eea04-2ac9-4662-b2fe-19d350d9d5bd
# ╠═d989a23c-8e15-4364-9217-3a78b719bec8
