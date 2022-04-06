### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 98f1970c-b5a1-11ec-2585-8de5d582401d
using SparseArrays, LinearAlgebra, Plots, SpinSymmetry, KrylovKit, Pkg

# ╔═╡ b07f37db-cfb4-4db1-87d3-dc008f6eebdd
begin
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 04017e18-de2e-4b50-a42f-e9fa716632bc
function krylov_from0_alternative(H,t,ψ,tmax=4)
	N_max_steps = floor(abs(t)/tmax)
	sgn = t >= 0 ? 1.0 : -1.0
	δt = sgn * tmax
	for i in 1:N_max_steps
		ψ = krylov_step(H,δt,ψ)
	end
	t_res = t - N_max_steps * δt #signed
	ψ = krylov_step(H,t_res,ψ)
	return ψ
end

# ╔═╡ c31e203e-6265-4082-8cd9-cf692cb72e5b
function otoc_old(H,A,B,t::Float64,ψ,δt=0.1)
	state = B*ψ
	state = krylov_from0(H,-t,state,δt)
	state = A*state
	state = krylov_from0(H,t,state,δt)
	state = B*state
	state = krylov_from0(H,-t,state,δt)
	state = A*state
	state = krylov_from0(H,t,state,δt)
	return real(dot(ψ,state))
end

# ╔═╡ 85705744-2608-4535-ba7b-1d8601e85abf
function otoc_old(H,A,B,trange::AbstractRange{Float64},ψ,δt=0.1)
	res = zeros(length(trange))
	ψl_tmp = krylov_step(H,-trange[1],ψ)
	ψr_tmp = krylov_step(H,-trange[1],B*ψ)
	for (ti, t) in enumerate(trange)
		state_l = B*krylov_from0(H,t,A*ψl_tmp,δt)
		state_r = krylov_from0(H,t,A*ψr_tmp,δt)
		res[ti] = real(dot(state_l,state_r))
		if ti != length(trange)
			ψl_tmp = krylov_step(H,-δt,ψl_tmp)
			ψr_tmp = krylov_step(H,-δt,ψr_tmp)
		end
	end
	return res
end

# ╔═╡ 9f153f87-a0e4-4c91-aa8d-bc6061cd8e63
begin
	N = 10
	H = xxz(N,6) + field_term(2.0,N)
	ψ0 = normalize!(ones(2^N))
	δt = 0.01
	trange = 0:δt:2
	i = div(N,2)
	A = single_spin_op(σz,i,N)
	B = single_spin_op(σz,2,N)
	t = 5.5
end

# ╔═╡ 23b3d522-972c-4c35-9fff-51e913382e2a
md"## krylov_from0"

# ╔═╡ f41553ef-f3d3-4c27-a0cd-1aae4b1a9e6b
@elapsed krylov_from0_alternative(H,t,ψ0,1) #faster as δt = 0.1

# ╔═╡ 506c0646-9675-4693-9567-27197a3c50d5
@elapsed krylov_from0(H,t,ψ0,δt) #much slower for δt = 0.01

# ╔═╡ 3b0cc39f-c0b4-4322-bb0a-ef42dd29b731
norm(krylov_from0_alternative(H,t,ψ0,1)-krylov_from0(H,t,ψ0,δt))

# ╔═╡ 66593664-6936-4274-9cb3-cedbc0e92954
md"## otoc(t)"

# ╔═╡ 1f8739f1-66f6-4381-8df9-e5738b8fa0b2
#faster now that tmax is fixed!

# ╔═╡ adf626a2-12eb-412a-834a-80a85271aca9
@elapsed otoc(H,A,B,t,ψ0,1)

# ╔═╡ 2481fba9-8779-459e-a656-40114a4b42a3
@elapsed otoc_old(H,A,B,t,ψ0,δt)

# ╔═╡ a9fbfc0f-f37b-40b5-ad50-2899195aaea8
norm(otoc(H,A,B,t,ψ0,1)-otoc_old(H,A,B,t,ψ0,δt))

# ╔═╡ 1376e368-fd6d-473f-a5c2-7f53f51ec931
md"## otoc(trange)"

# ╔═╡ d3b9be9d-de6a-4f5e-871f-4e550b2f33ea
#faster again! Almost half the runtime for δt = 0.1 but that highly depends on disorder strength! Factor of 10 compared to δt = 0.01

# ╔═╡ 3dd7686e-282f-4171-9459-87db1a35b6cc
@elapsed otoc(H,A,B,trange,ψ0,1)

# ╔═╡ 1b130302-17b2-4ea7-abce-35992b79f1db
@elapsed otoc_old(H,A,B,trange,ψ0,δt*10)

# ╔═╡ 6ee78e3f-bee0-4f36-9279-853951b04f8e
# ╠═╡ disabled = true
#=╠═╡
norm(otoc(H,A,B,trange,ψ0)-otoc_old(H,A,B,trange,ψ0,δt))
  ╠═╡ =#

# ╔═╡ 50915b7a-5397-4c5c-8b87-d8e340660df9
1422.89495826/182.483490198

# ╔═╡ 0d98ef57-669f-4654-9f05-e26ea055f60a


# ╔═╡ Cell order:
# ╠═98f1970c-b5a1-11ec-2585-8de5d582401d
# ╠═b07f37db-cfb4-4db1-87d3-dc008f6eebdd
# ╠═04017e18-de2e-4b50-a42f-e9fa716632bc
# ╠═c31e203e-6265-4082-8cd9-cf692cb72e5b
# ╠═85705744-2608-4535-ba7b-1d8601e85abf
# ╠═9f153f87-a0e4-4c91-aa8d-bc6061cd8e63
# ╠═23b3d522-972c-4c35-9fff-51e913382e2a
# ╠═f41553ef-f3d3-4c27-a0cd-1aae4b1a9e6b
# ╠═506c0646-9675-4693-9567-27197a3c50d5
# ╠═3b0cc39f-c0b4-4322-bb0a-ef42dd29b731
# ╠═66593664-6936-4274-9cb3-cedbc0e92954
# ╠═1f8739f1-66f6-4381-8df9-e5738b8fa0b2
# ╠═adf626a2-12eb-412a-834a-80a85271aca9
# ╠═2481fba9-8779-459e-a656-40114a4b42a3
# ╠═a9fbfc0f-f37b-40b5-ad50-2899195aaea8
# ╠═1376e368-fd6d-473f-a5c2-7f53f51ec931
# ╠═d3b9be9d-de6a-4f5e-871f-4e550b2f33ea
# ╠═3dd7686e-282f-4171-9459-87db1a35b6cc
# ╠═1b130302-17b2-4ea7-abce-35992b79f1db
# ╠═6ee78e3f-bee0-4f36-9279-853951b04f8e
# ╠═50915b7a-5397-4c5c-8b87-d8e340660df9
# ╠═0d98ef57-669f-4654-9f05-e26ea055f60a
