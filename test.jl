### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ ef776e9a-8c59-11ec-3b85-d712121bcb49
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ 8e3e79cf-6f64-4a57-bf54-064f6a13aa38
begin
	import Dates
	using SpinSymmetry, Random, LinearAlgebra
	using LightCones
	using KrylovKit
end

# ╔═╡ 5956500c-1878-436a-8beb-5aa3b5b71828
function krylov_step(H,δt,ψ)
	return exponentiate(H,im*δt,ψ;ishermitian=true)[1]
end

# ╔═╡ 8bc73f9b-8692-4f01-9850-0f3f913c1013
function krylov_from0(H,t,ψ,δt) 
	#Assume t multiple of timestep
	T = abs(round(t/δt,digits=0))
	sgn = t >= 0 ? 1.0 : -1.0
	N_unit_steps = floor(abs(t))
	δt_unit = sgn * 1.0
	for i in 1:N_unit_steps
		ψ = krylov_step(H,δt_unit,ψ)
	end
	N_short_steps = round((abs(t)-N_unit_steps)/δt,digits=0)
	δt_short = sgn*δt
	for i in 1:N_short_steps
		ψ = krylov_step(H,δt_short,ψ)
	end
	return ψ
end

# ╔═╡ e6b85e66-5364-43ca-941d-22e645c6bf1f
function otoc(H,A,B,t::Float64,ψ,δt=0.1)
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

# ╔═╡ 87080c9f-2b6a-430b-bcca-fee9438a8331
function otoc(H,A,B,trange::AbstractRange{Float64},ψ,δt=0.1)
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

# ╔═╡ 58d80e4d-2d21-4bf5-9660-aa0c2cc7d46c
function otoc_spat(H,opi,opj,trange::AbstractRange{Float64},ψ,N,δt=0.1)
	res = zeros(length(trange),N)
	Threads.@threads for j in 1:N
		single_spin_opj = single_spin_op(opj,j,N)
		res[:,j]=otoc(H,opi,single_spin_opj,trange,ψ,δt)
	end
	return res
end

# ╔═╡ 3d06f1d1-f7e8-404b-a769-809556985ebc
begin
	N = 12
	δt = 0.1
	H = xxz(N,6)
	ψ0 = normalize!(ones(2^N))
	
	op1 = single_spin_op(σz,5,N)
	op2 = single_spin_op(σz,1,N)
	
	trange = 0:δt:5
	
	corr = zeros(Float64,length(trange))
end

# ╔═╡ a1ab5e37-90e8-416f-8621-3d23f5ed9677
@time Threads.@threads for (ti,t) in collect(enumerate(trange))
    corr[ti] = 2-2*otoc(H, op1, op2, t, ψ0)
end

# ╔═╡ 1a8ffbd8-5201-40b7-bc9a-872faabefff0
begin
	i = 3
	σzi = single_spin_op(σz,i,N)
	
	otocs2 = zeros(length(trange),N)
	otocs2 = otoc_spat(H,σzi,σz,trange,ψ0,N,δt)
end

# ╔═╡ Cell order:
# ╠═ef776e9a-8c59-11ec-3b85-d712121bcb49
# ╠═8e3e79cf-6f64-4a57-bf54-064f6a13aa38
# ╠═5956500c-1878-436a-8beb-5aa3b5b71828
# ╠═8bc73f9b-8692-4f01-9850-0f3f913c1013
# ╠═e6b85e66-5364-43ca-941d-22e645c6bf1f
# ╠═58d80e4d-2d21-4bf5-9660-aa0c2cc7d46c
# ╠═87080c9f-2b6a-430b-bcca-fee9438a8331
# ╠═3d06f1d1-f7e8-404b-a769-809556985ebc
# ╠═a1ab5e37-90e8-416f-8621-3d23f5ed9677
# ╠═1a8ffbd8-5201-40b7-bc9a-872faabefff0
