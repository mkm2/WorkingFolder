### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 6e7494e2-84f6-11ec-1055-3dcdceddff31
import Pkg

# ╔═╡ f7bbc64a-e3cc-407b-b419-98e09d38db4e
begin
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ fc58fc6b-893e-4ffc-b5c3-634784b2a564
using SparseArrays, LinearAlgebra, Plots, KrylovKit, SpinSymmetry

# ╔═╡ 406b55d9-1dc8-45f2-b31e-b3c1b455d346
begin
	N = 8
	H = xxz(N,6)
	ψ0 = normalize!(ones(2^N))
	δt = 0.1
	trange = 0:δt:5
end

# ╔═╡ d1a45b7d-02f1-4b71-98e0-1817766bfb7d
function otoc_old(H,A,B,t,ψ)
	state = B*ψ
	state = exponentiate(H,-im*t,state)[1]
	state = A*state
	state = exponentiate(H,im*t,state)[1]
	state = B*state
	state = exponentiate(H,-im*t,state)[1]
	state = A*state
	state = exponentiate(H,im*t,state)[1]
	return real(dot(ψ,state))
end

# ╔═╡ dc33c272-3ca2-4785-ad3d-cb13a06cd72a
begin
	function krylov_step(H,δt,ψ)
		return exponentiate(H,im*δt,ψ;ishermitian=true)[1]
	end
	
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
end

# ╔═╡ 18532a7a-6a9d-4a2a-9e67-92225e4e1332
function incremental_krylovkit(H, ψ0, ts)
	times = zeros(length(ts))
	ret = Vector{Vector{ComplexF64}}(undef, length(ts))
	times[1] = @elapsed ret[1] = KrylovKit.exponentiate(H, -im*ts[1], ψ0)[1]
	δts = ts[2:end] .- ts[1:end-1]
	for (i,δt) in enumerate(δts)
		times[i+1] = @elapsed ret[i+1] = KrylovKit.exponentiate(H, -im*δt, ret[i])[1]
	end
	ret, times
end

# ╔═╡ ee29bd24-a780-4ad9-9093-5a2d42240936
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

# ╔═╡ 9f36b726-5b16-457e-827b-d67275389a25
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



# ╔═╡ 545941f5-0883-4f0d-bb59-cf519518b81d
begin
	op1 = single_spin_op(σz,5,N)
	op2 = single_spin_op(σz,1,N)
end

# ╔═╡ 612b6190-aa21-4b12-ae52-031996f99ae0
begin
	corr = zeros(Float64,length(trange))
	@elapsed for (ti,t) in enumerate(trange)
	    corr[ti] = 2-2*otoc_old(H, op1, op2, t, ψ0)
	end
end

# ╔═╡ a004f30f-f5fe-4e78-88eb-4ac59b392f4f
begin
	corr2 = zeros(Float64,length(trange))
	@elapsed for (ti,t) in enumerate(trange)
	    corr2[ti] = 2-2*otoc(H, op1, op2, t, ψ0)
	end
end

# ╔═╡ 3201ef85-9f9f-43f6-b8f1-833a1b1e6d79
@elapsed corr3=2*ones(length(trange))-2*otoc(H,op1,op2,trange,ψ0,δt)

# ╔═╡ 54717cba-349e-41cd-a574-4ae0561e697e
norm(corr3-corr)

# ╔═╡ c8757717-6a9b-450d-b1b6-c364b7a1b94c
begin
	@info @elapsed otoc(H,op1,op2,6.0,ψ0,δt)
	@info @elapsed otoc_old(H,op1,op2,6,ψ0)
	norm(otoc(H,op1,op2,6.0,ψ0,δt)-otoc_old(H,op1,op2,6,ψ0))
end

# ╔═╡ 5b649dce-d4c8-46de-8a0f-7c8c87e3285a
begin
	@info @elapsed krylov_from0(H,-16,ψ0,δt)
	@info @elapsed exponentiate(H,-im*16,ψ0)[1]
	norm(krylov_from0(H,6,ψ0,δt)-exponentiate(H,im*6,ψ0)[1])
end

# ╔═╡ 006eb8c7-2373-4c1f-ba41-6674234787bb
md"otoc\_old läuft mit ``O_0 = 4L_t x0``.

otoc(t) läuft mit ``O_0 \cdot 1.5 = 4L_tx``.

otoc(trange) spart Zeit und läuft mit ``O_0 \cdot 0.8 = 2L_t(x+y)``."

# ╔═╡ 240b710c-1a3a-4487-bcd9-3a28e3c871bd
md"## Spatial - Try different implementations to beat old one"

# ╔═╡ 8804cb51-a24a-4bf3-bdd7-84a0ad9a6986
function otoc_spat_old(H,opi,opj,i,t,ψ,N) #opj in single-particle Hilbert space
	σiUψ = opi * exponentiate(H,-im*t,ψ)[1]
	UdσiUψ = exponentiate(H,im*t,σiUψ)[1]
	C = zeros(N)
	Threads.@threads for j in 1:N
		single_spin_opj = single_spin_op(opj,j,N)
		state_r = opi*exponentiate(H,-im*t,single_spin_opj*ψ)[1]
		state_r = exponentiate(H,im*t,state_r)[1]
		state_l = single_spin_opj*UdσiUψ
		C[j] = real(dot(state_l,state_r))  #Note opi, opj self-adjoint!
	end
	return C
end

# ╔═╡ c249e433-7e2a-44df-88f6-f126d520987c
begin
	function otoc_spat1(H,opi,opj,t::Float64,ψ,N,δt=0.1) #opj in single-particle Hilbert space
		σiUψ = opi * krylov_from0(H,-t,ψ,δt)
		UdσiUψ = krylov_from0(H,t,σiUψ,δt)
		res = zeros(N)
		Threads.@threads for j in 1:N
			single_spin_opj = single_spin_op(opj,j,N)
			state_r = opi*krylov_from0(H,-t,single_spin_opj*ψ,δt)
			state_r = krylov_from0(H,t,state_r,δt)
			state_l = single_spin_opj*UdσiUψ
			res[j] = real(dot(state_l,state_r))  #Note opi, opj self-adjoint!
		end
		return res
	end
	
	function otoc_spat2(H,opi,opj,trange::AbstractRange{Float64},ψ,N,δt=0.1) #opj in single-particle Hilbert space
		res = zeros(length(trange),N)
		ψl_tmp_orig = krylov_step(H,-trange[1],ψ)
		σiψl_tmp_orig = opi*ψl_tmp_orig
		Threads.@threads for j in 1:N
			ψl_tmp = ψl_tmp_orig
			σiψl_tmp = σiψl_tmp_orig
			single_spin_opj = single_spin_op(opj,j,N)
			ψr_tmp = krylov_step(H,-trange[1],single_spin_opj * ψ)
			for (ti,t) in collect(enumerate(trange))
				state_l = single_spin_opj * krylov_from0(H,t,σiψl_tmp,δt)
				state_r = krylov_from0(H,t,opi*ψr_tmp,δt)
				res[ti,j] = real(dot(state_l,state_r))
				if ti!=length(trange)
					ψl_tmp = krylov_step(H,-δt,ψl_tmp)
					σiψl_tmp = opi*ψl_tmp
					ψr_tmp = krylov_step(H,-δt,ψr_tmp)
				end
			end
		end
		return res
	end
end

# ╔═╡ 8efb255b-6fdf-4326-8895-86970a79fe32
function otoc_spat3(H,opi,opj,trange::AbstractRange{Float64},ψ,N,δt=0.1)
	res = zeros(N,length(trange))
	ψl_tmp = krylov_step(H,-trange[1],ψ)
	σiψl_tmp = opi*ψl_tmp
	ψrs_tmp = zeros(ComplexF64,2^N,N)
	for j in 1:N
		single_spin_opj = single_spin_op(opj,j,N)
		ψrs_tmp[:,j] = krylov_step(H,-trange[1],single_spin_opj*ψ)
	end
	for (ti,t) in enumerate(trange)
		Threads.@threads for j in 1:N
			single_spin_opj = single_spin_op(opj,j,N)
			state_l = single_spin_opj * krylov_from0(H,t,σiψl_tmp,δt)
			state_r = krylov_from0(H,t,opi*ψrs_tmp[:,j],δt)
			res[j,ti] = real(dot(state_l,state_r))
		end
		if ti!=length(trange)
			print(ψl_tmp)
			ψl_tmp = krylov_step(H,-δt,ψl_tmp)
			σiψl_tmp = opi*ψl_tmp
			for j in 1:N
				ψrs_tmp[:,j] = krylov_step(H,-δt,ψrs_tmp[:,j])
			end
		end
	end
	return res
end

# ╔═╡ bac27bc1-ca34-4e7f-b55d-d5a50b0b1ec3
function otoc_spat4(H,opi,opj,trange::AbstractRange{Float64},ψ,N,δt=0.1)
	res = zeros(length(trange),N)
	Threads.@threads for j in 1:N
		single_spin_opj = single_spin_op(opj,j,N)
		res[:,j]=otoc(H,opi,single_spin_opj,trange,ψ0,δt)
	end
	return res
end

# ╔═╡ 8a70b208-fe35-4da7-b658-5d13a7b66b64
function otoc_spat_jfixed(H,opi,opj,trange::AbstractRange{Float64},ψ,N,δt=0.1)
	res = zeros(N,length(trange))
	ψr_tmp = krylov_step(H,-trange[1],opj * ψ)
	ψl_tmp = krylov_step(H,-trange[1],ψ)
	for (ti,t) in collect(enumerate(trange))
		Threads.@threads for i in 1:N
			single_spin_opi = single_spin_op(opi,i,N)
			state_r = krylov_from0(H,t,single_spin_opi * ψr_tmp,δt)
			state_l = opj * krylov_from0(H,t,single_spin_opi * ψl_tmp,δt)
			res[i,ti] = real(dot(state_l,state_r))
		end
		if ti!=length(trange)
			ψr_tmp = krylov_step(H,-δt,ψr_tmp)
			ψl_tmp = krylov_step(H,-δt,ψl_tmp)
		end
	end
	return res
end

# ╔═╡ 71ebadc4-e4ec-49c9-84a6-030f2aa345c2
begin
	i = 3
	σzi = single_spin_op(σz,i,N)
	trange2 = 0:δt:6
	otocs_old = zeros(length(trange2),N)
	@elapsed for (ti,t) in collect(enumerate(trange2))
		otocs_old[ti,:] = otoc_spat_old(H,σzi,σz,i,t,ψ0,N)
	end
end

# ╔═╡ 097a013f-0e6b-4e44-81bd-73f0e5a447fc
begin
	otocs1 = zeros(length(trange2),N)
	@elapsed for (ti,t) in collect(enumerate(trange2))
		otocs1[ti,:] = otoc_spat1(H,σzi,σz,t,ψ0,N,δt)
	end
	#corrs = 2*ones(length(trange),N)-2*otocs
	#heatmap(1:N, trange, corrs, c = :viridis)
end

# ╔═╡ e324f7a3-0930-4ade-9652-d57731278a47
begin
	otocs2 = zeros(length(trange2),N)
	@elapsed otocs2 = otoc_spat2(H,σzi,σz,trange2,ψ0,N,δt)
end

# ╔═╡ f68a0cd9-544c-402e-a962-8730531fd847
begin
	otocs3 = zeros(N,length(trange2))
	@elapsed otocs3 = otoc_spat3(H,σzi,σz,trange2,ψ0,N,δt)
end

# ╔═╡ 3e538a4a-004b-4a1c-8ac5-f49e371faad7
begin
	otocs4 = zeros(length(trange2),N)
	@elapsed otocs4 = otoc_spat4(H,σzi,σz,trange2,ψ0,N,δt)
end

# ╔═╡ 61cfa3a5-2616-4c28-8c17-f669fd82865f
begin
	σzj = single_spin_op(σz,3,N)
	otocs5 = zeros(N,length(trange2))
	@elapsed otocs5 = otoc_spat_jfixed(H,σz,σzj,trange2,ψ0,N,δt)
end

# ╔═╡ 7b38ff87-b919-4904-a282-b2122d39407d
norm(otocs_old-otocs4)

# ╔═╡ b27e90ec-f64d-491a-8a33-f6a993bc0113
md"Vor Threading:

otoc\_spat\_old läuft mit ``O_0/2*(N+1)``.

otoc\_spat1 mit ``O_0^s \cdot 1.5``.

otoc\_spat2 mit ``O_0^s \cdot 1.5 + 2L_t (-x+Ny) + (N+1)y``.

otoc\_spat3 mit ``O_0^s \cdot 1.5 + 2L_t (-x+\frac{N+1}{2}y) + (N+1)y``.

otoc\_spat4 mit ``N O_0 \cdot 0.8``, was für ``N>2`` langsamer als otoc\_spat\_old ist.

otoc\_spat\_fixedj läuft mit ``O_0 \cdot 1.5 + 2L_t (-x+y)``.

Selbst mit der optimistischen Abschätzung ``2L_tx \approx 0.75 O_0 `` kommt man zu keiner Ersparnis, sobald ``N>2``."

# ╔═╡ e3c5532b-f57b-448d-91d7-4d0f8eb3bc29
md"Nach Threading:
otoc\_spat\_old ist am schnellsten, aber otoc\_spat2 und otoc\_spat4 sind nur wenig langsamer, Faktor Runtime ca. 1.05 - 1.2" 

# ╔═╡ Cell order:
# ╠═6e7494e2-84f6-11ec-1055-3dcdceddff31
# ╠═f7bbc64a-e3cc-407b-b419-98e09d38db4e
# ╠═fc58fc6b-893e-4ffc-b5c3-634784b2a564
# ╠═406b55d9-1dc8-45f2-b31e-b3c1b455d346
# ╠═d1a45b7d-02f1-4b71-98e0-1817766bfb7d
# ╠═dc33c272-3ca2-4785-ad3d-cb13a06cd72a
# ╠═18532a7a-6a9d-4a2a-9e67-92225e4e1332
# ╠═ee29bd24-a780-4ad9-9093-5a2d42240936
# ╠═9f36b726-5b16-457e-827b-d67275389a25
# ╠═545941f5-0883-4f0d-bb59-cf519518b81d
# ╠═612b6190-aa21-4b12-ae52-031996f99ae0
# ╠═a004f30f-f5fe-4e78-88eb-4ac59b392f4f
# ╠═3201ef85-9f9f-43f6-b8f1-833a1b1e6d79
# ╠═54717cba-349e-41cd-a574-4ae0561e697e
# ╠═c8757717-6a9b-450d-b1b6-c364b7a1b94c
# ╠═5b649dce-d4c8-46de-8a0f-7c8c87e3285a
# ╠═006eb8c7-2373-4c1f-ba41-6674234787bb
# ╠═240b710c-1a3a-4487-bcd9-3a28e3c871bd
# ╠═8804cb51-a24a-4bf3-bdd7-84a0ad9a6986
# ╠═c249e433-7e2a-44df-88f6-f126d520987c
# ╠═8efb255b-6fdf-4326-8895-86970a79fe32
# ╠═bac27bc1-ca34-4e7f-b55d-d5a50b0b1ec3
# ╠═8a70b208-fe35-4da7-b658-5d13a7b66b64
# ╠═71ebadc4-e4ec-49c9-84a6-030f2aa345c2
# ╠═097a013f-0e6b-4e44-81bd-73f0e5a447fc
# ╠═e324f7a3-0930-4ade-9652-d57731278a47
# ╠═f68a0cd9-544c-402e-a962-8730531fd847
# ╠═3e538a4a-004b-4a1c-8ac5-f49e371faad7
# ╠═61cfa3a5-2616-4c28-8c17-f669fd82865f
# ╠═7b38ff87-b919-4904-a282-b2122d39407d
# ╠═b27e90ec-f64d-491a-8a33-f6a993bc0113
# ╠═e3c5532b-f57b-448d-91d7-4d0f8eb3bc29
