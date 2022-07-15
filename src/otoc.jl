module OTOC
using SparseArrays, LinearAlgebra, Plots
using SpinSymmetry
using KrylovKit
using ..LightCones

export  otoc, otoc_spat, krylov_from0, krylov_step


#Krylov Propagation

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

function krylov_from0_alternative(H,t,ψ,tmax=1.0)
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

#OTOCs

function otoc(H,A,B,t::Float64,ψ,tmax=1.0)
	state = B*ψ
	state = krylov_from0_alternative(H,-t,state,tmax)
	state = A*state
	state = krylov_from0_alternative(H,t,state,tmax)
	state = B*state
	state = krylov_from0_alternative(H,-t,state,tmax)
	state = A*state
	state = krylov_from0_alternative(H,t,state,tmax)
	return real(dot(ψ,state))
end

function otoc(H,A,B,trange::Union{AbstractRange{Float64},Vector{Float64}},ψ,tmax=1.0)
	res = zeros(length(trange))
	ψl_tmp = krylov_step(H,-trange[1],ψ)
	ψr_tmp = krylov_step(H,-trange[1],B*ψ)
	for (ti, t) in enumerate(trange)
		state_l = B*krylov_from0_alternative(H,t,A*ψl_tmp,tmax)
		state_r = krylov_from0_alternative(H,t,A*ψr_tmp,tmax)
		res[ti] = real(dot(state_l,state_r))
		if ti != length(trange)
			δt = trange[ti+1] - trange[ti]
			ψl_tmp = krylov_step(H,-δt,ψl_tmp)
			ψr_tmp = krylov_step(H,-δt,ψr_tmp)
		end
	end
	return res
end

###OTOCs computed over spatial indices


#Single time

function otoc_spat(H,A,b,t::Float64,ψ,N,tmax=1.0) #b in single-particle Hilbert space
	σiUψ = A * krylov_from0_alternative(H,-t,ψ,tmax)
	UdσiUψ = krylov_from0_alternative(H,t,σiUψ,tmax)
	res = zeros(N)
	Threads.@threads for j in 1:N
		B = single_spin_op(b,j,N)
		state_r = A*krylov_from0_alternative(H,-t,B*ψ,tmax)
		state_r = krylov_from0_alternative(H,t,state_r,tmax)
		state_l = B*UdσiUψ
		res[j] = real(dot(state_l,state_r))  #Note A, b self-adjoint!
	end
	return res
end

function otoc_spat(H,A,b,t::Float64,ψ,N,k,tmax=1.0) #b in single-particle Hilbert space
	σiUψ = A * krylov_from0_alternative(H,-t,ψ,tmax)
	UdσiUψ = krylov_from0_alternative(H,t,σiUψ,tmax)
	res = zeros(N)
	Threads.@threads for j in 1:N
		B = symmetrize_operator(single_spin_op(b,j,N),N,k)
		state_r = A*krylov_from0_alternative(H,-t,B*ψ,tmax)
		state_r = krylov_from0_alternative(H,t,state_r,tmax)
		state_l = B*UdσiUψ
		res[j] = real(dot(state_l,state_r))  #Note A, b self-adjoint!
	end
	return res
end

#Time Range

function calc_otoc(H,A,b,j,trange::Union{AbstractRange{Float64},Vector{Float64}},ψ,N,tmax=1.0)
	B = single_spin_op(b,j,N)
	return otoc(H,A,B,trange,ψ,tmax)
end

function calc_otoc(H,A,b,j,trange::Union{AbstractRange{Float64},Vector{Float64}},ψ,N,k,tmax=1.0)
	B = symmetrize_operator(single_spin_op(b,j,N),N,k)
	return otoc(H,A,B,trange,ψ,tmax)
end

function otoc_spat(H,A,b,trange::Union{AbstractRange{Float64},Vector{Float64}},ψ,N,tmax=1.0)
	res = zeros(length(trange),N)
	@sync for j in 1:N
		Threads.@spawn res[:,j]=calc_otoc(H,A,b,j,trange,ψ,N,tmax)
	end
	return res
end

function otoc_spat(H,A,b,trange::Union{AbstractRange{Float64},Vector{Float64}},ψ,N,tmax=1.0)
	res = zeros(length(trange),N)
	@sync for j in 1:N
		Threads.@spawn res[:,j]=calc_otoc(H,A,b,j,trange,ψ,N,k,tmax)
	end
	return res
end

function otoc_spat!(res,H,A,b,trange::Union{AbstractRange{Float64},Vector{Float64}},ψ,N,tmax=1.0)
	@sync for j in 1:N
		Threads.@spawn res[:,j]=calc_otoc(H,A,b,j,trange,ψ,N,tmax)
	end
	return res
end

function otoc_spat!(res,H,A,b,trange::Union{AbstractRange{Float64},Vector{Float64}},ψ,N,k,tmax=1.0)
	@sync for j in 1:N
		Threads.@spawn res[:,j]=calc_otoc(H,A,b,j,trange,ψ,N,k,tmax)
	end
	return res
end

#####################################
###Old implementations for comparison


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

function otoc_spat_old(H,opi,opj,i,t,ψ,N) #opj in single-particle Hilbert space
	σiUψ = opi * exponentiate(H,-im*t,ψ)[1]
	UσiUψ = exponentiate(H,im*t,σiUψ)[1]
	C = zeros(N)
	Threads.@threads for j in 1:N
		single_spin_opj = single_spin_op(opj,j,N)
		state_r = opi*exponentiate(H,-im*t,single_spin_opj*ψ)[1]
		state_r = exponentiate(H,im*t,state_r)[1]
		state_l = single_spin_opj*UσiUψ
		C[j] = real(dot(state_l,state_r))  #Note opi, opj self-adjoint!
	end
	return C
end

end #module