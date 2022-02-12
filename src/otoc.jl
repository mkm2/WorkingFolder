module OTOC

using SparseArrays, LinearAlgebra, Plots
using SpinSymmetry
using KrylovKit
using ..LightCones

export  otoc, otoc_spat, krylov_from0, krylov_step


#Krylov

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



#OTOCs

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

###Spatial OTOCs

function otoc_spat(H,opi,opj,t::Float64,ψ,N,δt=0.1) #opj in single-particle Hilbert space
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

function otoc_spat(H,opi,opj,trange::AbstractRange{Float64},ψ,N,δt=0.1)
	res = zeros(length(trange),N)
	Threads.@threads for j in 1:N
		single_spin_opj = single_spin_op(opj,j,N)
		res[:,j]=otoc(H,opi,single_spin_opj,trange,ψ0,δt)
	end
	return res
end

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