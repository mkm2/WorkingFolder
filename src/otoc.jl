module OTOC
using SparseArrays, LinearAlgebra, Plots
using SpinSymmetry, Statistics
using KrylovKit
using ..LightCones

export  otoc, otoc_spat, krylov_from0, krylov_step, krylov_from0_alternative
export otoc_ed, otoc_edtr, otoc_edψ, otoc_spat_ed, otoc_spat_edtr, otoc_spat_edψ
export Diag_OTOC, Diag_OTOCψ, Diag_OTOCtr, DiagOTOC_multiplestates


ExtRange = Union{AbstractRange{Float64},Vector{Float64}}
TvExtRange = Union{Float64,ExtRange}
BMat = Union{Matrix{ComplexF64},SparseMatrixCSC{ComplexF64}}

##########################
### Krylov Propagation ###
##########################


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

function otoc(H,A,B,trange::ExtRange,ψ,tmax=1.0)
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

###########################################################
### Krylov OTOCs computed over multiple spatial indices ###
###########################################################

#Single time

function otoc_spat(H,A,b,t::Float64,ψ,N,tmax) #b in single-particle Hilbert space
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

function otoc_spat(H,A,b,t::Float64,ψ,N,k,tmax) #b in single-particle Hilbert space
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

function calc_otoc(H,A,b,j,trange::ExtRange,ψ,N,tmax)
	B = single_spin_op(b,j,N)
	return otoc(H,A,B,trange,ψ,tmax)
end

function calc_otoc(H,A,b,j,trange::ExtRange,ψ,N,k,tmax)
	B = symmetrize_operator(single_spin_op(b,j,N),N,k)
	return otoc(H,A,B,trange,ψ,tmax)
end

function otoc_spat(H,A,b,trange::ExtRange,ψ,N,tmax)
	res = zeros(length(trange),N)
	@sync for j in 1:N
		Threads.@spawn res[:,j]=calc_otoc(H,A,b,j,trange,ψ,N,tmax)
	end
	return res
end

function otoc_spat(H,A,b,trange::ExtRange,ψ,N,k,tmax)
	res = zeros(length(trange),N)
	@sync for j in 1:N
		Threads.@spawn res[:,j]=calc_otoc(H,A,b,j,trange,ψ,N,k,tmax)
	end
	return res
end

function otoc_spat!(res,H,A,b,trange::ExtRange,ψ,N,tmax)
	@sync for j in 1:N
		Threads.@spawn res[:,j]=calc_otoc(H,A,b,j,trange,ψ,N,tmax)
	end
	return res
end

function otoc_spat!(res,H,A,b,trange::ExtRange,ψ,N,k,tmax)
	@sync for j in 1:N
		Threads.@spawn res[:,j]=calc_otoc(H,A,b,j,trange,ψ,N,k,tmax)
	end
	return res
end



#############################
### Exact Diagonalisation ###
#############################

#Core functions
function Fψ(A::Matrix{ComplexF64},B::Matrix{ComplexF64},λs::Vector{Float64},t::Float64,ψ::Vector{ComplexF64}) #A,B already in eigenbasis
	eigmt = exp(-im*Diagonal(λs)*t)
	return real(ψ'*eigmt'*A*eigmt*B*eigmt'*A*eigmt*B*ψ) #much faster than trace!
end

function Ftr(A::Matrix{ComplexF64},B::Matrix{ComplexF64},λs::Vector{Float64},t::Float64) #A,B already in eigenbasis
	eigmt = exp(-im*Diagonal(λs)*t)
	return real(tr(eigmt'*A*eigmt*B*eigmt'*A*eigmt*B))
end

#Single time - Typicality
otoc_ed(A::Matrix{ComplexF64},B::Matrix{ComplexF64},λs::Vector{Float64},t::Float64,ψs::Vector{Vector{ComplexF64}},s::Int64) = mean(Fψ(A,B,λs,t,ψs[i]) for i in 1:s)

#Single time - Single Vector
otoc_edψ(A::Matrix{ComplexF64},B::Matrix{ComplexF64},λs::Vector{Float64},t::Float64,ψ::Vector{ComplexF64}) = Fψ(A,B,λs,t,ψ)

#Single time - Trace
otoc_edtr(A::Matrix{ComplexF64},B::Matrix{ComplexF64},λs::Vector{Float64},t::Float64) = Ftr(A,B,λs,t)

#Time Range - Typicality
function otoc_ed(A::Matrix{ComplexF64},B::AbstractArray{ComplexF64},λs::Vector{Float64},Q::Matrix{Float64},trange::ExtRange,ψs::Vector{Vector{ComplexF64}},s::Int64)
	res = zeros(length(trange))
	QdBQ = Q'B*Q
	for (ti,t) in enumerate(trange)
		res[ti] = otoc_ed(A,QdBQ,λs,t,ψs,s)
	end
	return res
end


#Time Range - Single Vector
function otoc_edψ(A::Matrix{ComplexF64},B::AbstractArray{ComplexF64},λs::Vector{Float64},Q::Matrix{Float64},trange::ExtRange,ψ::Vector{ComplexF64})
	res = zeros(length(trange))
	QdBQ = Q'B*Q
	for (ti,t) in enumerate(trange)
		res[ti] = otoc_edψ(A,QdBQ,λs,t,ψ)
	end
	return res
end


#Time Range - Trace
function otoc_edtr(A::Matrix{ComplexF64},B::AbstractArray{ComplexF64},λs::Vector{Float64},Q::Matrix{Float64},trange::ExtRange)
	res = zeros(length(trange))
	QdBQ = Q'B*Q
	for (ti,t) in enumerate(trange)
		res[ti] = otoc_edtr(A,QdBQ,λs,t)
	end
	return res
end


###########################################################
### ED OTOCs computed over multiple spatial indices ###
###########################################################

#Any Times - Typicality
function otoc_spat_ed(A::Matrix{ComplexF64},b::AbstractArray{ComplexF64},λs::Vector{Float64},Q::Matrix{Float64},ts::TvExtRange,N::Int64,s::Int64) #b=σ(xyz) in original basis
	res = zeros(length(ts),N)
	ψs = Vector{Vector{ComplexF64}}(undef,s)
	for ind in 1:s
		ψs[ind] = random_state(N)
	end
	@sync for j in 1:N
		Threads.@spawn res[:,j] .= otoc_ed(A,single_spin_op(b,j,N),λs,Q,ts,ψs,s)
	end
	return res
end
function otoc_spat_ed(A::Matrix{ComplexF64},b::AbstractArray{ComplexF64},λs::Vector{Float64},Q::Matrix{Float64},ts::TvExtRange,N::Int64,s::Int64,symsec::Int64,dim::Int64) #b=σ(xyz) in original basis
	res = zeros(length(ts),N)
	ψs = Vector{Vector{ComplexF64}}(undef,s)
	for ind in 1:s
		ψs[ind] = random_state(N,dim)
	end
	@sync for j in 1:N
		Threads.@spawn res[:,j] .= otoc_ed(A,symmetrize_operator(single_spin_op(b,j,N),N,symsec),λs,Q,ts,ψs,s)
	end
	return res
end

#Any Times - Single Vector
function otoc_spat_edψ(A::Matrix{ComplexF64},b::AbstractArray{ComplexF64},λs::Vector{Float64},Q::Matrix{Float64},ts::TvExtRange,ψ::Vector{ComplexF64},N::Int64) #b=σ(xyz) in original basis
	res = zeros(length(ts),N)
	@sync for j in 1:N
		Threads.@spawn res[:,j] .= otoc_edψ(A,single_spin_op(b,j,N),λs,Q,ts,ψ)
	end
	return res
end
function otoc_spat_edψ(A::Matrix{ComplexF64},b::AbstractArray{ComplexF64},λs::Vector{Float64},Q::Matrix{Float64},ts::TvExtRange,ψ::Vector{ComplexF64},N::Int64,symsec::Int64) #b=σ(xyz) in original basis
	res = zeros(length(ts),N)
	@sync for j in 1:N
		Threads.@spawn res[:,j] .= otoc_edψ(A,symmetrize_operator(single_spin_op(b,j,N),N,symsec),λs,Q,ts,ψ) #ψ already projected down
	end
	return res
end

#Any Times - Trace
function otoc_spat_edtr(A::Matrix{ComplexF64},b::AbstractArray{ComplexF64},λs::Vector{Float64},Q::Matrix{Float64},ts::TvExtRange,N::Int64) #b=σ(xyz) in original basis
	res = zeros(length(ts),N)
	@sync for j in 1:N
		Threads.@spawn res[:,j] .= otoc_edtr(A,single_spin_op(b,j,N),λs,Q,ts)
	end
	return res
end
function otoc_spat_edtr(A::Matrix{ComplexF64},b::AbstractArray{ComplexF64},λs::Vector{Float64},Q::Matrix{Float64},ts::TvExtRange,N::Int64,symsec::Int64) #b=σ(xyz) in original basis
	res = zeros(length(ts),N)
	@sync for j in 1:N
		Threads.@spawn res[:,j] .= otoc_edtr(A,symmetrize_operator(single_spin_op(b,j,N),N,symsec),λs,Q,ts)
	end
	return res
end

#####################################
### Diagonlize and Calculate OTOC ###
#####################################

#Any Times - Typicality
function Diag_OTOC(H::Matrix{ComplexF64},A::SparseMatrixCSC{ComplexF64,Int64},b::SparseMatrixCSC{ComplexF64,Int64},ts::TvExtRange,N::Int64,s::Int64)
	λs, Q = eigen!(H)
	Q = convert(Matrix{Float64},Q)
	logmsg("Diagonalized H.")
	QdAQ =  Q'*A*Q
	return otoc_spat_ed(QdAQ,b,λs,Q,ts,N,s)
end
function Diag_OTOC(H::Matrix{ComplexF64},A::SparseMatrixCSC{ComplexF64,Int64},b::SparseMatrixCSC{ComplexF64,Int64},ts::TvExtRange,N::Int64,s::Int64,symsec::Int64,dim::Int64)
	λs, Q = eigen!(H)
	Q = convert(Matrix{Float64},Q)
	logmsg("Diagonalized H.")
	QdAQ =  Q'*A*Q
	return otoc_spat_ed(QdAQ,b,λs,Q,ts,N,s,symsec,dim)
end

function DiagOTOC_multiplestates(H,A,b,ts,N,s)
	λs, Q = eigen!(H)
	Q = convert(Matrix{Float64},Q)
	logmsg("Diagonalized H.")
	QdAQ =  Q'*A*Q
	oto = zeros(length(ts),N,s)
	for i in 1:s
		ψ = random_state(N)
		oto[:,:,i] = otoc_spat_edψ(QdAQ,b,λs,Q,ts,ψ,N)
	end
	return oto
end

function DiagOTOC_multiplestates(H,A,b,ts,N,s,symsec,dim)
	λs, Q = eigen!(H)
	Q = convert(Matrix{Float64},Q)
	logmsg("Diagonalized H.")
	QdAQ =  Q'*A*Q
	oto = zeros(length(ts),N,s)
	for i in 1:s
		ψ = random_state(N,dim)
		oto[:,:,i] = otoc_spat_edψ(QdAQ,b,λs,Q,ts,ψ,N,symsec)
	end
	return oto
end

#Any Times - Single Vector
function Diag_OTOCψ(H::Matrix{ComplexF64},A::SparseMatrixCSC{ComplexF64,Int64},b::SparseMatrixCSC{ComplexF64,Int64},ts::TvExtRange,ψ::Vector{ComplexF64},N::Int64)
	λs, Q = eigen!(H)
	Q = convert(Matrix{Float64},Q)
	logmsg("Diagonalized H.")
	QdAQ =  Q'*A*Q
	Qdψ = Q'*ψ
	return otoc_spat_edψ(QdAQ,b,λs,Q,ts,Qdψ,N)
end
function Diag_OTOCψ(H::Matrix{ComplexF64},A::SparseMatrixCSC{ComplexF64,Int64},b::SparseMatrixCSC{ComplexF64,Int64},ts::TvExtRange,ψ::Vector{ComplexF64},N::Int64,symsec::Int64)
	λs, Q = eigen!(H)
	Q = convert(Matrix{Float64},Q)
	logmsg("Diagonalized H.")
	QdAQ =  Q'*A*Q
	Qdψ = Q'*ψ
	return otoc_spat_edψ(QdAQ,b,λs,Q,ts,Qdψ,N,symsec)
end

#Any Times - Trace
function Diag_OTOCtr(H::Matrix{ComplexF64},A::SparseMatrixCSC{ComplexF64,Int64},b::SparseMatrixCSC{ComplexF64,Int64},ts::TvExtRange,N::Int64)
	λs, Q = eigen!(H)
	Q = convert(Matrix{Float64},Q)
	logmsg("Diagonalized H.")
	QdAQ =  Q'*A*Q
	return otoc_spat_edtr(QdAQ,b,λs,Q,ts,N)
end
function Diag_OTOCtr(H::Matrix{ComplexF64},A::SparseMatrixCSC{ComplexF64,Int64},b::SparseMatrixCSC{ComplexF64,Int64},ts::TvExtRange,N::Int64,symsec::Int64)
	λs, Q = eigen!(H)
	Q = convert(Matrix{Float64},Q)
	logmsg("Diagonalized H.")
	QdAQ =  Q'*A*Q
	return otoc_spat_edtr(QdAQ,b,λs,Q,ts,N,symsec)
end

##########################################
### Old implementations for comparison ###
##########################################

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