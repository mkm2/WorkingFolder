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
using SparseArrays, LinearAlgebra, Plots, KrylovKit, SpinSymmetry, BenchmarkTools, PlutoUI

# ╔═╡ 2eb4c90c-83dd-426e-be43-5058fcee8178
html"""<style>main {max-width: 60%;}</style>"""

# ╔═╡ d21dd1f5-3d5f-4561-afbf-0777c38fe011
TableOfContents()

# ╔═╡ e850d5f9-b62a-4df7-8dd9-d7edaadf2b80
Ω0 = 100.

# ╔═╡ 570cb354-ad5f-4234-9470-1968c13fad26
begin
	abstract type Pulse end
	
	struct SlowPulse <: Pulse
		n::Vector{Float64} #Rotation vector
		α::Real #Rotation angle = Ωt
		Ω::Float64 #Rabi frequency Ω
	end
	function SlowPulse(axis::String,α::Real,Ω::Float64)
		if length(axis) == 1
			ax = axis[1]
			sgn = '+'
		else
			ax = axis[2]
			sgn = axis[1]
		end
		if sgn == '+'
			sign = 1.0
		elseif sgn == '-'
			sign = -1.0
		else
			throw("Wrong format. Use +x or -x.")
		end
		if ax == 'x'
			return SlowPulse(Vector{Float64}([sign*1,0,0]),α,Ω)
		elseif ax == 'y'
			return SlowPulse(Vector{Float64}([0,sign*1,0]),α,Ω)
		elseif ax == 'z'
			return SlowPulse(Vector{Float64}([0,0,sign*1]),α,Ω)
		else
			throw("Incorrect axis given.")
		end
	end

	struct FastPulse <: Pulse
		n::Vector{Float64} #Rotation vector
		α::Real #Rotation angle
	end
	function FastPulse(axis::String,α::Real)
		if length(axis) == 1
			ax = axis[1]
			sgn = '+'
		else
			ax = axis[2]
			sgn = axis[1]
		end
		if sgn == '+'
			sign = 1.0
		elseif sgn == '-'
			sign = -1.0
		else
			throw("Wrong format. Use +x or -x.")
		end
		if ax == 'x'
			return FastPulse(Vector{Float64}([sign*1,0,0]),α)
		elseif ax == 'y'
			return FastPulse(Vector{Float64}([0,sign*1,0]),α)
		elseif ax == 'z'
			return FastPulse(Vector{Float64}([0,0,sign*1]),α)
		else
			throw("Incorrect axis given.")
		end
	end

	duration(pulse::FastPulse) = 0.0
	duration(pulse::SlowPulse) = pulse.α/pulse.Ω
	
	struct Sequence
		pulses::Vector{Pulse}
		pulse_times::Vector{Float64}
		n_pulses::Int
		n_fast::Int
		n_slow::Int
		τs::Vector{Float64}
		n_τs::Int
		tc::Float64
	end
	Sequence(pulses::Vector{Pulse},τs::Vector{Float64}) = Sequence(pulses,duration.(pulses), length(pulses), count(p -> p isa FastPulse, pulses), count(p -> p isa SlowPulse, pulses), τs, length(τs), sum(cat(duration.(pulses),τs;dims=1)))
	Sequence(pulses::Vector{Pulse}) = Sequence(pulses,zeros(1+length(pulses)))
	is_valid(seq::Sequence) = (seq.n_τs == 0 || seq.n_τs == seq.n_pulses + 1) ? true : false
end

# ╔═╡ 43cd2361-c747-4309-afdc-abbbb5bc0894
begin
	function global_nσ(nx::Float64,ny::Float64,nz::Float64, N::Int)
		op = nσ(nx,ny,nz)
		if N == 1
			return op
		else
			op_full = spzeros(2^N,2^N)
			for i in 1:N
				op_full += single_spin_op(op,i,N)
			end
			return op_full
		end
	end

	function nσ(nx::Float64,ny::Float64,nz::Float64)
		return nx*σx + ny*σy + nz*σz
	end

	function rotation(pulse::FastPulse)
		return cos(pulse.α/2.0) * 𝟙(1) - im*sin(pulse.α/2.0) * nσ(pulse.n[1],pulse.n[2],pulse.n[3])
	end

	function rotation(pulse::FastPulse,N::Int)
		op = rotation(pulse)
		if N == 1
			return op
		else
			return kron((op for i in 1:N)...)
		end
	end
	
	function hamiltonian(pulse::SlowPulse,N::Int)
		return pulse.Ω/2.0 * global_nσ(pulse.n[1],pulse.n[2],pulse.n[3],N)
	end

	function rotation(op::SparseMatrixCSC{ComplexF64},ϕ::Real,N::Int)
		return cos(ϕ/2.0) * 𝟙(N) - im*sin(ϕ/2.0) * op
	end
end

# ╔═╡ 944ec3d3-7076-4732-8637-fb94c7bc687d
round.(rotation(FastPulse("x",π/2),2)'* correlator(σx,σx,1,2,2) * rotation(FastPulse("x",π/2),2),digits=3)

# ╔═╡ 3ec6e41c-7106-4303-acff-45ceaff5bb46
single_spin_op(σx,2,2)

# ╔═╡ 7439aa08-1bc5-4d40-9483-66adab3f980c
correlator(σx,σx,1,2,2)

# ╔═╡ 56a74c99-df7f-4c61-8a57-801e232b2e27
round.(rotation(FastPulse("-y",π/2),2)'* single_spin_op(σz,1,2) * rotation(FastPulse("-y",π/2),2),digits=3)

# ╔═╡ 90f191c3-a1b5-48ca-ae8f-d272595009d0
σx*σz

# ╔═╡ 74003e3a-af07-4702-8093-06f57b390505
rotation(FastPulse("x",π/2),2)*[1,0,0,0]

# ╔═╡ be389577-7444-42be-8fe6-1a9fcfc2df42
begin
	function WAHUHA(τ1::Float64,τ2::Float64,τ3::Float64)
		pulses = Vector{Pulse}(undef,4)
		pulses[1]  = FastPulse("+x",π/2)
		pulses[2]  = FastPulse("-y",π/2)
		pulses[3]  = FastPulse("+y",π/2)
		pulses[4]  = FastPulse("-x",π/2)
		return Sequence(pulses,[τ1,τ2,2*τ3,τ2,τ1])
	end
	WAHUHA(τ1::Float64) = WAHUHA(τ1,2*τ1,2*τ1)

	function WAHUHA_ZF(τ1::Float64,τ2::Float64,τ3::Float64)
		pulses = Vector{Pulse}(undef,6)
		pulses[1]  = FastPulse("+x",π/2)
		pulses[2]  = FastPulse("-y",π/2)
		pulses[3] = FastPulse("+y",π/1)
		pulses[4]  = FastPulse("+y",π/2)
		pulses[5]  = FastPulse("x",π/2)
		pulses[6]  = FastPulse("-x",π/1)
		return Sequence(pulses,[τ1,τ2,τ3,τ3,τ2,τ1,0.])
	end
	WAHUHA_ZF(τ1::Float64) = WAHUHA_ZF(τ1,2*τ1,2*τ1)

	function WAHUHA_FR(τ1::Float64,τ2::Float64,τ3::Float64,τ4::Float64)
		pulses = Vector{Pulse}(undef,10)
		pulses[1] = FastPulse("+x",π/2)
		pulses[2] = FastPulse("-y",π/2)
		pulses[3] = FastPulse("y",π/1)
		pulses[4] = FastPulse("y",π/2)
		pulses[5] = FastPulse("x",π/2)
		pulses[6] = FastPulse("-x",π/2)
		pulses[7] = FastPulse("-y",π/2)
		pulses[8] = FastPulse("-y",π/1)
		pulses[9] = FastPulse("y",π/2)
		pulses[10]  = FastPulse("-x",π/2)
		return Sequence(pulses,[τ1,τ2,τ3,τ3,τ2,τ1+τ4+τ1,τ2,τ3,τ3,τ2,τ1])
	end
	WAHUHA_FR(τ1::Float64) = WAHUHA_FR(τ1,2*τ1,2*τ1,2*τ1)
	
	function Rhim71(Ω::Float64)
		pulses = Vector{Pulse}(undef,4)
		pulses[1] = FastPulse("-y",π/2.0)
		pulses[2] = SlowPulse("x",π,Ω)
		pulses[3] = SlowPulse("-x",π,Ω)
		pulses[4] = FastPulse("+y",π/2.0)
		seq = Sequence(pulses,[0.,0.,0.,0.,0.])
	end

	function Rhim71_ZF(Ω::Float64)
		pulses = Vector{Pulse}(undef,5)
		pulses[1] = FastPulse("-y",π/2.0)
		pulses[2] = SlowPulse("x",π,Ω)
		pulses[3] = FastPulse("+x",π/1.0)
		pulses[4] = SlowPulse("-x",π,Ω)
		pulses[5] = FastPulse("+y",π/2.0)
		#pulses[6] = FastPulse("+z",π/1.0)
		seq = Sequence(pulses,[0.,0.,0.,0.,0.,0.])
	end

	function Rhim71_FR(Ω::Float64)
		pulses = Vector{Pulse}(undef,10)
		pulses[1] = FastPulse("-y",π/2.0)
		pulses[2] = SlowPulse("+x",π,Ω)
		pulses[3] = FastPulse("+x",π/1.0)
		pulses[4] = SlowPulse("-x",π,Ω)
		pulses[5] = FastPulse("-y",π/2.0)

		pulses[6] = FastPulse("+y",π/2.0)
		pulses[7] = SlowPulse("+x",π,Ω)
		pulses[8] = FastPulse("-x",π/1.0)
		pulses[9] = SlowPulse("-x",π,Ω)
		pulses[10] = FastPulse("+y",π/2.0)
		
		seq = Sequence(pulses,[0.,0.,0.,0.,0.,π/Ω,0.,0.,0.,0.,0.])
	end
	
end

# ╔═╡ 85886233-a2ee-49e8-8549-4cb7ad464f10
WAHUHAFR(1.0)

# ╔═╡ bd42a527-e3b0-42b0-a941-dccfe67123c2
begin
	function evolve_forward(Hint::SparseMatrixCSC{ComplexF64},t::Float64,ψ::Vector{ComplexF64},method::String)
		if method == "ED"
			λs, Q = eigen!(Matrix(Hint))
			return Q*exp(-im*Diagonal(λs)*t)*Q'*ψ
		elseif method == "Krylov"
			return exponentiate(Hint,-im*t,ψ;ishermitian=true)[1]
		else
			throw("Method $(method) not supported.")
		end
	end
	
	function floquet_drive(Hint::SparseMatrixCSC{ComplexF64},t::Float64,ψ::Vector{ComplexF64},N::Int,seq::Sequence,n::Int,method::String)
		### ED ###
		if method == "ED"
			#Set up operations
			if any(τ->τ>0,seq.τs)
				λs_f, Q_f = eigen!(Matrix(Hint)) 
			end
			rotations = Vector{SparseMatrixCSC{ComplexF64}}([spzeros(2^N,2^N) for k in 1:seq.n_fast])
			λs = Vector{Vector{Float64}}([zeros(2^N) for k in 1:seq.n_slow])
			Qs = Vector{Matrix{ComplexF64}}([zeros(2^N,2^N) for k in 1:seq.n_slow])
			k_fast = 0
			k_slow = 0
			for pulse in seq.pulses
				if pulse isa FastPulse
					k_fast += 1
					rotations[k_fast] = rotation(pulse,N)
				elseif pulse isa SlowPulse
					k_slow += 1
					λs[k_slow], Qs[k_slow] = eigen!(Matrix(Hint+hamiltonian(pulse,N)))
				end
			end
			
			#Apply pulses
			for iter in 1:n
				k_fast = 0
				k_slow = 0
				for (k,pulse) in enumerate(seq.pulses)
					if seq.τs[k] > 0
						ψ = Q_f*exp(-im*Diagonal(λs_f)*seq.τs[k])*Q_f'*ψ
					end
					if pulse isa FastPulse
						k_fast += 1
						ψ = rotations[k_fast] * ψ
					elseif pulse isa SlowPulse
						k_slow += 1
						ψ = Qs[k_slow]*exp(-im*Diagonal(λs[k_slow])*seq.pulse_times[k])*Qs[k_slow]'*ψ
					end
				end
				if seq.τs[seq.n_τs] > 0
					ψ = Q_f*exp(-im*Diagonal(λs_f)*seq.τs[seq.n_τs])*Q_f'*ψ
				end
			end
			return ψ
			

		### KRYLOV ###
		elseif method == "Krylov"
			throw("Not yet.")
		else
			throw("Method $(method) not supported.")
		end
	end
end

# ╔═╡ ae6b7339-01f4-4eb4-852a-6f836d240a52
function echo(Hint::SparseMatrixCSC{ComplexF64},t::Float64,ψ::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String)
	if sequence_name == "WAHUHA"
		seq = WAHUHA(t/(2*n)) #5t = n*tc; tc = 10*τ1 => τ1 = t/(2*n)
	elseif sequence_name == "WAHUHA_ZF"
		seq = WAHUHA_ZF(t/(2*n)) #5t = n*tc; tc = 10*τ1 => τ1 = t/(2*n)
	elseif sequence_name == "WAHUHA_FR"
			seq = WAHUHA_FR(t/(2*n)) #11t = n*tc; tc = 22*τ1 => τ1 = t/(2*n)
	elseif sequence_name == "Rhim71"
		seq = Rhim71(π*n/t) #2t = n*tc; tc = 2*π/Ω => Ω = π*n/t
	elseif sequence_name == "Rhim71_ZF"
		seq = Rhim71_ZF(π*n/t) #2t = n*tc; tc = 2*π/Ω => Ω = π*n/t
	elseif sequence_name == "Rhim71_FR"
		seq = Rhim71_FR(π*n/t) #5t = n*tc; tc = 5*π/Ω => Ω = π*n/t
	else
		throw("Unknown sequence.")
	end
	ψ = evolve_forward(Hint,t,ψ,method)
	ψ = floquet_drive(Hint,t,ψ,N,seq,n,method)
	return ψ
end

# ╔═╡ ec1c5839-e606-4365-b60c-befd5149914c
begin
	function echo(Hint::SparseMatrixCSC{ComplexF64},A::SparseMatrixCSC{ComplexF64},ϕ::Float64,t::Float64,ψ::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String)
		if sequence_name == "WAHUHA"
		seq = WAHUHA(t/(2*n)) #5t = n*tc; tc = 10*τ1 => τ1 = t/(2*n)
	elseif sequence_name == "WAHUHA_ZF"
		seq = WAHUHA_ZF(t/(2*n)) #5t = n*tc; tc = 10*τ1 => τ1 = t/(2*n)
	elseif sequence_name == "WAHUHA_FR"
			seq = WAHUHA_FR(t/(2*n)) #11t = n*tc; tc = 22*τ1 => τ1 = t/(2*n)
	elseif sequence_name == "Rhim71"
		seq = Rhim71(π*n/t) #2t = n*tc; tc = 2*π/Ω => Ω = π*n/t
	elseif sequence_name == "Rhim71_ZF"
		seq = Rhim71_ZF(π*n/t) #2t = n*tc; tc = 2*π/Ω => Ω = π*n/t
	elseif sequence_name == "Rhim71_FR"
		seq = Rhim71_FR(π*n/t) #5t = n*tc; tc = 5*π/Ω => Ω = π*n/t
	else
		throw("Unknown sequence.")
	end
		ψ = evolve_forward(Hint,t,ψ,method)
		ψ = exp(-im*ϕ/2*Matrix(A)) * ψ
		ψ = floquet_drive(Hint,t,ψ,N,seq,n,method)
		return ψ
	end
	echo(Hint::SparseMatrixCSC{ComplexF64},A::SparseMatrixCSC{ComplexF64},t::Float64,ψ::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String) = echo(Hint,A,π/1.,t,ψ,sequence_name,n,N,method)
end

# ╔═╡ 0b1ea4bf-5a83-4523-8438-da7e08f65e1e
begin
	function echo(Hint::SparseMatrixCSC{ComplexF64},Hf,A::SparseMatrixCSC{ComplexF64},ϕ::Float64,t::Float64,ψ::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String)
		if sequence_name == "WAHUHA"
		seq = WAHUHA(t/(2*n)) #5t = n*tc; tc = 10*τ1 => τ1 = t/(2*n)
	elseif sequence_name == "WAHUHA_ZF"
		seq = WAHUHA_ZF(t/(2*n)) #5t = n*tc; tc = 10*τ1 => τ1 = t/(2*n)
	elseif sequence_name == "WAHUHA_FR"
			seq = WAHUHA_FR(t/(2*n)) #11t = n*tc; tc = 22*τ1 => τ1 = t/(2*n)
	elseif sequence_name == "Rhim71"
		seq = Rhim71(π*n/t) #2t = n*tc; tc = 2*π/Ω => Ω = π*n/t
	elseif sequence_name == "Rhim71_ZF"
		seq = Rhim71_ZF(π*n/t) #2t = n*tc; tc = 2*π/Ω => Ω = π*n/t
	elseif sequence_name == "Rhim71_FR"
		seq = Rhim71_FR(π*n/t) #5t = n*tc; tc = 5*π/Ω => Ω = π*n/t
	else
		throw("Unknown sequence.")
	end
		ψ = evolve_forward(Hint+Hf,t,ψ,method)
		ψ = exp(-im*ϕ/2*Matrix(A)) * ψ
		ψ = floquet_drive(Hint+Hf,t,ψ,N,seq,n,method)
		return ψ
	end
	echo(Hint::SparseMatrixCSC{ComplexF64},Hf,A::SparseMatrixCSC{ComplexF64},t::Float64,ψ::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String) = echo(Hint,Hf,A,π/1.,t,ψ,sequence_name,n,N,method)
end

# ╔═╡ 7e0de451-bb15-49da-84a0-38d575da658b
begin
	N = 9
	i = div(N,2)+1
	H = convert(SparseMatrixCSC{ComplexF64},xxz(N,6))
	A = convert(SparseMatrixCSC{ComplexF64},single_spin_op(σx,i,N))
	B = convert(SparseMatrixCSC{ComplexF64},σx)
	ψ0 = [0,1.0+0.0*im,0,0]
	δt = 0.1
	trange = 0:δt:5
end

# ╔═╡ 4aa08654-873f-4443-9930-e86c373579a2
md"# Fidelity Test WAHUHA"

# ╔═╡ df835e06-c963-4d57-b524-33aa8291aec0
begin
	state = kron(leftx,rightx)#,up,rightx,rightx,rightx)
end

# ╔═╡ aec39b24-c065-4152-a2d1-b281fd4359c7
function apply(state,τ)
	state_tmpr = evolve_forward(H,τ,state,"ED")
	state_tmpr = rotation(FastPulse("+x",π/2),N)*state_tmpr
	state_tmpr = evolve_forward(H,2*τ,state_tmpr,"ED")
	state_tmpr = rotation(FastPulse("-y",π/2),N)*state_tmpr
	state_tmpr = evolve_forward(H,4*τ,state_tmpr,"ED")
	state_tmpr = rotation(FastPulse("+y",π/2),N)*state_tmpr
	state_tmpr = evolve_forward(H,2*τ,state_tmpr,"ED")
	state_tmpr = rotation(FastPulse("-x",π/2),N)*state_tmpr
	state_tmpr = evolve_forward(H,τ,state_tmpr,"ED")
	return state_tmpr
end

# ╔═╡ aa97d91e-9300-4d53-94ca-b2451d6c436d
md"# Fidelity Test Rhim71"

# ╔═╡ 4cd80398-edcd-404e-9065-c3088e84ce99


# ╔═╡ 2f0a8c3e-fafa-4420-b69c-6c00066891c1
function applyRhim(state,Ω,t,N)
		state_tmpr = rotation(FastPulse("-y",π/2),N)*state
		state_tmpr = evolve_forward(H+hamiltonian(SlowPulse("x",π,Ω),N),π/Ω,state_tmpr,"ED")
		state_tmpr = evolve_forward(H+hamiltonian(SlowPulse("-x",π,Ω),N),π/Ω,state_tmpr,"ED")
		state_tmpr = rotation(FastPulse("+y",π/2),N)*state_tmpr
	return state_tmpr
end

# ╔═╡ 4f027f64-e663-4a7a-a7ff-c4f497dc8dbc
H

# ╔═╡ 99b66eac-dbda-492f-b863-59748d1fd1f2
md" # Test OTOC"

# ╔═╡ effd5287-4d0d-41a1-9983-d4a229d68f8e
f=field_term(12.,N)

# ╔═╡ 37ab30ee-0206-4060-89c3-2fefc8fe1dd9
function measure_all(B::SparseMatrixCSC{ComplexF64},ψ::Vector,N::Int)
	res = zeros(ComplexF64,N)
	for i in 1:N
		res[i] = ψ'single_spin_op(B,i,N)*ψ
	end
	return res
end

# ╔═╡ 093ea9de-4163-4712-8f16-ec83fb1146b3
state_o = random_state(N)#kron(leftx,leftx,leftx,leftx,leftx,leftx,leftx,leftx,leftx)

# ╔═╡ efe0885b-41d5-4252-857e-1191785d5e51
begin
	ts = 0:1:8*π/1.
	res = zeros(length(ts))
	res2 = zeros(length(ts))
	n = 2000
	for (i,t) in enumerate(ts)
		res[i] = abs(state_o'evolve_forward(H,-t,state_o,"ED"))^2#real(magnetisation(σx,evolve_forward(H,-t,state,"ED"),N))#norm(ψ1'echo(H,A,t,ψ1,"WAHUHA",1,N,"ED"))^2
	
		τ = t/(2*n)
		seq = WAHUHA(τ)
		@debug t
		state_tmpr = echo(H+field_term(6.0,N),t,state_o,"WAHUHA_FR",n,N,"ED")
		res2[i] = abs(state_o'state_tmpr)^2# real(magnetisation(σx,state_tmpr,N))
	end
end

# ╔═╡ bdb7ac4a-136e-40a6-90ca-a8da36f5bdb6
begin
	plot(ts,res)
	plot!(ts,res2,label="Floquet")
end

# ╔═╡ b2213312-64c6-4d66-bf27-0fbf75e6867c
begin
	tsRhim = 0.1:1:8*π/1.
	resRhim2 = zeros(length(tsRhim))
	resRhim = zeros(length(tsRhim))
	nRhim = 2000
	for (i,t) in enumerate(tsRhim)
		resRhim[i] = abs(state_o'evolve_forward(H,-t,state_o,"ED"))^2#real(magnetisation(σx,evolve_forward(H,-t,state,"ED"),N))#norm(ψ1'echo(H,A,t,ψ1,"WAHUHA",1,N,"ED"))^2

		@debug t
		seq = Rhim71(π*nRhim/t)
		state_tmpr = echo(H+field_term(6.0,N),t,state_o,"Rhim71_FR",nRhim,N,"ED")
		resRhim2[i] = abs(state_o'state_tmpr)^2# real(magnetisation(σx,state_tmpr,N))
	end
end

# ╔═╡ 55f1c138-9d73-407e-922a-58cc8e7191f1
begin
	plot(tsRhim,resRhim)
	plot!(tsRhim,resRhim2,label="Floquet")
end

# ╔═╡ 24f376db-77a3-46ae-a738-f58dac76b57b
A2 = convert(SparseMatrixCSC{ComplexF64},single_spin_op(σz,i,N))

# ╔═╡ d3ec883e-cca0-4cdf-b644-0916d67d0ae8
begin
	tso = 0.1:0.5:5*π/1.
	reso = zeros(length(tso))
	reso2 = zeros(length(tso))
	resoto = zeros(ComplexF64,length(tso),N)
	no = 2000
	for (i,t) in enumerate(tso)
		reso[i] = abs(state_o'evolve_forward(H,-t,state_o,"ED"))^2#real(magnetisation(σx,evolve_forward(H,-t,state,"ED"),N))#norm(ψ1'echo(H,A,t,ψ1,"WAHUHA",1,N,"ED"))^2
		@debug t
		state_tmpr = echo(H,f,A2,π/1.,t,state_o,"Rhim71_FR",no,N,"ED")
		reso2[i] = abs(state_o'state_tmpr)^2# real(magnetisation(σx,state_tmpr,N))
		resoto[i,:] = 2*(ones(N)+real(measure_all(B,state_tmpr,N)))
	end
end

# ╔═╡ 12f52048-18fb-4852-887b-cab28af71f6b
resoto[6,:]

# ╔═╡ 59bf7abb-12a1-4649-a3e2-e77ea5b260c8
begin
	plot(tso,reso,ylims=[0,1])
	plot!(tso,reso2,label="Floquet")
end

# ╔═╡ 368f34e9-0e42-43b3-a27e-bb56857d88ec
begin
	λs, Q = eigen!(Matrix(H+f))
	Q = convert(Matrix{Float64},Q)
	QdAQ =  Q'*A2*Q
	oto_corr = 2*(ones(length(tso),N)-otoc_spat_edψ(QdAQ,B,λs,Q,tso,Q'state_o,N))
end

# ╔═╡ 6dd2eb77-9e38-48cb-a145-caa64839898b
begin
	k = 5
	plot(tso,real.(resoto[:,k]),label="Floquet",xlim=[0,25])
	plot!(tso,oto_corr[:,k],label="Exact")
end

# ╔═╡ 71c7c94c-50ef-4347-bec0-cbc18088f9c2
plot(tso,oto_corr[:,4],xlim=(0,5))

# ╔═╡ 401ce41e-7e35-40db-ad06-7a5a31c691aa


# ╔═╡ Cell order:
# ╠═e2b065a0-1fb2-11ed-300e-d79a41c0557a
# ╠═65b9ceeb-38b3-496a-a10a-58d322f7a5bb
# ╠═02725372-99f4-4af9-a9e8-d81444b5af98
# ╠═2eb4c90c-83dd-426e-be43-5058fcee8178
# ╠═d21dd1f5-3d5f-4561-afbf-0777c38fe011
# ╠═e850d5f9-b62a-4df7-8dd9-d7edaadf2b80
# ╠═570cb354-ad5f-4234-9470-1968c13fad26
# ╠═43cd2361-c747-4309-afdc-abbbb5bc0894
# ╠═944ec3d3-7076-4732-8637-fb94c7bc687d
# ╠═3ec6e41c-7106-4303-acff-45ceaff5bb46
# ╠═7439aa08-1bc5-4d40-9483-66adab3f980c
# ╠═56a74c99-df7f-4c61-8a57-801e232b2e27
# ╠═90f191c3-a1b5-48ca-ae8f-d272595009d0
# ╠═74003e3a-af07-4702-8093-06f57b390505
# ╠═be389577-7444-42be-8fe6-1a9fcfc2df42
# ╠═85886233-a2ee-49e8-8549-4cb7ad464f10
# ╠═bd42a527-e3b0-42b0-a941-dccfe67123c2
# ╠═ae6b7339-01f4-4eb4-852a-6f836d240a52
# ╠═ec1c5839-e606-4365-b60c-befd5149914c
# ╠═0b1ea4bf-5a83-4523-8438-da7e08f65e1e
# ╠═7e0de451-bb15-49da-84a0-38d575da658b
# ╠═4aa08654-873f-4443-9930-e86c373579a2
# ╠═df835e06-c963-4d57-b524-33aa8291aec0
# ╠═aec39b24-c065-4152-a2d1-b281fd4359c7
# ╠═efe0885b-41d5-4252-857e-1191785d5e51
# ╠═bdb7ac4a-136e-40a6-90ca-a8da36f5bdb6
# ╠═aa97d91e-9300-4d53-94ca-b2451d6c436d
# ╠═4cd80398-edcd-404e-9065-c3088e84ce99
# ╠═2f0a8c3e-fafa-4420-b69c-6c00066891c1
# ╠═4f027f64-e663-4a7a-a7ff-c4f497dc8dbc
# ╠═b2213312-64c6-4d66-bf27-0fbf75e6867c
# ╠═55f1c138-9d73-407e-922a-58cc8e7191f1
# ╠═99b66eac-dbda-492f-b863-59748d1fd1f2
# ╠═effd5287-4d0d-41a1-9983-d4a229d68f8e
# ╠═37ab30ee-0206-4060-89c3-2fefc8fe1dd9
# ╠═093ea9de-4163-4712-8f16-ec83fb1146b3
# ╠═24f376db-77a3-46ae-a738-f58dac76b57b
# ╠═d3ec883e-cca0-4cdf-b644-0916d67d0ae8
# ╠═12f52048-18fb-4852-887b-cab28af71f6b
# ╠═59bf7abb-12a1-4649-a3e2-e77ea5b260c8
# ╠═6dd2eb77-9e38-48cb-a145-caa64839898b
# ╠═368f34e9-0e42-43b3-a27e-bb56857d88ec
# ╠═71c7c94c-50ef-4347-bec0-cbc18088f9c2
# ╠═401ce41e-7e35-40db-ad06-7a5a31c691aa
