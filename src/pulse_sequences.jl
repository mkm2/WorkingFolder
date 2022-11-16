module PulseSequences
using SparseArrays, LinearAlgebra, Plots
using SpinSymmetry, Statistics
using KrylovKit
using ..LightCones

export Pulse, SlowPulse, FastPulse, Sequence
export duration, is_valid, get_sequence
export hamiltonian, rotation, nσ, global_nσ
export WAHUHA, WAHUHA_ZF, WAHUHA_FR, Rhim71, Rhim71_ZF, Rhim71_FR


##########################
### Implement Pulses   ###
##########################

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

##########################
### Implement Sequence ###
##########################

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


########################################
### Pulse Rotations and Hamiltonians ###
########################################

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


#######################
### Named Sequences ###
#######################


function WAHUHA(τ1::Float64,τ2::Float64,τ3::Float64) #Simple WAHUHA
    pulses = Vector{Pulse}(undef,4)
    pulses[1]  = FastPulse("+x",π/2)
    pulses[2]  = FastPulse("-y",π/2)
    pulses[3]  = FastPulse("+y",π/2)
    pulses[4]  = FastPulse("-x",π/2)
    return Sequence(pulses,[τ1,τ2,2*τ3,τ2,τ1])
end
WAHUHA(τ1::Float64) = WAHUHA(τ1,2*τ1,2*τ1)

function WAHUHA_ZF(τ1::Float64,τ2::Float64,τ3::Float64) #WAHUHA removing fields
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

function WAHUHA_FR(τ1::Float64,τ2::Float64,τ3::Float64,τ4::Float64) #WAHUHA reversing fields
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

function Rhim71(Ω::Float64)  #Simple magic echo
    pulses = Vector{Pulse}(undef,4)
    pulses[1] = FastPulse("-y",π/2.0)
    pulses[2] = SlowPulse("x",π,Ω)
    pulses[3] = SlowPulse("-x",π,Ω)
    pulses[4] = FastPulse("+y",π/2.0)
    return Sequence(pulses,[0.,0.,0.,0.,0.])
end

function Rhim71_ZF(Ω::Float64) #Magic echo removing fields
    pulses = Vector{Pulse}(undef,5)
    pulses[1] = FastPulse("-y",π/2.0)
    pulses[2] = SlowPulse("x",π,Ω)
    pulses[3] = FastPulse("+x",π/1.0)
    pulses[4] = SlowPulse("-x",π,Ω)
    pulses[5] = FastPulse("+y",π/2.0)
    #pulses[6] = FastPulse("+z",π/1.0)
    return Sequence(pulses,[0.,0.,0.,0.,0.,0.])
end

function Rhim71_FR(Ω::Float64) #Magic echo reversing fields
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
    
    return Sequence(pulses,[0.,0.,0.,0.,0.,π/Ω,0.,0.,0.,0.,0.])
end


function get_sequence(sequence_name::String,t::Float64,n::Int)
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
    return seq
end

end #module