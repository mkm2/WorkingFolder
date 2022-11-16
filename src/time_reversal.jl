module TimeReversal
using SparseArrays, LinearAlgebra, Plots
using SpinSymmetry, Statistics
using KrylovKit
using Reexport
using ..LightCones

@reexport using ..PulseSequences

export evolve_forward, perturb, floquet_drive, echo

ExtRange = Union{AbstractRange{Float64},Vector{Float64}}
TvExtRange = Union{Float64,ExtRange}

######################
### Evolve Forward ###
######################

#Single Time

function evolve_forward(H::SparseMatrixCSC{Float64},t::Float64,ψ0::Vector{ComplexF64},method::String,tmax::Float64=1.0)
    if method == "ED"
        λs, Q = eigen!(Matrix(H))
        return Q*exp(-im*Diagonal(λs)*t)*Q'*ψ0
    elseif method == "Krylov"
        return krylov_from0_alternative(H,-t,ψ0,tmax)
    else
        throw("Method $(method) not supported.")
    end
end

#Time Range - ED

function evolve_forward(λs::Vector{Float64},Q::Matrix{Float64},trange::ExtRange,ψ0::Vector{ComplexF64},N::Int)
    res = zeros(ComplexF64,2^N,length(trange))
    for (i,t) in enumerate(trange)
        res[:,i] = Q*exp(-im*Diagonal(λs)*t)*Q'*ψ0 
    end
    return res
end

#Time Range - Krylov
function evolve_forward(H::SparseMatrixCSC{Float64},trange::ExtRange,ψ0::Vector{ComplexF64},N::Int)
    res = zeros(ComplexF64,2^N,length(trange))
    res[:,1] = krylov_step(H,-trange[1],ψ0)
    for i in 2:length(trange)
        δt = trange[i] - trange[i-1]
        res[:,i] = krylov_step(H,-δt,res[:,i-1])
    end
    return res
end

###############
### Perturb ###
###############

#Single Time

function perturb(A::SparseMatrixCSC{ComplexF64},ϕ::Float64,ψ::Vector{ComplexF64})
    return exp(-im*ϕ/2*Matrix(A)) * ψ
end

#Time Range

function perturb(A::SparseMatrixCSC{ComplexF64},ϕ::Float64,ψs::Matrix{ComplexF64})
    res = zeros(ComplexF64,size(ψs))
    for i in 1:size(ψs)[2]
        res[:,i] = exp(-im*ϕ/2*Matrix(A)) * ψs[:,i]
    end
    return res
end

#####################
### Floquet Drive ###
#####################

#Single Time

function floquet_drive(H::SparseMatrixCSC{Float64},ψ::Vector{ComplexF64},N::Int,seq::Sequence,n::Int,method::String,tmax::Float64=1.0)
    ### ED ###
    if method == "ED"
        #Set up operations
        if any(τ->τ>0,seq.τs)
            λs_f, Q_f = eigen!(Matrix(H))
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
                λs[k_slow], Qs[k_slow] = eigen!(Matrix(H+hamiltonian(pulse,N)))
            end
        end
        
        #Apply pulses
        for _iter in 1:n
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
        #Set up operations
        rotations = Vector{SparseMatrixCSC{ComplexF64}}([spzeros(2^N,2^N) for k in 1:seq.n_fast])
        k_fast = 0
        for pulse in seq.pulses
            if pulse isa FastPulse
                k_fast += 1
                rotations[k_fast] = rotation(pulse,N)
            end
        end

        #Apply pulses
        for _iter in 1:n
            k_fast = 0
            for (k,pulse) in enumerate(seq.pulses)
                if seq.τs[k] > 0
                    ψ = krylov_from0_alternative(H,-seq.τs[k],ψ,tmax)
                end
                if pulse isa FastPulse
                    k_fast += 1
                    ψ = rotations[k_fast] * ψ
                elseif pulse isa SlowPulse
                    ψ = krylov_from0_alternative(H+hamiltonian(pulse,N),-seq.pulse_times[k],ψ,tmax)
                end
            end
            if seq.τs[seq.n_τs] > 0
                ψ = krylov_from0_alternative(H,-seq.τs[seq.n_τs],ψ,tmax)
            end
        end
        return ψ
    else
        throw("Method $(method) not supported.")
    end
end

#Single Time - Free ED done

function floquet_drive(H::SparseMatrixCSC{Float64},λs_f::Vector{Float64},Q_f::Matrix{Float64},ψ::Vector{ComplexF64},N::Int,seq::Sequence,n::Int)
    ### ED ###
    #Set up operations
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
            λs[k_slow], Qs[k_slow] = eigen!(Matrix(H+hamiltonian(pulse,N)))
        end
    end
    
    #Apply pulses
    for _iter in 1:n
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
end


##############################
### Echo - No Perturbation ###
##############################

#Single Time

function echo(H::SparseMatrixCSC{Float64},t::Float64,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String)
	seq = get_sequence(sequence_name,t,n)
	ψ = evolve_forward(H,t,ψ0,method)
	ψ = floquet_drive(H,ψ,N,seq,n,method)
	return ψ
end

#Time Range

function echo(H::SparseMatrixCSC{Float64},trange::ExtRange,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String,tmax::Float64=1.0)
    ψs = zeros(2^N,length(trange))
    if method == "ED"
        λs, Q = eigen!(Matrix(H))
        Q = convert(Matrix{Float64},Q)
        logmsg("Diagonalized H.")
        ψs = evolve_forward(λs,Q,trange,ψ0,N)
        for (i,t) in enumerate(trange)
            seq = get_sequence(sequence_name,t,n)
            ψs[:,i] = floquet_drive(H,λs,Q,ψs[:,i],N,seq,n)
        end
    elseif method == "Krylov"
        ψs = evolve_forward(H,trange,ψ0,N)
        for (i,t) in enumerate(trange)
            seq = get_sequence(sequence_name,t,n)
            ψs[:,i] = floquet_drive(H,ψs[:,i],N,seq,n,method,tmax)
        end
    end
    return ψs
end

#############################################
### Echo - Perturbation by exp(-im*ϕ/2*A) ###
#############################################

#Single Time

function echo(H::SparseMatrixCSC{Float64},A::SparseMatrixCSC{ComplexF64},ϕ::Float64,t::Float64,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String)
    seq = get_sequence(sequence_name,t,n)
    ψ = evolve_forward(H,t,ψ0,method)
    ψ = perturb(A,ϕ,ψ)
    ψ = floquet_drive(H,ψ,N,seq,n,method)
    return ψ
end
echo(H::SparseMatrixCSC{Float64},A::SparseMatrixCSC{ComplexF64},t::Float64,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String) = echo(H,A,π/1.,t,ψ0,sequence_name,n,N,method) #ϕ=π


#Time Range

function echo(H::SparseMatrixCSC{Float64},A::SparseMatrixCSC{ComplexF64},ϕ::Float64,trange::ExtRange,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String,tmax::Float64=1.0)
    ψs = zeros(2^N,length(trange))
    if method == "ED"
        λs, Q = eigen!(Matrix(H))
        Q = convert(Matrix{Float64},Q)
        logmsg("Diagonalized H.")
        ψs = evolve_forward(λs,Q,trange,ψ0,N)
        ψs = perturb(A,ϕ,ψs)
        for (i,t) in enumerate(trange)
            seq = get_sequence(sequence_name,t,n)
            ψs[:,i] = floquet_drive(H,λs,Q,ψs[:,i],N,seq,n)
        end
    elseif method == "Krylov"
        ψs = evolve_forward(H,trange,ψ0,N)
        ψs = perturb(A,ϕ,ψs)
        for (i,t) in enumerate(trange)
            seq = get_sequence(sequence_name,t,n)
            ψs[:,i] = floquet_drive(H,ψs[:,i],N,seq,n,method,tmax)
        end
    end
    return ψs
end
echo(H::SparseMatrixCSC{Float64},A::SparseMatrixCSC{ComplexF64},trange::ExtRange,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String,tmax::Float64=1.0) = echo(H,A,π/1.,trange,ψ0,sequence_name,n,N,method,tmax) #ϕ=π


end #module