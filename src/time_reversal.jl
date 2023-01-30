module TimeReversal
using SparseArrays, LinearAlgebra, Plots
using SpinSymmetry, Statistics
using KrylovKit
using Reexport
using ..LightCones

@reexport using ..PulseSequences

export evolve_forward, perturb, floquet_drive, echo, fidelity_timereversal, oto_commutator_timereversal

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

function perturb(A::SparseMatrixCSC{ComplexF64},ψ::Vector{ComplexF64}) # with ϕ=π
    return -im * A * ψ
end

#Time Range

function perturb(A::SparseMatrixCSC{ComplexF64},ψs::Matrix{ComplexF64})# with ϕ=π
    res = zeros(ComplexF64,size(ψs))
    for i in 1:size(ψs)[2]
        res[:,i] = -im * A * ψs[:,i]
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




#Single Time Krylov - Optimized

function floquet_drive(H::SparseMatrixCSC{Float64},ψ::Vector{ComplexF64},N::Int,seq::Sequence,n::Int,rotations::Vector{Matrix{ComplexF64}},proto_hamiltonians::Vector{SparseMatrixCSC{ComplexF64}},tmax::Float64=1.0)
### KRYLOV ###
    #Apply pulses
    for _iter in 1:n
        k_fast = 0
        k_slow = 0
        for (k,pulse) in enumerate(seq.pulses)
            if seq.τs[k] > 0
                ψ = krylov_from0_alternative(H,-seq.τs[k],ψ,tmax)
            end
            if pulse isa FastPulse
                k_fast += 1
                ψ = rotations[k_fast] * ψ
            elseif pulse isa SlowPulse
                k_slow += 1
                ψ = krylov_from0_alternative(H + pulse.Ω/2.0*proto_hamiltonians[k_slow],-seq.pulse_times[k],ψ,tmax)
            end
        end
        if seq.τs[seq.n_τs] > 0
            ψ = krylov_from0_alternative(H,-seq.τs[seq.n_τs],ψ,tmax)
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
    if ϕ == π/1.
        ψ = perturb(A,ψ)
    else
        ψ = perturb(A,ϕ,ψ)
    end
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
        if ϕ == π/1.
            ψs = perturb(A,ψs)

        else
            ψs = perturb(A,ϕ,ψs)
        end
        for (i,t) in enumerate(trange)
            seq = get_sequence(sequence_name,t,n)
            ψs[:,i] = floquet_drive(H,λs,Q,ψs[:,i],N,seq,n)
        end
    elseif method == "Krylov"
        ψs = evolve_forward(H,trange,ψ0,N)
        if ϕ == π/1.
            ψs = perturb(A,ψs)

        else
            ψs = perturb(A,ϕ,ψs)
        end        
        for (i,t) in enumerate(trange)
            seq = get_sequence(sequence_name,t,n)
            ψs[:,i] = floquet_drive(H,ψs[:,i],N,seq,n,method,tmax)
        end
    end
    return ψs
end
echo(H::SparseMatrixCSC{Float64},A::SparseMatrixCSC{ComplexF64},trange::ExtRange,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String,tmax::Float64=1.0) = echo(H,A,π/1.,trange,ψ0,sequence_name,n,N,method,tmax) #ϕ=π


#More efficient than saving states.

##################
### Fidelities ###
##################

function fidelity_timereversal(H::SparseMatrixCSC{Float64},t::Float64,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String)
    seq = get_sequence(sequence_name,t,n)
    ψ = evolve_forward(H,t,ψ0,method)
    ψ = floquet_drive(H,ψ,N,seq,n,method)
    return fidelity(ψ0,ψ)
end

function fidelity_timereversal(H::SparseMatrixCSC{Float64},trange::ExtRange,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String,tmax::Float64=1.0)
    fidelities = zeros(length(trange))
    if method == "ED"
        λs, Q = eigen!(Matrix(H))
        Q = convert(Matrix{Float64},Q)
        logmsg("Diagonalized H.")
        for (i,t) in enumerate(trange)
            seq = get_sequence(sequence_name,t,n)
            ψ = evolve_forward(λs,Q,t,ψ0,N)
            ψ = floquet_drive(H,λs,Q,ψ,N,seq,n)
            fidelities[i] = fidelity(ψ0,ψ)
        end
    elseif method == "Krylov"

        seq = get_sequence(sequence_name,trange[1],n)
        ψ_forward = krylov_from0_alternative(H,-trange[1],ψ0,tmax)
        ψ = floquet_drive(H,ψ_forward,N,seq,n,method,tmax)
        fidelities[1] = fidelity(ψ0,ψ)
        for (i,t) in enumerate(trange)
            if i == 1
                continue
            end
            seq = get_sequence(sequence_name,t,n)
            δt = trange[i] - trange[i-1]
            ψ_forward = krylov_step(H,-δt,ψ_forward)
            ψ = floquet_drive(H,ψ_forward,N,seq,n,method,tmax)
            fidelities[i] = fidelity(ψ0,ψ)
        end
    end
    return fidelities
end



#Optimized for Krylov and time ranges

function fidelity_timereversal(H::SparseMatrixCSC{Float64},trange::ExtRange,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,rotations::Vector{Matrix{ComplexF64}},proto_hamiltonians::Vector{SparseMatrixCSC{ComplexF64}},tmax::Float64=1.0)
    fidelities = zeros(length(trange))
    seq = get_sequence(sequence_name,trange[1],n)
    ψ_forward = krylov_from0_alternative(H,-trange[1],ψ0,tmax)
    ψ = floquet_drive(H,ψ_forward,N,seq,n,rotations,proto_hamiltonians,tmax)
    fidelities[1] = fidelity(ψ0,ψ)
    for (i,t) in enumerate(trange)
        if i == 1
            continue
        end
        seq = get_sequence(sequence_name,t,n)
        δt = trange[i] - trange[i-1]
        ψ_forward = krylov_step(H,-δt,ψ_forward)
        ψ = floquet_drive(H,ψ_forward,N,seq,n,rotations,proto_hamiltonians,tmax)
        fidelities[i] = fidelity(ψ0,ψ)
    end
    return fidelities
end



#######################
### OTO Commutators ###
#######################


function oto_commutator_timereversal(H::SparseMatrixCSC{Float64},A::SparseMatrixCSC{ComplexF64},B::SparseMatrixCSC{ComplexF64},ϕ::Float64,t::Float64,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String)
    signs = signs_of_eigenstate(B,ψ0,N)
    seq = get_sequence(sequence_name,t,n)
    ψ = evolve_forward(H,t,ψ0,method)
    if ϕ == π/1.
        ψs = perturb(A,ψ)

    else
        ψs = perturb(A,ϕ,ψ)
    end
    ψ = floquet_drive(H,ψ,N,seq,n,method)
    return otoc_by_eigenstate_measurement(B,ψ,signs,N)

end
oto_commutator_timereversal(H::SparseMatrixCSC{Float64},A::SparseMatrixCSC{ComplexF64},B::SparseMatrixCSC{ComplexF64},t::Float64,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String) = echo(H,A,B,π/1.,t,ψ0,sequence_name,n,N,method) #ϕ=π

function oto_commutator_timereversal(H::SparseMatrixCSC{Float64},A::SparseMatrixCSC{ComplexF64},B::SparseMatrixCSC{ComplexF64},ϕ::Float64,trange::ExtRange,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String,tmax::Float64=1.0)
    oto_commutators = zeros(length(trange),N)
    signs = signs_of_eigenstate(B,ψ0,N)
    if method == "ED"
        λs, Q = eigen!(Matrix(H))
        Q = convert(Matrix{Float64},Q)
        logmsg("Diagonalized H.")
        for (i,t) in enumerate(trange)
            seq = get_sequence(sequence_name,t,n)
            ψ = evolve_forward(λs,Q,t,ψ0,N)
            if ϕ == π/1.
                ψs = perturb(A,ψ)
    
            else
                ψs = perturb(A,ϕ,ψ)
            end
            ψ = floquet_drive(H,λs,Q,ψ,N,seq,n)
            oto_commutators[i,:] = otoc_by_eigenstate_measurement(B,ψ,signs,N)
        end
    elseif method == "Krylov"
        seq = get_sequence(sequence_name,trange[1],n)
        ψ_forward = krylov_from0_alternative(H,-trange[1],ψ0,tmax)
        if ϕ == π/1.
            ψ_perturbed = perturb(A,ψ_forward)
        else
            ψ_perturbed = perturb(A,ϕ,ψ_forward)
        end
        ψ = floquet_drive(H,ψ_perturbed,N,seq,n,method,tmax)
        oto_commutators[1,:] = otoc_by_eigenstate_measurement(B,ψ,signs,N)
        for (i,t) in enumerate(trange)
            if i == 1
                continue
            end
            seq = get_sequence(sequence_name,t,n)
            δt = trange[i] - trange[i-1]
            ψ_forward = krylov_step(H,-δt,ψ_forward)
            if ϕ == π/1.
                ψ_perturbed = perturb(A,ψ_forward)
            else
                ψ_perturbed = perturb(A,ϕ,ψ_forward)
            end
            ψ = floquet_drive(H,ψ_perturbed,N,seq,n,method,tmax)
            oto_commutators[i,:] = otoc_by_eigenstate_measurement(B,ψ,signs,N)
        end
    end
    return oto_commutators
end
oto_commutator_timereversal(H::SparseMatrixCSC{Float64},A::SparseMatrixCSC{ComplexF64}, B::SparseMatrixCSC{ComplexF64}, trange::ExtRange,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,method::String,tmax::Float64=1.0) = oto_commutator_timereversal(H,A,B,π/1.,trange,ψ0,sequence_name,n,N,method,tmax) #ϕ=π




#Optimized for Krylov and time ranges
function oto_commutator_timereversal(H::SparseMatrixCSC{Float64},A::SparseMatrixCSC{ComplexF64},B::SparseMatrixCSC{ComplexF64},ϕ::Float64,trange::ExtRange,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,rotations::Vector{Matrix{ComplexF64}},proto_hamiltonians::Vector{SparseMatrixCSC{ComplexF64}},tmax::Float64=1.0)
    oto_commutators = zeros(length(trange),N)
    signs = signs_of_eigenstate(B,ψ0,N)
    seq = get_sequence(sequence_name,trange[1],n)
    ψ_forward = krylov_from0_alternative(H,-trange[1],ψ0,tmax)
    if ϕ == π/1.
        ψ_perturbed = perturb(A,ψ_forward)
    else
        ψ_perturbed = perturb(A,ϕ,ψ_forward)
    end
    ψ = floquet_drive(H,ψ_perturbed,N,seq,n,rotations,proto_hamiltonians,tmax)
    oto_commutators[1,:] = otoc_by_eigenstate_measurement(B,ψ,signs,N)
    for (i,t) in enumerate(trange)
        if i == 1
            continue
        end
        seq = get_sequence(sequence_name,t,n)
        δt = trange[i] - trange[i-1]
        ψ_forward = krylov_step(H,-δt,ψ_forward)
        if ϕ == π/1.
            ψ_perturbed = perturb(A,ψ_forward)
        else
            ψ_perturbed = perturb(A,ϕ,ψ_forward)
        end
        ψ = floquet_drive(H,ψ_perturbed,N,seq,n,rotations,proto_hamiltonians,tmax)
        oto_commutators[i,:] = otoc_by_eigenstate_measurement(B,ψ,signs,N)
    end
    return oto_commutators
end
oto_commutator_timereversal(H::SparseMatrixCSC{Float64},A::SparseMatrixCSC{ComplexF64}, B::SparseMatrixCSC{ComplexF64}, trange::ExtRange,ψ0::Vector{ComplexF64},sequence_name::String,n::Int,N::Int,rotations::Vector{Matrix{ComplexF64}},proto_hamiltonians::Vector{SparseMatrixCSC{ComplexF64}},tmax::Float64=1.0) = oto_commutator_timereversal(H,A,B,π/1.,trange,ψ0,sequence_name,n,N,rotations,proto_hamiltonians,tmax) #ϕ=π


end #module