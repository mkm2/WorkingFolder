module Basics

using SparseArrays, LinearAlgebra, Plots
using SpinSymmetry, XXZNumerics
using KrylovKit
using Reexport
using ..LightCones

@reexport using ..GeomPos

export σplus, σminus, σz, σx, ⊗, Δ, 𝟙
export chainJ, correlator, single_spin_op,xxz
export field_term, random_state, random_product_state
export magnetisation

const σplus = sparse([1],[2],[1],2,2)
const σminus = sparse([2],[1],[1],2,2)
const σz = spdiagm([1,-1])
const σx = sparse([1,2],[2,1],[1,1])
const σy = sparse([1,2],[2,1],[-im,+im])
const ⊗ = kron

const Δ = -2

speye(k) = spdiagm(ones(k))
𝟙(N) = speye(2^N)

function chainJ(N,α=6)
    Symmetric(diagm([k=>k^-float(α) * ones(N-k) for k in 1:N-1]...))
end

correlator(op1, i,j,N) = correlator(op1,op1,i,j,N)
function correlator(op1, op2,i,j,N)
    i > j && return correlator(op2,op1,j,i,N)
    return 𝟙(i-1) ⊗ op1 ⊗ 𝟙(j-i-1) ⊗ op2 ⊗ 𝟙(N-j)
end
single_spin_op(op, k::Integer, N::Integer) = 𝟙(k-1) ⊗ op ⊗ 𝟙(N-k)


xxz(J::AbstractMatrix) = xxz(J, size(J,1))
xxz(N::Int, α=6) = xxz(chainJ(N,α), N)
function xxz(J::AbstractMatrix, N)
    res = spzeros(Float64, 2^N, 2^N)
    for i in 1:size(J,1)
        for j in i+1:size(J,2)
            res += J[i,j]*(2*correlator(σplus,σminus,i,j,N) + 2*correlator(σminus,σplus,i,j,N) + Δ*correlator(σz, i,j,N))
        end
    end
    return res
end

function field_term(h::Float64, N::Int)
    res = spzeros(Float64, 2^N, 2^N)
    hs = -h*ones(N) + 2*h*rand(Float64,N) #uniform distribution in [-h,+h]
    for i in 1:N
        res += hs[i]*single_spin_op(σz,i,N)
    end
    return res
end

hamiltonian_from_positions(pd::PositionData,shot::Int) = hamiltonian_from_positions(pd.data[:,shot],geometry_from_density(pd.descriptor.geometry,pd.descriptor.density,pd.descriptor.system_size,1))

function hamiltonian_from_positions(positions::Vector{Float64},geometry::Union{Box, BoxPBC, NoisyChain, NoisyChainPBC};α=6)
    interaction = PowerLaw(α)
    return xxz(interaction_matrix(interaction,distance_matrix(geometry,positions)))
end

function random_state(N::Int)
	return normalize!(randn(ComplexF64,2^N))
end

function random_product_state(N::Int)
	gen = (random_state(1) for i in 1:N)
	return kron(gen...)
end

function magnetisation(σ,ψ,N)
	S = 0
	for i in 1:N
		σ_i = single_spin_op(σ,i,N)
		S += dot(ψ,σ_i,ψ)
	end
	return S
end

end #module