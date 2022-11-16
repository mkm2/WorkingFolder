module Basics

using SparseArrays, LinearAlgebra, Plots
using SpinSymmetry
using KrylovKit
using Reexport
using ..LightCones

@reexport using ..Geom
@reexport using ..Interactions
@reexport using ..GeomPos

export σplus, σminus, σz, σx, σy, ⊗, Δ, 𝟙
export up, down, rightx, leftx, neel_state, neel_x_state
export chainJ, chainJ_pbc, correlator, single_spin_op
export xxz, xxz_pbc, xyz, xyz_pbc, field_term, hamiltonian_from_positions, const_field
export nearest_neighbourJ, nearest_neighbourJ_pbc
export random_state, random_product_state, random_bit_state, random_bitstring_state, random_basis_state
export magnetisation, fidelity, measure_at_j, measure_all, sign_of_eigenstate, signs_of_eigenstate, otoc_by_eigenstate_measurement


const σplus = sparse([1],[2],[1.0],2,2)
const σminus = sparse([2],[1],[1.0],2,2)
const σz = spdiagm([1.0,-1.0])
const σx = sparse([1,2],[2,1],[1.0,1.0])
const σy = sparse([1,2],[2,1],[-im,+im])
const ⊗ = kron

const Δ = -2

const up = convert(Vector{ComplexF64},[1.0, 0.0])
const down = convert(Vector{ComplexF64},[0.0, 1.0])
const rightx = convert(Vector{ComplexF64},[1.0, 1.0]/sqrt(2.0))
const leftx = convert(Vector{ComplexF64},[1.0, -1.0]/sqrt(2.0))

speye(k) = spdiagm(ones(k))
𝟙(N) = speye(2^N)

function chainJ(N,α=6)
    Symmetric(diagm([k=>k^-float(α) * ones(N-k) for k in 1:N-1]...))
end

chainJ_pbc(N,α=6) = interaction_matrix(PowerLaw(α),distance_matrix(RegularChainPBC(N,1.0),[k*1.0 for k in 1:N]))

function nearest_neighbourJ(N)
    Symmetric(diagm(1=>ones(N-1)))
end

function nearest_neighbourJ_pbc(N)
    Symmetric(diagm(1=>ones(N-1), N-1=>ones(1)))
end

correlator(op1, i,j,N) = correlator(op1,op1,i,j,N)
function correlator(op1, op2,i,j,N)
    i > j && return correlator(op2,op1,j,i,N)
    return 𝟙(i-1) ⊗ op1 ⊗ 𝟙(j-i-1) ⊗ op2 ⊗ 𝟙(N-j)
end
single_spin_op(op, k::Integer, N::Integer) = 𝟙(k-1) ⊗ op ⊗ 𝟙(N-k)


xxz(J::AbstractMatrix) = xxz(J, size(J,1))
xxz(N::Int, α=6) = xxz(chainJ(N,α), N)
xxz_pbc(N::Int,α=6) = xxz(chainJ_pbc(N,α), N)
function xxz(J::AbstractMatrix, N)
    res = spzeros(Float64, 2^N, 2^N)
    for i in 1:size(J,1)
        for j in i+1:size(J,2)
            res += J[i,j]*(2*correlator(σplus,σminus,i,j,N) + 2*correlator(σminus,σplus,i,j,N) + Δ*correlator(σz, i,j,N)) #σxσx+σyσy-2σzσz
        end
    end
    return res
end

xyz(J::AbstractMatrix,Jx,Jy,Jz) = xyz(J,Jx,Jy,Jz,size(J,1))
xyz(N::Int,Jx,Jy,Jz,α=6) = xyz(chainJ(N,α),Jx,Jy,Jz,N)
xyz_pbc(N::Int,Jx,Jy,Jz,α=6) = xyz(chainJ_pbc(N,α),Jx,Jy,Jz,N)
function xyz(J::AbstractMatrix, Jx, Jy, Jz, N) #xxz = xyz with Jx = Jy = 1, Jz = Δ = -2
    res = spzeros(Float64, 2^N, 2^N)
    for i in 1:size(J,1)
        for j in i+1:size(J,2)
            res += J[i,j]*(Jx*correlator(σx,σx,i,j,N) + Jy*correlator(σy,σy,i,j,N) + Jz*correlator(σz, i,j,N))
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

function field_term(h::Float64, N::Int, k::Int)
    res = spzeros(Float64, 2^N, 2^N)
    hs = -h*ones(N) + 2*h*rand(Float64,N) #uniform distribution in [-h,+h]
    for i in 1:N
        res += hs[i]*single_spin_op(σz,i,N)
    end
    return symmetrize_operator(res,N,k)
end

function const_field(h::Float64, N::Int)
    res = spzeros(Float64, 2^N, 2^N)
    for i in 1:N
        res += h*single_spin_op(σz,i,N)
    end
    return res
end

hamiltonian_from_positions(pd::PositionData,shot::Int) = hamiltonian_from_positions(pd.data[:,shot],geometry_from_density(pd.descriptor.geometry,pd.descriptor.ρ,pd.descriptor.system_size,1))

function hamiltonian_from_positions(positions::Vector{Float64},geometry::Union{Box, BoxPBC, NoisyChain, NoisyChainPBC};α=6)
    interaction = PowerLaw(α)
    return xxz(interaction_matrix(interaction,distance_matrix(geometry,positions)))
end

function random_state(N::Int, d::Int) #N particles, restricted to d dimensions
    return normalize!(randn(ComplexF64,d))
end
random_state(N::Int) = random_state(N, 2^N)

function random_product_state(N::Int)
	gen = (random_state(1) for i in 1:N)
	return kron(gen...)
end
function random_product_state(N::Int, k::Int) #sector of k spins |↑⟩
    x = random_product_state(N)
    return normalize!(symmetrize_state(x,N,k))
end

function random_bitstring_state(N::Int, d::Int) #simply choose random basis state
    x = zeros(ComplexF64,d)
    ind = rand(1:d)
    x[ind] = 1
    return x
end
random_bitstring_state(N::Int) = random_bitstring_state(N, 2^N)

function random_basis_state(axis::String)
    if axis == "z"
        return rand() < 0.5 ? up : down
    elseif axis == "x"
        return rand() < 0.5 ? leftx : rightx
    else
        throw("Axis not implemented.")
    end
end

function random_bitstring_state(N::Int, axis::String) #simply choose random basis state
    gen = (random_basis_state(axis) for i in 1:N)
    return kron(gen...)
end

function neel_state(N::Int)
    N%2 == 0 || throw("N has to be even to create a Néel state.")
    return kron((i%2 == 1 ? up : down for i in 1:N)...)
end

function neel_x_state(N::Int)
    N%2 == 0 || throw("N has to be even to create a Néel state.")
    return kron((i%2 == 1 ? rightx : leftx for i in 1:N)...)
end

#function random_bit_state()
#    x = Vector{ComplexF64}(undef,2)
#    tmp = convert(ComplexF64,rand(Bool,1)[1])
#    x[1] = tmp
#    x[2] = 1-tmp
#    return x
#end
#function random_bitstring_state(N::Int) old version
#    gen = (random_bit_state() for i in 1:N)
#    return kron(gen...)
#end

function magnetisation(σ,ψ,N)
	S = 0
	for i in 1:N
		σ_i = single_spin_op(σ,i,N)
		S += dot(ψ,σ_i,ψ)
	end
	return S
end

function fidelity(ψ1::Vector,ψ2::Vector)
    return abs(ψ1'ψ2)^2
end

function fidelity(ψ1::Vector,ψ2s::Matrix)
    res = zeros(size(ψ2s)[2])
    for i in 1:size(ψ2s)[2]
        res[i] = abs(ψ1'ψ2s[:,i])^2
    end
    return res
end

function measure_at_j(B,ψ,j)
    return ψ'single_spin_op(B,j,N)*ʋ
end

function measure_all(B,ψ,N)
	res = zeros(ComplexF64,N)
	for j in 1:N
		res[j] = ψ'single_spin_op(B,j,N)*ψ
	end
	return res
end

function signs_of_eigenstate(B,ψ,N)
    signs = zeros(N)
    for i in 1:N
        signs[i] = sign_of_eigenstate(single_spin_op(B,i,N),ψ)
    end
    return signs
end

function sign_of_eigenstate(B,ψ)
    if isapprox(ψ'B*ψ,1.0)
        return +1.0
    elseif isapprox(ψ'B*ψ,-1.0) 
        return -1.0
    else
        throw("ψ is not an eigenstate of B to an eigenvalue of magnitude 1.")
    end
end

function otoc_by_eigenstate_measurement(B,ψ::Vector,signs,N)
    return 2*(ones(N)-signs .* real.(measure_all(B,ψ,N)))
end

function otoc_by_eigenstate_measurement(B,ψs::Matrix,signs,N)
    res = zeros(size(ψs)[2],N)
    for i in 1:size(ψs)[2]
        res[i,:] = otoc_by_eigenstate_measurement(B,ψs[:,i],signs,N)
    end
    return res
end

end #module