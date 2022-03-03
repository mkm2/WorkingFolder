using SparseArrays, LinearAlgebra, Plots
using SpinSymmetry
using KrylovKit

const σplus = sparse([2],[1],[1],2,2)
const σminus = sparse([1],[2],[1],2,2)
const σz = spdiagm([1,-1])
const σx = sparse([1,2],[2,1],[1,1])
const ⊗ = kron

const Δ = 2

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

function otoc(H,A,B,t,ψ)
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

N = 8
H = xxz(N,6)
ψ0 = normalize!(ones(2^N))

function computeJ(posdata, shot, rhoIndex)
    rho = posdata.ρs[rhoIndex]
    J = interaction_matrix(
            PowerLaw(6), 
            geometry_from_density(posdata.geometry, rho, posdata.system_size,1),
            posdata.data[:,:,shot,rhoIndex])
end

function hamiltonian(J)
    H = xxz(J, -0.7)
    return symmetrize_operator(H, N, div(N-1,2))
end

op1 = single_spin_op(σz,5,N)
op2 = single_spin_op(σz,1,N)

print("tr")

let trange = 0:0.1:5,
	p = plot(; xlabel="time t", ylabel="<σ_i(t)σ_jσ_i(t)σ_j>", legend=nothing)
	plot!(trange, 2*ones(length(trange))-2*otoc.(Ref(H), Ref(op1), Ref(op2), trange, Ref(ψ0)))
    savefig(p,"test.png")
end
