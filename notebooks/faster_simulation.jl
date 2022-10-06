### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 7b2269d0-3eaa-11ed-177e-ef37a8ec1b99
import Pkg

# ╔═╡ 3f215687-ed44-4eff-930f-aedd80a16f67
begin
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ de019eca-125a-4bd2-8ecb-9722d0337ef1
using SparseArrays, LinearAlgebra, Plots, KrylovKit, SpinSymmetry, BenchmarkTools

# ╔═╡ 82328eb1-9305-4b3a-8c34-4dab9c06e15b
html"""<style>main {max-width: 60%;}</style>"""

# ╔═╡ 4f9f2a5b-ba47-42e7-917a-f57b32f051d2
begin
	N = 8
	H = convert(SparseMatrixCSC{ComplexF64},xxz(N,6))
	ψ0 = random_state(N)
	δt = 0.1
	trange = 0:δt:5
end

# ╔═╡ 2ff95cd7-8e1c-4581-82a9-e323f4d6e967


# ╔═╡ Cell order:
# ╠═7b2269d0-3eaa-11ed-177e-ef37a8ec1b99
# ╠═3f215687-ed44-4eff-930f-aedd80a16f67
# ╠═de019eca-125a-4bd2-8ecb-9722d0337ef1
# ╠═82328eb1-9305-4b3a-8c34-4dab9c06e15b
# ╠═4f9f2a5b-ba47-42e7-917a-f57b32f051d2
# ╠═2ff95cd7-8e1c-4581-82a9-e323f4d6e967
