### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 12015c66-5b99-43b6-bc9e-2a0b9ea41ae7
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 923436d4-70ec-11ed-1d80-592332204a67
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 1bce0382-b112-4bb8-83f4-b6f8127ab823
TableOfContents()

# ╔═╡ 27d3d294-f772-482a-8e6a-ce9d8f523466
begin
	geom = :box_pbc
	N = 13
	SHOTS = 100000
	DISORDER_PARAM= 0.6
	DISORDER_PARAMS = [1.99,1.8,1.5,1.2,0.9,0.6]
	s = 2/DISORDER_PARAM
	σ = 1.5*(s-1)
end

# ╔═╡ faf69760-1e0f-43cf-890b-e1b4ab2e7e5f
md"# Play"

# ╔═╡ c6b64851-cf7c-4e2a-be9d-c856e9e4d7aa
s

# ╔═╡ 45829f06-50a9-4dff-b499-d12da0ce2c77
σ

# ╔═╡ eb1a46d3-a207-4805-a3f0-6b59040312ad
begin
	desc = PositionDataDescriptor(geom,N,SHOTS,DISORDER_PARAM)
	pd = create(desc)
end

# ╔═╡ d41bf3e3-6efd-456c-84eb-9404631391cb
histogram(pd.data[1,:],bins=100)

# ╔═╡ 0d343411-1b72-4692-ade3-c7f8572ece7b
histogram(pd.data[2,:],bins=100)

# ╔═╡ 40e7da8a-3217-4e61-a87c-9de2a498202b
histogram(pd.data[3,:],bins=100)

# ╔═╡ 7f9a04f4-39a0-4625-ac05-7c2414b71a94
histogram(pd.distances[1,2,:],bins=100)

# ╔═╡ 4458344c-fcf4-497c-ad36-7d2889fe095c
histogram(pd.distances[1,3,:],bins=100)

# ╔═╡ 7c91bc4d-cdea-42de-9f32-91174d8b7c7b
md"# Generate for plots"

# ╔═╡ 37d62929-33c4-474e-852e-2e4d7117c1c1
path = pwd()

# ╔═╡ 361a1972-0b76-4a11-b779-74ff567b2539
pd.descriptor

# ╔═╡ aa7981f6-85aa-44de-bbfe-a68a0ce69456
function save_positions(posdata,datapath)
	jldopen(datapath, "w") do file
		file["descriptor"] = posdata.descriptor
		file["positions"] = posdata.data
		file["distances"] = posdata.distances
	end
	print("Saved positions.\n")
end

# ╔═╡ 80dcbc04-d8e8-4e03-bb5d-394ab7061855
begin
	geom1 = :noisy_chain_pbc
	geom2 = :box_pbc
end

# ╔═╡ 7071fa55-7b46-4ef5-99ff-a6009a692960
for ρ in DISORDER_PARAMS
	print(ρ)
	desc1 = PositionDataDescriptor(geom1,N,SHOTS,ρ)
	pd1 = create(desc1)
	save_positions(pd1,path*"/positions_ncpbc_$(ρ).jld2")
	print(desc1)
end

# ╔═╡ b47beffd-dfbb-498c-a9b7-c584c932ef56
reverse(DISORDER_PARAMS)

# ╔═╡ 5ae4c351-21ec-45e6-89ca-818345f821d1
for ρ in reverse(DISORDER_PARAMS)
	print(ρ)
	desc2 = PositionDataDescriptor(geom2,N,SHOTS,ρ)
	pd2 = create(desc2)
	save_positions(pd2,path*"/positions_boxpbc_$(ρ).jld2")
end

# ╔═╡ Cell order:
# ╠═923436d4-70ec-11ed-1d80-592332204a67
# ╠═12015c66-5b99-43b6-bc9e-2a0b9ea41ae7
# ╠═1bce0382-b112-4bb8-83f4-b6f8127ab823
# ╠═27d3d294-f772-482a-8e6a-ce9d8f523466
# ╠═faf69760-1e0f-43cf-890b-e1b4ab2e7e5f
# ╠═c6b64851-cf7c-4e2a-be9d-c856e9e4d7aa
# ╠═45829f06-50a9-4dff-b499-d12da0ce2c77
# ╠═eb1a46d3-a207-4805-a3f0-6b59040312ad
# ╠═d41bf3e3-6efd-456c-84eb-9404631391cb
# ╠═0d343411-1b72-4692-ade3-c7f8572ece7b
# ╠═40e7da8a-3217-4e61-a87c-9de2a498202b
# ╠═7f9a04f4-39a0-4625-ac05-7c2414b71a94
# ╠═4458344c-fcf4-497c-ad36-7d2889fe095c
# ╠═7c91bc4d-cdea-42de-9f32-91174d8b7c7b
# ╠═37d62929-33c4-474e-852e-2e4d7117c1c1
# ╠═361a1972-0b76-4a11-b779-74ff567b2539
# ╠═aa7981f6-85aa-44de-bbfe-a68a0ce69456
# ╠═80dcbc04-d8e8-4e03-bb5d-394ab7061855
# ╠═7071fa55-7b46-4ef5-99ff-a6009a692960
# ╠═b47beffd-dfbb-498c-a9b7-c584c932ef56
# ╠═5ae4c351-21ec-45e6-89ca-818345f821d1
