### A Pluto.jl notebook ###
# v0.19.17

using Markdown
using InteractiveUtils

# ╔═╡ d097b376-ba6c-403d-bc81-26f4ea5134c7
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 99afc594-0961-435a-8792-6e738d6106a8
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 85117564-27be-4e3d-adae-0a11aa893ec4
TableOfContents()

# ╔═╡ fef0fdcb-aa8c-4d25-9371-97a447e5c694
hs = [12]

# ╔═╡ 1eeaa9ff-8b31-4f10-ab40-052e46a6bf71
path = pwd()

# ╔═╡ 268d61b7-f119-41c6-8ac7-d2f1c04ffd73
div(17,2)+1

# ╔═╡ c5d6759e-efce-4fbc-806f-92f03d992295
folders = ["/h = 12/"]

# ╔═╡ a44aa14d-b972-47d4-9738-5c7b4a6083c0
begin
	hds = []
	for s in readdir(path*folders[1])
		if !occursin("otoc",s) && !occursin("combined",s) && !occursin("77",s)
			tmp = load(path*folders[1]*s,"hs")
			for i in 1:25
					append!(hds,[tmp[i][9]])
			end
		end
	end
end

# ╔═╡ 79345997-8d5a-4c42-aa0a-b4c6934269f5
minimum(hds)

# ╔═╡ 51db9505-5a08-4658-9ae9-610dea242884
histogram(hds)

# ╔═╡ ec037253-9b2e-489c-bb42-6186fbcabbe1
load(path*folders[1]*"7851832_N17_RS.jld2","hs")

# ╔═╡ 5cc93d17-e70d-48ae-867b-955f7d293fd3


# ╔═╡ Cell order:
# ╠═99afc594-0961-435a-8792-6e738d6106a8
# ╠═d097b376-ba6c-403d-bc81-26f4ea5134c7
# ╠═85117564-27be-4e3d-adae-0a11aa893ec4
# ╠═fef0fdcb-aa8c-4d25-9371-97a447e5c694
# ╠═1eeaa9ff-8b31-4f10-ab40-052e46a6bf71
# ╠═268d61b7-f119-41c6-8ac7-d2f1c04ffd73
# ╠═c5d6759e-efce-4fbc-806f-92f03d992295
# ╠═a44aa14d-b972-47d4-9738-5c7b4a6083c0
# ╠═79345997-8d5a-4c42-aa0a-b4c6934269f5
# ╠═51db9505-5a08-4658-9ae9-610dea242884
# ╠═ec037253-9b2e-489c-bb42-6186fbcabbe1
# ╠═5cc93d17-e70d-48ae-867b-955f7d293fd3
