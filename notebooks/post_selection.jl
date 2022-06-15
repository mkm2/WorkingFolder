### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 7e334382-ebfe-11ec-1ab8-c1aa60d9fde2
using LinearAlgebra, Random, Statistics

# ╔═╡ f1b77e9c-a51f-462a-90fa-1e72ea64dce4
function draw_random()
		return normalize!(randn(ComplexF64,2))
end

# ╔═╡ 0a84b101-8e69-48b6-99c4-5af6c508d636
begin
	S = Int(1e6)
	pfail=0
	for i in 1:S
		a,b = draw_random()
		pfail += 0.5*norm(a-b)^2
	end
	pfail = pfail/S
end

# ╔═╡ 4a54ee11-ac4b-4d89-80c7-20fd02fb771b
a,b = draw_random()

# ╔═╡ 90a9400c-6d87-4691-945d-9e2887243647
(norm(a+b)^2+norm(a-b)^2)/2

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
"""

# ╔═╡ Cell order:
# ╠═7e334382-ebfe-11ec-1ab8-c1aa60d9fde2
# ╠═f1b77e9c-a51f-462a-90fa-1e72ea64dce4
# ╠═0a84b101-8e69-48b6-99c4-5af6c508d636
# ╠═4a54ee11-ac4b-4d89-80c7-20fd02fb771b
# ╠═90a9400c-6d87-4691-945d-9e2887243647
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
