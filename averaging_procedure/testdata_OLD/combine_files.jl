### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ b7108bea-80ee-4bd0-87ba-4be0201c15ac
begin
	import Pkg
	Pkg.add(url="https://github.com/mkm2/LightCones.git")
	using LightCones
end

# ╔═╡ 5983335e-4d64-11ed-2bec-0f03ae137366
using LinearAlgebra,Plots,JLD2,Statistics,PlutoUI

# ╔═╡ 1163b6a2-fbac-44d3-82e1-787f5d946b25
function combine_files(files,path,new_file)
	jobids = Vector{String}(undef,length(files))
	params = Vector{Any}(undef,length(files))
	data = Vector{Array{Float64,4}}(undef,length(files))
	#positiondata = Vector{PositionData}(undef, length(files))
	for (i,f) in enumerate(files)
		jobids[i] = load(path*f,"jobid")
		params[i] = load(path*f,"params")
		#positiondata[i] = load(path*f,"positiondata")
		data[i] = load(path*f,"data")
	end
	jldopen(path*new_file, "w") do file
        file["data"] = data
		#file["positiondata"] = positiondata
        file["params"] = params
        file["jobid"] = jobids
    end
end

# ╔═╡ 00955629-b15c-495c-8ab2-dbf515092afc
path = pwd()*"/"

# ╔═╡ dcdcf323-262b-485d-9932-187f9803acc3
md"# data ref"

# ╔═╡ dd7aa8e3-720d-405a-9401-7e95ee3daa79
begin
	l_ref0 = cat(["7296575_N13_RS.jld2"],["73474$(i)_N13_RS.jld2" for i in 46:54],dims=1)
	l_ref4 = cat(["7296576_N13_RS.jld2"],["73474$(i)_N13_RS.jld2" for i in 55:63],dims=1)
	l_ref12 = cat(["7296577_N13_RS.jld2"],["73474$(i)_N13_RS.jld2" for i in 64:68],dims=1)
	append!(l_ref12,["73474$(i)_N13_RS.jld2" for i in 70:73])
	combine_files(l_ref0,path,"combined_file_N13_ref0.jld2")
	combine_files(l_ref4,path,"combined_file_N13_ref4.jld2")
	combine_files(l_ref12,path,"combined_file_N13_ref12.jld2")
end

# ╔═╡ a43770ec-f2a8-4b04-acf2-773136832fa6
md"# data h"

# ╔═╡ 84a61afb-1043-4f82-9304-0383482b143d
begin
	l_h0 = cat(["7301823_N13_BS.jld2"],["73474$(i)_N13_BS.jld2" for i in 74:77],dims=1)
	l_h4 = cat(["7301824_N13_BS.jld2"],["73474$(i)_N13_BS.jld2" for i in 78:81],dims=1)
	l_h12 = cat(["7301825_N13_BS.jld2"],["73474$(i)_N13_BS.jld2" for i in 82:85],dims=1)
	combine_files(l_h0,path,"combined_file_N13_h0.jld2")
	combine_files(l_h4,path,"combined_file_N13_h4.jld2")
	combine_files(l_h12,path,"combined_file_N13_h12.jld2")
end

# ╔═╡ 4c5ce4a8-8101-45de-9e56-f144cefbb58e


# ╔═╡ c208ac9e-5978-4bbd-8698-73ae59e2de7d
md"# data psi"

# ╔═╡ ffb55a80-91ed-4cdc-b3e8-457ee78022aa
begin
	l_psi0 = cat(["7273738_N13_BS.jld2"],["73474$(i)_N13_BS.jld2" for i in 86:89],dims=1)
	l_psi4 = cat(["7273739_N13_BS.jld2"],["73474$(i)_N13_BS.jld2" for i in 90:93],dims=1)
	#l_psi12 = cat(["7273740_N13_BS.jld2"],["73474$(i)_N13_BS.jld2" for i in 64:68],dims=1)
	#append!(l_psi12,["73474$(i)_N13_RS.jld2" for i in 70:73])
	combine_files(l_psi0,path,"combined_file_N13_psi0.jld2")
	combine_files(l_psi4,path,"combined_file_N13_psi4.jld2")
	#combine_files(l_psi12,path,"combined_file_N13_psi12.jld2")
end

# ╔═╡ 75a46ea4-637d-443a-8098-87c66fdb607f


# ╔═╡ Cell order:
# ╠═5983335e-4d64-11ed-2bec-0f03ae137366
# ╠═b7108bea-80ee-4bd0-87ba-4be0201c15ac
# ╠═1163b6a2-fbac-44d3-82e1-787f5d946b25
# ╠═00955629-b15c-495c-8ab2-dbf515092afc
# ╠═dcdcf323-262b-485d-9932-187f9803acc3
# ╠═dd7aa8e3-720d-405a-9401-7e95ee3daa79
# ╠═a43770ec-f2a8-4b04-acf2-773136832fa6
# ╠═84a61afb-1043-4f82-9304-0383482b143d
# ╠═4c5ce4a8-8101-45de-9e56-f144cefbb58e
# ╠═c208ac9e-5978-4bbd-8698-73ae59e2de7d
# ╠═ffb55a80-91ed-4cdc-b3e8-457ee78022aa
# ╠═75a46ea4-637d-443a-8098-87c66fdb607f
