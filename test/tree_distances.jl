### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 868207ec-e19a-11ee-1ec3-eb94d07803c4
begin
    import Pkg;
    # activate a temporary environment
    Pkg.activate(mktempdir());
    Pkg.add([
        Pkg.PackageSpec(name="Plots"),
        Pkg.PackageSpec(name="PlutoUI"),
        Pkg.PackageSpec(name="Monstera", url="https://github.com/eascarrunz/Monstera.jl"),
    ]);
    using Plots, PlutoUI, Monstera;
end

# ╔═╡ aaeb74ba-312f-4ea4-83cc-2477676c4d94
using Statistics

# ╔═╡ a345c526-268a-455a-be15-f070c9155ce5
import PlutoUI: combine

# ╔═╡ fa7cd028-bd59-454c-8549-b87de708947c
@bind tree_params confirm(
	combine() do bind
		md"""
		- Number of taxa (`ntax`): $(bind(Scrubbable(8:2048; default=16)))
		- Number of random NNIs to apply (`n_nni`): $(bind(Scrubbable(1:512; default=6)))
		"""
	end
)

# ╔═╡ 2f36a28b-21af-41e5-bad8-afb4a5da100c
ntax, n_nni = tree_params;

# ╔═╡ a87fecf2-8ca7-4a05-bc46-17ff96af1125
taxonset = TaxonSet(ntax);

# ╔═╡ 0d0505b4-e353-45e6-8833-ddbfc5ab2347
initial_tree = randtree(RoundaboutTree, taxonset, false)

# ╔═╡ cf46de24-1b9c-41e2-87ef-6aea74839da8
md"""
**Initial tree:**
"""

# ╔═╡ 8f85c0a7-81df-4e7e-a9b1-9b3df29aab08
textplot(initial_tree)

# ╔═╡ a8c7e337-0481-4954-8700-c001a4cc6014


# ╔═╡ 02165576-7205-499b-a1fd-d0ce4e02fbfe


# ╔═╡ 2801207b-4fcc-42b7-bf71-66ac40e2916d


# ╔═╡ e1cb348e-4de3-4fac-b8aa-6495f9c2aaa5


# ╔═╡ 0cef6173-7749-43bb-9af6-9a702fac5b50


# ╔═╡ 2473159f-0c21-407f-8632-5a2be195c5f0


# ╔═╡ c7e2b74a-0dbd-4db2-85e2-fe0deaf8888d


# ╔═╡ bac15aea-5031-4308-b8be-09ad946628dc


# ╔═╡ 8b4821bd-a535-463f-9881-d10e8156e7d3
function rand_nni!(tree, info=false)
    branch = rand(tree.branches)
    left_node, right_node = nodes_flanking(branch)
    while isouter(left_node) || isouter(right_node)
        branch = rand(tree.branches)
        left_node, right_node = nodes_flanking(branch)
    end
    left_children = children(branch => left_node) |> collect
    right_children = children(branch => right_node) |> collect
    left_target = rand(left_children)
    right_target = rand(right_children)
    info && @info "Swapping $(last(left_target))) and $(last(right_target))"
    swapnodes!(left_target, right_target)
end

# ╔═╡ 3aaf490f-409e-4583-8d49-5469009b4836
begin
	trees = [clone(initial_tree) for _ in 1:3]
	for tree in trees
		for _ in 1:n_nni
			rand_nni!(tree)
		end
	end
end

# ╔═╡ 27a8dc83-4a8d-4d0a-b8e2-a25145d47218
textplot(trees[1])

# ╔═╡ 187b96e4-7ca1-4afd-96ff-af337381a5ae
textplot(trees[2])

# ╔═╡ b7bb5a7e-1903-4633-9de8-7b0866ad63a9
textplot(trees[3])

# ╔═╡ 3b753518-d8ab-438b-89b1-ab9436530afe
d = distance(trees, :ci)

# ╔═╡ 82ceee22-3d4c-4180-83a7-dbc08537643a
md"""
Mean distance between trees: **$(Statistics.mean(d))**
"""

# ╔═╡ 061dd972-d42f-4a15-bd25-36649818fa07
length(trees)

# ╔═╡ Cell order:
# ╟─868207ec-e19a-11ee-1ec3-eb94d07803c4
# ╠═aaeb74ba-312f-4ea4-83cc-2477676c4d94
# ╠═a345c526-268a-455a-be15-f070c9155ce5
# ╟─fa7cd028-bd59-454c-8549-b87de708947c
# ╟─82ceee22-3d4c-4180-83a7-dbc08537643a
# ╠═2f36a28b-21af-41e5-bad8-afb4a5da100c
# ╠═a87fecf2-8ca7-4a05-bc46-17ff96af1125
# ╟─0d0505b4-e353-45e6-8833-ddbfc5ab2347
# ╟─cf46de24-1b9c-41e2-87ef-6aea74839da8
# ╟─8f85c0a7-81df-4e7e-a9b1-9b3df29aab08
# ╠═3aaf490f-409e-4583-8d49-5469009b4836
# ╠═27a8dc83-4a8d-4d0a-b8e2-a25145d47218
# ╠═187b96e4-7ca1-4afd-96ff-af337381a5ae
# ╠═b7bb5a7e-1903-4633-9de8-7b0866ad63a9
# ╠═a8c7e337-0481-4954-8700-c001a4cc6014
# ╠═02165576-7205-499b-a1fd-d0ce4e02fbfe
# ╠═2801207b-4fcc-42b7-bf71-66ac40e2916d
# ╠═3b753518-d8ab-438b-89b1-ab9436530afe
# ╠═e1cb348e-4de3-4fac-b8aa-6495f9c2aaa5
# ╠═0cef6173-7749-43bb-9af6-9a702fac5b50
# ╠═2473159f-0c21-407f-8632-5a2be195c5f0
# ╠═061dd972-d42f-4a15-bd25-36649818fa07
# ╠═c7e2b74a-0dbd-4db2-85e2-fe0deaf8888d
# ╠═bac15aea-5031-4308-b8be-09ad946628dc
# ╟─8b4821bd-a535-463f-9881-d10e8156e7d3
