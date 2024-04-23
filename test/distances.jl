metric_list = [:ci, :rf]
metric_kwargs = Dict()

tree_file = "../data/hashrftrees/567-taxa-10-trees.tre"
trees = Newick.read(tree_file, Vector{RoundaboutTree}, TaxonSet())
ntree = length(trees)

for metric in metric_list
    D = distance(trees, metric)

    for j in 1:ntree, i in 1:j
        @test D[i, j] ≈ distance(trees[i], trees[j], metric)

        if i == j
            @test D[i, j] == 0.0
        else
            @test D[i, j] ≥ 0.0
        end
    end
end

@testset "Clustering Information Distance diffonly on/off" begin
    Don = distance(trees, :ci; diffonly=true)
    Doff = distance(trees, :ci; diffonly=false)

    for j in 1:ntree, i in 1:j
        if Don[i, j] == 0.0
            @test isapprox(Doff[i, j], 0.0, atol=1e-12)
        else
            @test Don[i, j] ≈ Doff[i, j]
            @test Don[i, j] ≈ distance(trees[i], trees[j], :ci; diffonly=false)
        end
    end
end

