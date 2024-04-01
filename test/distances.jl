metric_list = [:ci, :rf]

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
