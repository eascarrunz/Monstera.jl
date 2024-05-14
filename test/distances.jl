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


function test_against_reference_implementation(tree_file, D_file, mode; atol=0.0)
    trees = Newick.read(tree_file, Vector{RoundaboutTree}, TaxonSet())
    ntree = length(trees)
    
    D = distance(trees, mode)
    DT = eltype(typeof(D))
    Dref = readdlm(D_file, ' ', DT)
    
    for j in 1:ntree, i in 1:(j-1)
        @test isapprox(D[i, j], Dref[i, j], atol=atol)
    end
end

@testset "Robinson-Foulds Distance against reference implementation" begin
    @testset "RF on Binary trees" begin
        test_against_reference_implementation(
            "../data/Dtest/random_binary_trees.nwk",
            "../data/Dtest/D_RF_random_binary.txt",
            :rf
        )
    end

    @testset "RF on Polytomous trees" begin
        test_against_reference_implementation(
            "../data/Dtest/random_polytomous_trees.nwk",
            "../data/Dtest/D_RF_random_polytomous.txt",
            :rf
        )
    end
end

@testset "Clustering Information Distance against reference implementation" begin
    @testset "MCI on Binary trees" begin
        test_against_reference_implementation(
            "../data/Dtest/random_binary_trees.nwk",
            "../data/Dtest/D_MCI_random_binary.txt",
            :ci
        )
    end

    @testset "MCI on Polytomous trees" begin
        test_against_reference_implementation(
            "../data/Dtest/random_polytomous_trees.nwk",
            "../data/Dtest/D_MCI_random_polytomous.txt",
            :ci
            )
    end
end



