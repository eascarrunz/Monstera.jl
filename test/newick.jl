for Ttype in CONCRETE_TREE_TYPES
    @testset "$Ttype" begin
        @testset "Wikipedia sample Newick trees" begin
            for (i, line) in enumerate(readlines("../data/Newick_wikipedia.nwk"))
                if i == 5
                    # The fifth line contains a tree with a branch at the root (not supported)
                    @test_throws ErrorException Newick.parse(Ttype, line)
                else
                    tree = Newick.parse(Ttype, line)
                    @test Newick.string(tree) == line
                end
            end
        end


    end
end
