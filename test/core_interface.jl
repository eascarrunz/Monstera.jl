for Ttype in CONCRETE_TREE_TYPES
    @testset "Matching type functions" begin
        Ntype = nodetype(Ttype)
        Btype = branchtype(Ttype)

        @test Ntype <: AbstractNode
        @test Btype <: AbstractBranch
        @test nodetype(Btype) ≡ Ntype
        @test branchtype(Ntype) ≡ Btype
        @test treetype(Ntype) ≡ Ttype
        @test treetype(Btype) ≡ Ttype
    end

    taxonset = TaxonSet()
    tree = Ttype(taxonset, 8)

    @testset "Node properties" begin
        for propname in (:id, :taxon, :label)
            @test hasproperty(tree.nodes[1], propname)
        end
    end

    @testset "Branch properties" begin
        for propname in (:id, :length)
            @test hasproperty(tree.branches[1], propname)
        end
    end

    @testset "Tree properties" begin
        for propname in (:taxonset, :root, :nodes, :branches)
            @test hasproperty(tree, propname)
        end

        @test length(tree.nodes) == 8
        @test length(tree.branches) == 7
    end

    @testset "Neighbours and children" begin
        @test isempty(neighbours(tree.nodes[1]))
        @test isempty(children(tree.nodes[1]))
        @test isempty(children(nothing => tree.nodes[1]))
    
        link!(tree.nodes[1], tree.branches[1], tree.nodes[2])
        @test only(neighbours(tree.nodes[1])) == (tree.branches[1] => tree.nodes[2])
        @test only(children(tree.nodes[1])) == (tree.branches[1] => tree.nodes[2])
        @test only(neighbours(tree.nodes[2])) == (tree.branches[1] => tree.nodes[1])
        @test isempty(children(tree.nodes[2]))
        @test only(children(nothing => tree.nodes[2])) == (tree.branches[1] => tree.nodes[1])

        link!(tree.nodes[2], tree.branches[2], tree.nodes[3])
        @test only(neighbours(tree.nodes[1])) == (tree.branches[1] => tree.nodes[2])
        @test length(neighbours(tree.nodes[2])) == 2
        @test only(neighbours(tree.nodes[3])) == (tree.branches[2] => tree.nodes[2])
        @test nodes_flanking(tree.branches[1]) == (tree.nodes[1], tree.nodes[2])
        @test nodes_flanking(tree.branches[2]) == (tree.nodes[2], tree.nodes[3])

        @test ! hasparent(tree.nodes[1])
        @test hasparent(tree.nodes[2])
        @test hasparent(tree.nodes[3])
    end
end