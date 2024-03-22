
clonable_fields(::Type{<: AbstractNode}) = (:id, :taxon, :label)
clonable_fields(::Type{<: AbstractBranch}) = (:id, :length)
clonable_fields(::Type{<: AbstractTree}) = (:taxonset, )


function clone_topology(branchnode, tree2)
    _, node = branchnode

    for cbranchnode in children(branchnode)
        cbranch, cnode = cbranchnode
        link!(tree2.nodes[node.id], tree2.branches[cbranch.id], tree2.nodes[cnode.id])

        clone_topology(cbranchnode, tree2)        
    end

    return nothing
end


"""
    clone(tree)

Create a copy of a tree.

The copied tree has the same topology and properites, but its contents are fully 
independent. Only the taxon set is shared between the original tree and the clone.
"""
function clone(tree1::T) where T <: AbstractTree
    tree2 = T(length(tree1.nodes))

    clone_topology(nothing => tree1.root, tree2)

    for field in clonable_fields(T)
        setfield!(tree2, field, getfield(tree1, field))
    end

    for field in clonable_fields(nodetype(T))
        for (node1, node2) in zip(tree1.nodes, tree2.nodes)
            setfield!(node2, field, getfield(node1, field))
        end
    end

    for field in clonable_fields(branchtype(T))
        for (branch1, branch2) in zip(tree1.branches, tree2.branches)
            setfield!(branch2, field, getfield(branch1, field))
        end
    end

    tree2.root = tree2.nodes[tree1.root.id]

    return tree2
end


function _randtree(T::Type{<:AbstractTree}, n, rooted)
    N = 2 * n - 2 + rooted  # Total number of nodes
    tree = T(N)
    buds = Set{nodetype(T)}()
    i = 0    # Last node connected to the tree

    if n < 2
        push!(buds, tree.nodes[1])
        
        return tree, buds
    end

    rooted = rooted | (n < 4) 

    link!(tree.nodes[1], tree.branches[1], tree.nodes[2])
    link!(tree.nodes[1], tree.branches[2], tree.nodes[3])
    if rooted
        push!(buds, tree.nodes[2:3]...)
        i = 3
    else
        N < 3 && error("cannot create unrooted tree with fewer than 3 outer nodes")
        link!(tree.nodes[1], tree.branches[3], tree.nodes[4])
        push!(buds, tree.nodes[2:4]...)
        i = 4
    end

    j = i - 1             # Last branch connected to the tree

    while i < N
        i += 2
        j += 2
        bud = rand(buds)
        childnode1, childnode2 = tree.nodes[(i - 1):i]
        childbranch1, childbranch2 = tree.branches[(j - 1):j]
        link!(bud, childbranch1, childnode1)
        link!(bud, childbranch2, childnode2)
        delete!(buds, bud)
        push!(buds, childnode1, childnode2)
    end

    tree.root = tree.nodes[1]

    return tree, buds
end



"""
    randtree(T, taxonset, rooted=true)
    randtree(T, n, rooted=true)

Create a random binary tree of type `T` with the given taxon set, or the given number `n` of
 outer nodes.

If `rooted` is `false`, the tree will be unrooted, i.e. it will have a basal trichotomy
(this requires at least 3 outer nodes).

```jldoctest
Random.seed!(1)
textplot(randtree(RoundaboutTree, 5, false))

# output

┌────○(2)
│
│    ┌───────○(5)
○(1)─○(3)
│    └───────○(6)
│
│    ┌───────○(7)
└────○(4)
     └───────○(8)
```
"""
function randtree(T::Type{<:AbstractTree}, taxonset::TaxonSet, rooted=true)::T
    n = length(taxonset)
    tree, buds = _randtree(T, n, rooted)
    for (id, bud) in enumerate(buds)
        bud.taxon = Int32(id)
        bud.label = key(taxonset, id)
    end

    tree.taxonset = taxonset

    return tree
end

randtree(T::Type{<:AbstractTree}, n, rooted=true)::T = first(_randtree(T, n, rooted))

