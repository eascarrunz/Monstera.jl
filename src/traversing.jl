"""
    pretraverse(f, tree)
    pretraverse(f, branch => node)
    pretraverse(f, node)

Traverse a tree (or subtree) in preorder, calling the function `f` on each branch-node pair
as they are visited.

Returns the return value of the last call of `f`.
"""
function pretraverse(f::Function, branchnode)
    result = f(branchnode)

    for cbranchnode in children(branchnode)
        result = pretraverse(f, cbranchnode)
    end

    return result
end

pretraverse(f::Function, node::AbstractNode) = pretraverse(f, nothing => node)

pretraverse(f::Function, tree::AbstractTree) = pretraverse(f, tree.root)

function _cladesize(branchnode, i)
    i += 1
    for cbranchnode in children(branchnode)
        i = _cladesize(cbranchnode, i)
    end

    return i
end

function _cladesize_outer(branchnode, i)
    i += isouter(last(branchnode))
    for cbranchnode in children(branchnode)
        i = _cladesize_outer(cbranchnode, i)
    end

    return i
end

"""
    cladesize(branch => node, outer=false)
    cladesize(node, outer=false)
    cladesize(tree, outer=false)

Return the number of nodes in a clade, optionally counting only the outer nodes.

The clade size of a tree is computed from the tree's root. Unlike the `size` function,
cladesize counts that are directly or indirectly connected to the root.
"""
cladesize(branchnode, outer=false) =
    outer ? _cladesize_outer(branchnode, 0) : _cladesize(branchnode, 0)

cladesize(node::AbstractNode, outer=false) = cladesize(nothing => node, outer)

cladesize(tree::AbstractTree, outer=false) = cladesize(tree.root, outer)



