struct Patristic end

function _patristic(branchnode, node2, d, found_node2)
    branch, node = branchnode

    if found_node2 || node ≡ node2
        found_node2 = true
    else
        for cbranchnode in children(branchnode)
            d, found_node2 = _patristic(cbranchnode, node2, d, found_node2)
        end
    end

    d += found_node2 && ! isnothing(branch) && branch.length

    return d, found_node2
end


"""
Return the patristric distance between two nodes.
"""
function patristic(node1::T, node2::T) where T <: AbstractNode
    d, found_node2 = _patristic(nothing => node1, node2, 0.0, false)

    found_node2 || error("a path from $node1 to $node2 does not exist")

    return d
end


function patristic(tree::T, taxaonly=true) where T <: AbstractTree
    if taxaonly
        n = length(tree.taxonset)

        targets = Matrix{Union{Nothing,nodetype(T)}}(nothing, n, n)
        fill!(targets, nothing)

        for node in tree.nodes
            if node.taxon > 0
                targets[node.taxon] = node
            end
        end
    else
        n = length(tree.nodes)
        targets = tree.nodes
    end

    D = zeros(n, n)

    for j in 1:n
        if isnothing(targets[j])
            D[j, :] = NaN
            continue
        end
        
        for i in 1:j
            D[i, j] = patristic(targets[i], targets[j])
        end
    end

    return D
end
