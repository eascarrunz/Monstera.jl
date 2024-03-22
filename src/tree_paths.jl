function _find_path(branchnode, target, path_nodes, path_branches, found_node1)
    branch, node = branchnode
    found_node1 = found_node1 | (node ≡ target)

    for cbranchnode in children(branchnode)
        found_node1 && break
        
        found_node1 =
        _find_path(cbranchnode, target, path_nodes, path_branches, found_node1)
    end

    if found_node1
        push!(path_nodes, node)
        push!(path_branches, branch)
    end

    return found_node1
end



"""
    find_path(node1, node2)

Find the path from `node1` to `node2` in a tree.

Returns a tuple with the path as a vector of nodes and a vector of branches. Raises an error
if no path is found.

# Examples

```jldoctest
julia> tree = Newick.parse(RoundaboutTree, "((5,2),((4,1),3));");

# julia> textplot(tree, outerlabels=false)
#      ┌───────○(3)
# ┌────○(2)
# │    └───────○(4)
# ○(1)
# │            ┌───────○(7)
# │    ┌───────○(6)
# └────○(5)    └───────○(8)
#      │
#      └───────○(9)

# Find the path between node 4 and node 6

julia> path_nodes, path_branches = find_path(tree.nodes[4], tree.nodes[6]);

julia> path_nodes
5-element Vector{RoundaboutNode}:
 RoundaboutNode #4 (no taxon) - "2"
 RoundaboutNode #2 (no taxon)
 RoundaboutNode #1 (no taxon)
 RoundaboutNode #5 (no taxon)
 RoundaboutNode #6 (no taxon)

julia> path_branches
4-element Vector{RoundaboutBranch}:
 RoundaboutBranch #3: 2 ○⎯⎯⎯○ 4 (length NaN)
 RoundaboutBranch #1: 1 ○⎯⎯⎯○ 2 (length NaN)
 RoundaboutBranch #4: 1 ○⎯⎯⎯○ 5 (length NaN)
 RoundaboutBranch #5: 5 ○⎯⎯⎯○ 6 (length NaN)
```
"""
function find_path(node1::N, node2::N) where N <: AbstractNode
    B = branchtype(N)
    path_nodes = Vector{N}()
    path_branches = Vector{B}()

    if node1 ≢ node2
        found_node1 = false
        for branchnode in neighbours(node2)
            found_node1 =
                _find_path(branchnode, node1, path_nodes, path_branches, found_node1)

            if found_node1
                push!(path_nodes, node2)
                break
            end
        end

        found_node1 || error("No path found between $node1 and $node2")
    end

    return (nodes = path_nodes, branches = path_branches)
end


