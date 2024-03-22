#=
This file contains core functions, including function prototypes that defines the signatures
and docstrings for the functions of the interface for tree (and node and branch) types.

Every tree (and node and branch) type must implement the methods that are defined
tautologically in this file. The tautological definitions exist to bind docstrings to the
methods and to make the type signatures explicit.
=#

abstract type AbstractTree end

"""
The `AbstractNode` type is the supertype of node types. Subtypes of `AbstractNode` must
implement the public properties listed below.

# Public properties
- `id`: integer number that uniquely identifies a node in a tree. Not user-settable.
- `taxon`: integer number that identifies the taxon associated to a node. User-settable.
- `label`: string to display for this node when printing, plotting, etc. User-settable.
"""
abstract type AbstractNode end


"""
The `AbstractBranch` type is the supertype of branch types. Subtypes of `AbstractBranch`
must implement the public properties listed below.

# Public properties
- `id`: integer number that uniquely identifies a node in a tree. Not user-settable.
- `length`: branch length as a `Float64`. Branches with `NA` values are considered to
have no length. User-settable.
"""
abstract type AbstractBranch end


"""
    treetype(x) -> type <: AbstractTree
    nodetype(x) -> type <: AbstractNode
    branchtype(x) -> type <: AbstractBranch

Return the node, branch, or tree type that matches another type.

This is an auxiliary function used to select matching types in constructors.
"""
function treetype end, function nodetype end, function branchtype end

treetype(::Type{AbstractNode}) = AbstractTree
treetype(::Type{AbstractBranch}) = AbstractTree
nodetype(::Type{AbstractTree}) = AbstractNode
nodetype(::Type{AbstractBranch}) = AbstractNode
branchtype(::Type{AbstractTree}) = AbstractBranch
branchtype(::Type{AbstractNode}) = AbstractBranch


const AbstractBranchNodePair = Pair{<:AbstractBranch,<:AbstractNode}
const AbstractNothingNodePair = Pair{<:Nothing,<:AbstractNode}


"""
    size(tree)

Return the number of branches and nodes in a tree.
"""
Base.size(tree::AbstractTree) = (length(tree.branches), length(tree.nodes))



"""
    neighbours(node) -> iterator

Return an iterator of the branches and nodes that are connected to a given `node`.

See also [`children`](@ref).
"""
neighbours(node::AbstractNode) = neighbours(node)


"""
    children(branch => node) -> iterator
    children(node) -> iterator

Return an iterator of the branches and nodes that are children of a `node`.

See also [`neighbours`](@ref).
"""
children(branchnode::AbstractBranchNodePair) = children(branchnode)
children(nothingnode::Pair{Nothing,<:AbstractNode}) = neighbours(last(nothingnode))
children(node::AbstractNode) = Iterators.filter(
    x -> node ≡ first(nodes_flanking(first(x))),
    neighbours(node)
    )


"""
    isinner(node) -> Bool

Check whether a node is inner (a.k.a. "internal").

Inner nodes have a degree of 2 or greater.

See also [`isouter`](@ref).
"""
isinner(node::AbstractNode) = ! isouter(node)


"""
    isouter(node) -> Bool

Check whether a node is outer (a.k.a. "external", "tip", or "leaf").

Outer nodes always have a degree of 1. This can include the tree root and isolated nodes.

See also [`isinner`](@ref).
"""
isouter(node::AbstractNode) = length(neighbours(node)) < 2


"""
    nodes_flanking(branch)

Return a tuple with the nodes connected to a branch.

When nodes are missing from one side of the branch, `nothing` will take their place in the tuple.

See also [`branch_between`](@ref).
"""
nodes_flanking(branch::AbstractBranch) = nodes_flanking(branch)


"""
    branch_between(node1, node2)

Return the branch that connects node1 to node2.

See also [`nodes_flanking`](@ref).
"""
function branch_between(node1::AbstractNode, node2::AbstractNode)
    for (branch, node) in neighbours(node1)
        node ≡ node2 && return branch
    end

    error("node1 is not connected to node2")
end



"""
    link!(node1, branch, node2)

Connect `node1` to `node2` through a `branch`.

Linking order is left-to-right, so that `node1` is connected on the left side of the
branch and `node2` is connected on the right side of the branch. The new link is added to
the tail of the link list of each node.

The `graft!` function can also be used to link nodes without need to pass a branch object.

See also [`graft!`, `unlink!`](@ref)
"""
link!(node1::AbstractNode, branch::AbstractBranch, node2::AbstractNode)


"""
unlink!(branch, node)
unlink!(branch => node)

Disconnect a branch from a node.

See also [`link!`](@ref)
"""
unlink!(branch::AbstractBranch, node::AbstractNode) = unlink!(branch, node)
unlink!(branchnode::Pair{TB,TN}) where {TB<:AbstractBranch, TN<:AbstractNode} =
    unlink!(first(branchnode), last(branchnode))


"""
    unlink!(node1, node2) -> Branch

Disconnect node1 from node2 and return the branch that was between them.
"""
function unlink!(node1::AbstractNode, node2::AbstractNode)
    branch = branch_between(node1, node2)
    unlink!(branch, node1)
    unlink!(branch, node2)

    return branch
end


"""
    hasparent(node)

Check whether a node has a parent.
"""
function hasparent(node)
    for neighbranchnode in neighbours(node)
        first(nodes_flanking(first(neighbranchnode))) ≢ node && return true
    end

    return false
end


"""
    parent(node)
    parent(branch => node)

Return the parent of a node.

Throws an error if the node has no parent.
"""
Base.parent(node::AbstractNode) = 
    only(Iterators.filter(x -> first(nodes_flanking(first(x))) ≢ node, neighbours(node)))

function Base.parent(branchnode::AbstractBranchNodePair)
    branch, node = branchnode
    lnode, rnode = nodes_flanking(branch)

    return node ≡ lnode ? rnode : lnode
end


"""
    swapnodes!(branch1 => node1, branch2 => node2) -> nothing

Exchange nodes between two pairs of nodes and branches.
"""
function swapnodes!(
    branchnode1::Pair{<:AbstractBranch,<:AbstractNode},
    branchnode2::Pair{<:AbstractBranch,<:AbstractNode}
    )
    branch1, node1 = branchnode1
    branch2, node2 = branchnode2

    node1a, node1b = nodes_flanking(branch1)
    node2a, node2b = nodes_flanking(branch2)

    pnode1 = node1 ≡ node1a ? node1b : node1a
    pnode2 = node2 ≡ node2a ? node2b : node2a

    unlink!(pnode1, node1)
    unlink!(pnode2, node2)

    link!(pnode1, branch1, node2)
    link!(pnode2, branch2, node1)
    
    return nothing
end



##### NOT USED YET, OR NOT IMPLEMENTED YET #####

"""
    newnode!(tree, n=1)
    newbranch!(tree, n=1)

Create and return `n` new branches or nodes in a tree.
"""
(newnode!(tree::AbstractTree, n) = newnode!(tree,n)),
(newbranch!(tree::AbstractTree, n) = newbranch!(tree, n))

newnode!(tree::AbstractTree) = only(newnode!(tree, 1))
newbranch!(tree::AbstractTree) = only(newbranch!(tree, 1))


"""
    graft!(tree, graft_node, receiver_node, brlength=NaN)
    graft!(tree, graft_node, receiver_branch, newinner=false)

Graft a node (and the subtree to which it is connected) on a receiver branch or node.

Returns `nothing`, or a new branch-node pair when the `newinner` option is used (see below).

## Grafting on a node
To graft a node on another node is to link the two with a new branch. The length of the 
new branch can be set with the `brlength` argument. The function returns the new branch.

## Grafting on a branch
Say that you have a branch `br1` (1, in the figure below), connected to nodes `a` and `d`.
To graft a node `g` directly on branch `br1` means that `g` is placed in between `a` and
`d`, which results in the creation of a new branch (2, in the figure below). The function 
returns `nothing`.
                                      

```

      (a)              (e)                                   (a)                     (e)                    
        \\              /                                       \\                     /                      
         \\            /                                         \\                   /                       
          \\          /                                           \\                 /                        
           \\        /            `graft!(tree, g, br1)`           \\               /                         
(g)        (a)----(d)         -------------------------->         (a)----(g)----(d)                         
           /   1    \\                                             /   1      2    \\                         
          /          \\                                           /                 \\                        
         /            \\                                         /                   \\                       
        /              \\                                       /                     \\                      
      (c)              (f)                                   (c)                     (f)
```
"""
function graft!(
    tree::AbstractTree,
    graft_node::AbstractNode,
    receiver_node::AbstractNode,
    brlength=NaN
    )
    newbr = newbranch!(tree)
    link!(receiver_node, newbr, graft_node)
    newbr.length = brlength

    return newbr
end

function graft!(
    tree::AbstractTree,
    graft_node::AbstractNode,
    receiver_branch::AbstractBranch
    )
    br2 = newbranch!(tree)
    a, d = nodes_flanking(receiver_branch)
    unlink!(receiver_branch, d)
    link!(a, receiver_branch, graft_node)
    link!(graft_node, br2, d)

    return br2
end


"""
    delete!(tree, node...)
    delete!(tree, branch...)

Delete nodes or branches from a tree.

This operation removes elements from the list of nodes or branches of a tree.
"""


"""
    ungraft!(tree, branch => node)

Separate a node with a branch from the rest of the tree without leaving other nodes 
disconnected.
"""
function ungraft!(tree::AbstractTree, branchnode::Pair{<:AbstractBranch,<:AbstractNode})
    branch, node = branchnode
    d = degree(node)
    if d > 3
        o = newnode(tree)
        for (nbr, nnode) in neighbours(node)
            nbr ≡ branch && continue
            left_node, _ = nodes_flanking(nbr)
            unlink!(nbr, node)
            if left_node ≡ node
                link!(o, nbr, nnode)
            else
                link!(nnode, nbr, o)
            end
        end
    elseif d > 2
        (nbr1, nnd1), (nbr2, nnd2) =
            Iterators.filter(x -> first(x) ≢ branch, neighbours(node))    # Lazy filter
        unlink!(nbr1, node)
        unlink!(nbr2, node)
        link!(nnd1, nbr1, nnd2)
        delete!(tree, nbr2)
    elseif d > 1
        left_node, right_node = nodes_flanking(branch)
        left_node ≡ node ? unlink!(branch, right_node) : unlink!(branch, left_node)
    end

    error("cannot ungraft, as $(node) is only linked to $(branch)")
end


# Optional methods

"""
    cleanup!(tree)

Remove unconnected nodes and branches and re-index the connected elements.
"""
function cleanup! end


#TODO: implement `collapse!`
"""
    collapse!(tree, branch => node)

Eliminate a `branch` and `node` in the `tree`, and reconnect their children to their parent.

In the example below, we collapse branch "1" and node "b".
                                               
                                               
  (c)       (d)                                 
    \\       /                                   
     \\     /                                    
     3\\  4/                        (c)  (d)  (e)
       \\ /                           \\   |   /  
       (b)       (e)                  \\ 4|  /   
         \\       /                    3\\ |2/    
          \\     /      -------->        \\|/     
          1\\  2/                        (a)     
            \\ /                          |      
            (a)                          |      
             |                           |      
             |                           ⃨      
             |                                  
             ⃨                                 
                                                
"""
function collapse!(tree::AbstractTree, branchnode::AbstractBranchNodePair)
    branch, nodeb = branchnode
    nodea = parent(branchnode)
    unlink!(branch, nodea)
    cnodeb = children(branchnode)
    for (cbr, cnd) in cnodeb
        unlink!(cbr, cnd)
        unlink!(cbr, nodeb)
        link!(nodea, cbr, cnd)
    end

    return nothing
end


#TODO: implement `prune!`
"""
    prune!(branch1, node, branch2)


"""
function prune! end

