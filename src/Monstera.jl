#=
# Variable name conventions

Try not to abbreviate variable names any further than the following:

- Objects
    - `tree`: Tree
    - `node`: Node
        - `cnode`: Child node
        - `pnode`: Parent node
        - `node1`, `node2`, `node3`, ... : Any other nodes
    - `br`: Branch
    - `brnode`: Branch-node pair (Pair{<:AbstractBranch,<:AbstractNode})
        - `cbrnode`: Child branch-node pair
        - `pbrnode`: Parent branch-node pair
        - `brnode1`, `brnode2`, `brnode3`, ... : Any other branch-node pair
    - `bp`: Bipartition
    - `taxonset`: Set of taxa
    - `id`: ID
- Quantitites
    - `ntree`: Number of trees
    - `nnode`: Number of nodes
    - `nbr`: Number of branches
    - `ntax`: Number of taxa
    - `nbp`: Number of taxon bipartitions
    - `nc`: Number of children
    - `deg`: Degree of a node (number of neighbours)
    - `nk`: Number of characters (or sites)

Use underscores to improve readiblity, per Julia's style, but use it consistently among
related variables.

When it's needed to distinguish between a specialised collection type and a simple vector, 
use the `list` postfix for the specialised collection and `vec` for the vector.

These conventions are new and **I've just started using them in the code**.
=#

module Monstera

using DataStructures

pubpropertyerror(T::Type, key::Symbol) =
    error("`$(key) is not a public property of type $(T)")
setpropertyerror(T::Type, key::Symbol) =
    error("`$(key) cannot be set by the user for objects of type $(T)")

# Linear assignment problem solver for the clustering information distance
include("util/lap.jl")

include("taxa.jl")
export TaxonSet, add!, key

include("tree_interface.jl")
export treetype, nodetype, branchtype
export AbstractTree, AbstractNode, AbstractBranch
export isinner, isouter
export neighbours, children, parent, branch_between, nodes_flanking
export newnode!, newbranch!
export link!, unlink!, graft!, ungraft!, collapse!, swapnodes!

include("show.jl")

# Tree data structure based on circular lists of branches
include("RoundaboutTree/core.jl")
export RoundaboutTree, RoundaboutNode, RoundaboutBranch

include("RoundaboutTree/linking.jl")

include("traversing.jl")
export pretraverse, cladesize

include("plotting.jl")
export textplot

include("tree_builders.jl")
export clone, randtree

include("tree_paths.jl")
export find_path

include("bipartitions.jl")
export TaxonBipartition, bipartitions, bipartition_table, are_compatible, is_singleton

include("cluster_information.jl")

include("tree_distances.jl")
export distance

"""
Newick

Submodule for reading and writing trees in Newick format.

The functions for reading and writing Newick trees from and to files and IO streams are 
[`Newick.read`](@ref) and [`Newick.write`](@ref); their use is modelled after `read` and 
`write`from the `CSV` package.

Reading and writing Newick strings in Julia is done with [`Newick.parse`](@ref) and 
[`Newick.string`](@ref).
"""
module Newick
    using ..Monstera
    using ..Monstera: Int32

    include("Newick/write.jl")
    include("Newick/read.jl")
end

export Newick

end
