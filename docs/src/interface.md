# Interface

Tree data structures must implement their own types for trees, nodes, and branches, which ought to be subtypes of `AbstractTree`, `AbstractNode`, and `AbstractBranch`, respectively. As an example, let's call them `MyTree`, `MyNode`, and `MyBranch`.

```julia
struct MyBranch <: AbstractBranch
    # ...<Implementation details>...
end

struct MyNode <: AbstractNode
    # ...<Implementation details>...
end

struct MyTree <: AbstractTree
    # ...<Implementation details>...
end
```

The `Monstera` interface requires us to indicate that those types match each other by implementing the following methods:

```julia
treetype(::Type{MyNode}) = MyTree
treetype(::Type{MyBranch}) = MyTree

nodetype(::Type{MyTree}) = MyNode
nodetype(::Type{MyBranch}) = MyNode

branchtype(::Type{MyTree}) = MyBranch
branchtype(::Type{MyNode}) = MyBranch
```

## Node types

| Required public properties | Default value | Description | User settable |
|:--------------------|:--------------|:------------|:--------------|
| `id::PhyloInd` | `0` (should be set at tree construction) | Numerical identifier of the node | No |
| `taxon::PhyloInd` | `0` (= no taxon) | Numerical identifier of the taxon associated to the node | Yes |
| `label::String` | `""` | Label for display | Yes |


## Branch types

| Required public properties | Default value | Description | User settable |
|:--------------------|:--------------|:------------|:--------------|
| `id::PhyloInd` | `0` (should be set at tree construction) | Numerical identifier of the branch | No |
| `length::Float64` | `NaN` | Length of the branch | Yes |


## Tree types

| Required public properties | Default value | Description | User settable |
|:--------------------|:--------------|:------------|:--------------|
| `taxonset::Union{Nothing,TaxonSet}` | `nothing` | Set of taxa represented in the tree | Yes |
| `root::Union{Nothing,MyNode}` | `nothing` | Access point to the tree | Yes |
| `nodes` | - | Collection with settable indices listing all the nodes associated to the tree | No |
| `branches` | - | Collection with settable indices listing all the branches associated to the tree | No |

Tree types must have a constructor with the signature `MyTree(taxonset::Union{Nothing,TaxonSet}, n)`. The constructor initialises the tree with `n` nodes (`n` > 0) and `n - 1` branches.


## Mandatory methods

Mandatory methods pertaining to the interconnections among branches and nodes follow a.

| Required methods | Brief description |
|:-----------------|:------------------|
| `neighbours(node)` | Return a `branch => node` iterator of the neighbours of a node |
| `children(branch => node)` | Return a `branch => node` iterator of the children of `node` as seen coming from `branch` |
| `hasparent(node)` | Check whether the node as an intrinsic parent |
| `parent(branch => node)` | Return the parent of a `node`, as seen coming  from `branch` |
| `nodes_flanking(branch)` | Return a tuple with the nodes connected to a branch |
| `link!(node1, branch, node2)` | Link two nodes through a branch, in left-to-right direction |
| `unlink!(node1, node2)` | Unlink two nodes and return the branch that was connecting them |

