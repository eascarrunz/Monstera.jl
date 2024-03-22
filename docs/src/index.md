# Monstera.jl Documentation


```@meta
CurrentModule = Monstera
```


This is an experimental package for phylogenetics in Julia.

These are the main objectives of this package:

- Provide a general interface for basic operations with phylogenetic trees, allowing the 
experimentation with a variety of data structures to represent trees

- Allow as much flexibility as possible:
    - Taxa can live in inner or outer nodes
    - Non-splitting nodes are permitted

- Be ergonomic for interactive use

Only one base tree type is available for now: `RoundaboutTree` (see below).

## Quick Introduction

Let's start off loading a tree in Newick format. This is a little extract of the tree of the
Araceae from Haigh et al. (2022): 

```jldoctest quick_intro
julia> philo_string = "(((Philodendron,Adelonema),(Homalomena,Furtadoa)),(Cercestis,Culcasia),Monstera);";

```
We can create a tree object from that string using the `Newick.parse` method, and specifying
the tree type that we want. Only one tree type is available for now: `RoundaboutTree`. Optionally, we can create a `TaxonSet` from the labels of
the tree.

```jldoctest quick_intro
julia> philo_tree = Newick.parse(RoundaboutTree, philo_string, TaxonSet())
RoundaboutTree: 12 nodes, 11 branches
```
Let's take a look at it. The `textplot` function will print a diagram of the tree:

```julia-repl
julia> textplot(philo_tree)
             ┌───────●(4)Philodendron
     ┌───────○(3)
     │       └───────●(5)Adelonema
┌────○(2)
│    │       ┌───────●(7)Homalomena
│    └───────○(6)
│            └───────●(8)Furtadoa
○(1)
│    ┌───────●(10)Cercestis
├────○(9)
│    └───────●(11)Culcasia
│
└────●(12)Monstera
```
Nodes are represented by small circles, with numerical identifiers in parentheses. Next to that, node labels are shown. Nodes are objects, and they can be retrieved from the tree by
their IDs. For instance, let's take a look at nodes #3 and #4, and their main properties.

```jldoctest quick_intro
julia> philo_tree.nodes[3]
RoundaboutNode #3 (no taxon)
julia> philo_tree.nodes[3].id
3
julia> philo_tree.nodes[3].taxon
0
julia> philo_tree.nodes[3].label
""
julia> philo_tree.nodes[4]
RoundaboutNode #4 (taxon 1) - "Philodendron"
julia> philo_tree.nodes[4].id
4
julia> philo_tree.nodes[4].taxon
1
julia> philo_tree.nodes[4].label
"Philodendron"
```

Nodes can be associated to a taxon with its own numeric ID in a `TaxonSet`, or have no taxon
(interally represented with a `0`). Also, the diagram from `textplot` represents nodes with 
taxa as filled circles. The taxon set is property of the tree, accessible as `philo_tree.taxonset`. A single taxon set object is shared with groups of trees that represent the same
group of taxa.

Branches are have IDs, like nodes, and are accessible with the `branches` property.

```jldoctest quick_intro
julia> philo_tree.branches
11-element Vector{RoundaboutBranch}:
 RoundaboutBranch #1: 1 ○⎯⎯⎯○ 2 (length NaN)
 RoundaboutBranch #2: 2 ○⎯⎯⎯○ 3 (length NaN)
 RoundaboutBranch #3: 3 ○⎯⎯⎯● 4 (length NaN)
 RoundaboutBranch #4: 3 ○⎯⎯⎯● 5 (length NaN)
 RoundaboutBranch #5: 2 ○⎯⎯⎯○ 6 (length NaN)
 RoundaboutBranch #6: 6 ○⎯⎯⎯● 7 (length NaN)
 RoundaboutBranch #7: 6 ○⎯⎯⎯● 8 (length NaN)
 RoundaboutBranch #8: 1 ○⎯⎯⎯○ 9 (length NaN)
 RoundaboutBranch #9: 9 ○⎯⎯⎯● 10 (length NaN)
 RoundaboutBranch #10: 9 ○⎯⎯⎯● 11 (length NaN)
 RoundaboutBranch #11: 1 ○⎯⎯⎯● 12 (length NaN)
```
Little diagrams of the branches show the nodes to which they are attached. We have no branch
lengths in this example, so the "length" value is set to `NaN`.

### Orientation conventions and the `branch => node` notation

This is a central concept in the Monstera package. Tree objects can have an intrinsic 
orientation, so that branches are directed outward from the root node. Monstera also makes 
it possible to traverse the tree in any arbitrary direction by means of special conventions 
in function calls. For instance, if we have a function `f` that traverses a tree, or a 
sub-tree of a tree, we can call `f` in a number of ways to indicate exactly what part of the
 tree to traverse, and in what direction:


| Notation | Meaning |
| --- | --- |
| `f(tree)` | Traverse the entire tree starting from its root node |
| `f(node)` | Traverse the part of the tree that descends from a node, following to the intrinsic orientation of the tree |
| `f(nothing => node)` | Traverse the entire tree as if it were rooted on a node |
| `f(branch => node)` | Traverse the part of the tree that descends from a node, with the parent-child direction indicated by the branch-node pair |

The `textplot` function is one of the many functions that work in this way. It prints a 
diagram of the tree, or a sub-tree, using characters on a terminal window.  Let's see it in  action. First, we create a tree from a Newick string, and we plot it with 
`textplot(tree)`. The tree is plotted from left to right, from its root node.

```julia-repl
julia> tree = Newick.parse(RoundaboutTree, "(((A,B),C),(D,E),H);");

julia> textplot(tree)
             ┌───────○(4)A
     ┌───────○(3)
┌────○(2)    └───────○(5)B
│    │
│    └───────○(6)C
│
○(1) ┌───────○(8)D
├────○(7)
│    └───────○(9)E
│
└────○(10)H
```

As node #1 is the root node, `textplot(tree.nodes[1])` or `textplot(tree.root)` produces an 
identical result.

```julia-repl
julia> textplot(tree.nodes[1])
             ┌───────○(4)A
     ┌───────○(3)
┌────○(2)    └───────○(5)B
│    │
│    └───────○(6)C
│
○(1) ┌───────○(8)D
├────○(7)
│    └───────○(9)E
│
└────○(10)H
```

If we do the same with node #2 instead, we get the corresponding sub-tree.

```julia-repl
julia> textplot(tree.nodes[2])
     ┌───────○(4)A
┌────○(3)
○(2) └───────○(5)B
│
└────○(6)C
```

The `branch => node` notation gives us a more flexible way to plot sub-trees, but we need to
know which branches are connecting which nodes. We get the information we need just by 
showing the branch list on the REPL.

```julia-repl
julia> tree.branches
9-element Vector{RoundaboutBranch}:
 RoundaboutBranch #1: 1 ○⎯⎯⎯○ 2 (length NaN)
 RoundaboutBranch #2: 2 ○⎯⎯⎯○ 3 (length NaN)
 RoundaboutBranch #3: 3 ○⎯⎯⎯○ 4 (length NaN)
 RoundaboutBranch #4: 3 ○⎯⎯⎯○ 5 (length NaN)
 RoundaboutBranch #5: 2 ○⎯⎯⎯○ 6 (length NaN)
 RoundaboutBranch #6: 1 ○⎯⎯⎯○ 7 (length NaN)
 RoundaboutBranch #7: 7 ○⎯⎯⎯○ 8 (length NaN)
 RoundaboutBranch #8: 7 ○⎯⎯⎯○ 9 (length NaN)
 RoundaboutBranch #9: 1 ○⎯⎯⎯○ 10 (length NaN)
```

Since branch #1 connects node #1 (the root) to node #2, 
`textplot(tree.branches[1] => tree.nodes[2])` will print the same sub-tree as in the 
previous example.

```julia-repl
julia> textplot(tree.branches[1] => tree.nodes[2])
     ┌───────○(4)A
┌────○(3)
○(2) └───────○(5)B
│
└────○(6)C
```

With the `branch => node` notation, we can plot sub-trees with a completely different 
orientation, without actually changing the root of the tree. The key is that we consider the
given `branch` as coming from the parent of the `node`. Take a look at the following 
example:

```julia-repl
julia> textplot(tree.branches[2] => tree.nodes[2])
┌────○(6)C
│
○(2)         ┌───────○(8)D
│    ┌───────○(7)
└────○(1)    └───────○(9)E
     │
     └───────○(10)H
```

Branch #2 connects node #2 to node #3, so we are in effect plotting the sub-tree of node #2
as though node #3 were its parent. Hence, node #1, the true root of the tree, appears as
a descendant of node #2, and node #3 is left out.

The `branch => node` notation also gives us the ability to plot the entire tree as though 
any node were the root of the tree. Since a "proper" root node has no parent, we replace the
`branch` with `nothing`, like so:

```julia-repl
julia> textplot(nothing => tree.nodes[2])
     ┌───────○(4)A
┌────○(3)
│    └───────○(5)B
│
○(2)─○(6)C
│
│            ┌───────○(8)D
│    ┌───────○(7)
└────○(1)    └───────○(9)E
     │
     └───────○(10)H
```

The usefulness of having all these options for defining tree orientations might not be 
obvious in common use cases, but it has practical value. For instance, in many cases the
root of the tree might not be biologically meaningful, such as in trees inferred by 
maximum-likelihood without an explicit outgroup.

Also, some important algorithms require 
treating any arbitrary node as the root of the tree, as in maximum-likelihood branch length 
optimisation.

Finally, the `branch => node` notation also a useful *return value* of many
functions. The `children` function, for instance, returns an iterator of `branch => node` 
pairs: instead of just listing the children nodes, we get the corresponding branch and node
of each child.

```julia-repl
julia> collect(children(tree.branches[2] => tree.nodes[2]))
2-element Vector{Pair{RoundaboutBranch{RoundaboutNode}, RoundaboutNode}}:
 RoundaboutBranch #5 => RoundaboutNode #6 (no taxon) - "C"
 RoundaboutBranch #1 => RoundaboutNode #1 (no taxon)
```


## The `RoundaboutTree` type

This type implements a complex data structure based on circular linked lists of branches 
surrounding nodes: "roundabouts" of branches. The use of linked lists makes the construction
 and modification of trees very fast, as it avoids allocations associated to array-based 
 lists. This comes at the cost of complexity and indirections that slow down tree 
 traversals.

Each `RoundaboutBranch` object has a "left" and "right" side. They correspond to the 
"parent" and "child" sides of the branch, following the orientation in which the tree was 
constructed, for instance, from a Newick string.

Each side of a branch contain fields with the following information:
- A reference to the `RoundaboutNode` connected by this side, or `nothing` if there is no node yet (this should only happen in intermediate stages of tree construction or branch swapping)

- A reference to the next `RoundaboutBranch` connected to this side's node. If there are no
 other branches, this field points back to the present branch.

- An instance of the side type of the next `RoundaboutBranch` which is connected to this side's node.

Side types are `RoundaboutLeft`, `RoundaboutRight` and `RoundaboutNoSide`; subtypes of `RoundaboutSide`.

The `RoundaboutNode` objects contain a reference to the last `RoundaboutBranch` that was
connected to them. This provides quick access to both the head and tail of the circular 
list.

The `RoundaboutTree` data structure is a spin on Felsenstein's nodelet-based tree data 
structure (Felsenstein 2004), used in PHYLIP and in the Phylogenetic Likelihood Library used
by RAxML-NG. I replaced the pairs of "nodelets" pointing back to each other with single 
branches that contain two sides. I also made nodes concrete objects, insetad of emerging 
properties of cycles of nodelets. The advantage is that this structure fits a generic tree 
interface better, and it is conceptually simpler and easier to use in interactive 
programming.


