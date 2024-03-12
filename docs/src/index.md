# Monstera.jl Documentation

This is an experimental package for phylogenetics in Julia.

These are the main objectives of this package:

- Provide a general interface for basic operations with phylogenetic trees, allowing the 
experimentation with a variety of data structures to represent trees

- Allow as much flexibility as possible:
    - Taxa can live in inner or outer nodes
    - Non-splitting nodes are permitted

- Be ergonomic for interactive use

## Quick Introduction

Let's start off loading a tree in Newick format. This is a little extract of the tree of the
Araceae from Haigh et al. (2022): 

```julia-repl
julia> using Monstera
julia> philo_string = "(((Philodendron,Adelonema),(Homalomena,Furtadoa)),(Cercestis,Culcasia),Monstera);"
```
We can create a tree object from that string using the `Newick.parse` method, and specifying
the tree type that we want. Only one tree type is available for now: `RoundaboutTree`. Optionally, we can create a `TaxonSet` from the labels of
the tree.

```julia-repl
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

```julia-repl
julia> philo_tree.nodes[3]
RoundaboutNode #3 (no taxon)
julia> philo_tree.nodes[3].id
0x0003
julia> philo_tree.nodes[3].taxon
0x0000
julia> philo_tree.nodes[3].label
""

julia> philo_tree.nodes[4]
RoundaboutNode #4 (taxon 1) - "Philodendron"
julia> philo_tree.nodes[4].id
0x0004
julia> philo_tree.nodes[4].taxon
0x0001
julia> philo_tree.nodes[4].label
"Philodendron"
```

Nodes can be associated to a taxon with its own numeric ID in a `TaxonSet`, or have no taxon
(interally represented with a `0x0000`). Also, the diagram from `textplot` represents nodes with 
taxa as filled circles. The taxon set is property of the tree, accessible as `philo_tree.taxonset`. A single taxon set object is shared with groups of trees that represent the same
group of taxa.

Branches are have IDs, like nodes, and are accessible with the `branches` property.

```julia-repl
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

### The `branch => node` notation

The `branch => node` notation, i.e. a `Pair` containing a branch object and a node object, is a central concept in the Monstera package. It is used to indicate an origin and direction. It's easier to explain by example. Many functions take `branch => node` pairs as alternatives to trees. For instance, this is `textplot` from node #2:

```julia-repl
julia> textplot(philo_tree.nodes[2])
     ┌───────●(4)Philodendron
┌────○(3)
│    └───────●(5)Adelonema
○(2)
│    ┌───────●(7)Homalomena
└────○(6)
     └───────●(8)Furtadoa
```


Say, for instance, that we want 

`Pair` objects with a branch as the first element and a node as the second element (`branch => node`) are very common in Monstera. Let's create a tree and see an example with the `children` function.

```julia-repl
# Create the tree
julia> example_newick = "((A,B)C,(D,E)F,(G,(H,I)J)K)L;"
"((A,B)C,(D,E)F,(G,(H,I)J)K)L;"
julia> example_tree = Newick.parse(RoundaboutTree, example_newick, TaxonSet())
RoundaboutTree: 12 nodes, 11 branches

# Take a look at the nodes
julia> example_tree.nodes |> collect
 RoundaboutNode #1 (no taxon) - "L"
 RoundaboutNode #2 (no taxon) - "C"
 RoundaboutNode #3 (taxon 1) - "A"
 RoundaboutNode #4 (taxon 2) - "B"
 RoundaboutNode #5 (no taxon) - "F"
 RoundaboutNode #6 (taxon 3) - "D"
 RoundaboutNode #7 (taxon 4) - "E"
 RoundaboutNode #8 (no taxon) - "K"
 RoundaboutNode #9 (taxon 5) - "G"
 RoundaboutNode #10 (no taxon) - "J"
 RoundaboutNode #11 (taxon 6) - "H"
 RoundaboutNode #12 (taxon 7) - "I"
```
And take a quick look at the node numbers:

```julia-repl

```

We can check out the neighbours of a specific node, node #10 labelled `"J"`, with the `neighbours` function. It returns a lazy iterator, so we will collect it right away.

```julia-repl
julia> neighbours(example_tree.nodes[10]) |> collect
3-element Vector{Pair{RoundaboutBranch{RoundaboutNode}, RoundaboutNode}}:
 RoundaboutBranch #10 => RoundaboutNode #11 (taxon 6) - "H"
 RoundaboutBranch #11 => RoundaboutNode #12 (taxon 7) - "I"
  RoundaboutBranch #9 => RoundaboutNode #8 (no taxon) - "K"
```


Now let's take a closer look at the topology of the tree. What are the neighbours of the node labelled "Protostomia"? Let's fetch that node (we already know it's node #3):

```julia-repl
julia> node_protostomia = animal_tree.nodes[findfirst(x -> x.label == "Protostomia", animal_tree.nodes)]
RoundaboutNode #3 (no taxon) - "Protostomia"
```

And now we can use the `neighbours` function. It returns a lazy iterator, so we'll collect it right away.

```julia-repl
julia> neighbours(node_protostomia) |> collect
3-element Vector{Pair{RoundaboutBranch{RoundaboutNode}, RoundaboutNode}}:
 RoundaboutBranch #3 => RoundaboutNode #4 (no taxon) - "Arthropoda"
 RoundaboutBranch #6 => RoundaboutNode #7 (no taxon) - "Mollusca"
 RoundaboutBranch #2 => RoundaboutNode #2 (no taxon) - "Bilateria"
```

As each neighbouring node is connected to `node_protostomia` through a branch, `neighbours` gives us `Pair`s of branches and nodes. This `branch => node` convention is central to `Monstera`, because it allows us to work with specific orientations of the relationships in the tree. For instance, let's look at the children of `node_protostomia` (this also returns a lazy iterator).

```julia-repl
children(node_protostomia) |> collect
2-element Vector{Any}:
 RoundaboutBranch #3 => RoundaboutNode #4 (no taxon) - "Arthropoda"
 RoundaboutBranch #6 => RoundaboutNode #7 (no taxon) - "Mollusca"
```

As expected, the nodes labelled "Arthropoda" and "Mollusca" are children of `node_protostomia`, per the rooting of the tree in the original Newick string. But we can get the nodes that would be children of `node_protostomia` as seen from *any* direction; just use the `branch => node` convention.

```julia-repl
# Children as coming from branch #3
julia> children(animal_tree.branches[3] => node_protostomia) |> collect
2-element Vector{Pair{RoundaboutBranch{RoundaboutNode}, RoundaboutNode}}:
 RoundaboutBranch #6 => RoundaboutNode #7 (no taxon) - "Mollusca"
 RoundaboutBranch #2 => RoundaboutNode #2 (no taxon) - "Bilateria"

# Children as coming from branch #6
julia> children(animal_tree.branches[6] => node_protostomia) |> collect
2-element Vector{Pair{RoundaboutBranch{RoundaboutNode}, RoundaboutNode}}:
 RoundaboutBranch #2 => RoundaboutNode #2 (no taxon) - "Bilateria"
 RoundaboutBranch #3 => RoundaboutNode #4 (no taxon) - "Arthropoda"
```
This is very handy, because trees are not always meaningfully rooted. For instance, it is easy to write functions to traverse the tree from any point, in any direction.

```julia-repl
"""
Traverse the nodes of the tree in preorder, printing the node ID and label.
"""
function my_traverse(branch_node_pair)
    # Get the node from the `branch_node` pair
    node = last(branch_node_pair)

    # Print node ID and label
    println(node.id, ": ", node.label)

    # Use recursion to repeat to visit the children nodes
    for child_branch_node_pair in children(branch_node_pair)
        my_traverse(child_branch_node_pair)
    end
end
```

If we traverse the tree from `node_protostomia` as coming from branch #2:

```julia-repl
julia> my_traverse(animal_tree.branches[2] => node_protostomia)
3: Protostomia
4: Arthropoda
5: 🦋
6: 🦀
7: Mollusca
8: 🐌
9: 🦑
```

Or from branch 3

```julia-repl
julia> my_traverse(animal_tree.branches[3] => node_protostomia)
3: Protostomia
7: Mollusca
8: 🐌
9: 🦑
2: Bilateria
10: Vertebrata
11: 🦆
12: 🐟
1: Animalia
13: 🪼
```

The `Newick.string` method is an example of the use of the `branch => node` convention.



