# Functions
```@meta
CurrentModule = Monstera
```

```@contents
Pages = ["docstrings.md"]
```

## Managing sets of taxa
```@docs
TaxonSet
get(::TaxonSet, ::String)
key(::TaxonSet, ::Any...)
haskey
keys(::TaxonSet)
add!(::TaxonSet, ::String...)
```

## Creating trees
```@docs
clone
randtree
```

## Reading and writing Newick trees
```@docs
Newick
Newick.parse
Newick.read
Newick.string
Newick.write
Newick.WriterSettings
Newick.ReaderSettings
```

## Plotting trees
```@docs
textplot
```

## Abstract types
```@docs
AbstractNode
AbstractBranch
treetype
nodetype
branchtype
```

## Branches and nodes
```@docs
neighbours
children
isinner
isouter
nodes_flanking
branch_between
hasparent
parent
link!
unlink!
newnode!
newbranch!
swapnodes!
graft!
```

## Traversal
```@docs
pretraverse
cladesize
find_path
```

## Bipartitions
```@docs
TaxonBipartition
bipartitions
is_singleton
are_compatible
bipartition_record
```

## Tree distances
```@docs
distance
```
