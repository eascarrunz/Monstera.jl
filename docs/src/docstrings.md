# Functions
```@meta
CurrentModule = Monstera
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
Newick.parse
Newick.read
Newick.string
Newick.write
Newick.WriterSettings
Newick.ReaderSettings
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
graft!
```

## Bipartitions
```@docs
TaxonBipartition
bipartitions
is_singleton
are_compatible
bipartition_table
```

## Tree distances
```@docs
distance
```