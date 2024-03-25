# To-Do

## Repo

- [ ] Organise in branches
- [ ] CI: trigger testing and doc builds selectively

## Testing

- [ ] Improve test coverage:
  - [ ] Core
  - [ ] Newick
  - [ ] Bipartitions
  - [ ] Distances
- [ ] Use doctests

## Documentation

- [ ] Organise docstring pages
- [ ] Add examples (notebooks?)

## Functionality

- [ ] `tree[branch_id => node_id]` syntactic sugar
- [ ] Core functions
  - [ ] `reroot!` (change intrinsic branch orientation)
  - [ ] `collapse!`
  - [ ] `graft!`
  - [ ] `ungraft`
  - [ ] `rebuild!` (remove unlinked elements and re-assign ID numbers)
- [ ] Tree builders
  - [x] Random addition (`randtree`)
    - [ ] Replace method to extend `Base.rand` instead?
  - [ ] Balanced
  - [ ] Unbalanced
  - [ ] Star
- [ ] Implement neighbour-joining
  - [x] Naive method
  - [ ] Dynamic method (Clausen 2022)
- [ ] Bipartitions, distances, dissimilarities &c.
  - [ ] Node distance
  - [ ] Patristic distance
  - [ ] Reduce memory use of `BipartitionTable`
  - [ ] `is_refinement` method
  - [ ] Robinson-Foulds (Penny & Hendy 1985)
    - [x] Basic implementation
    - [ ] Normalisations
    - [ ] Weighted variants
  - [ ] Clustering information distance (Smith 2020)
    - [x] Basic implementation
    - [ ] Associated methods (common clustering information, entropy, &c.)
    - [ ] Normalisations
    - [ ] Reduce LAP allocations
    - [ ] Weighted variants?
  - [ ] Make distance functions follow the Distances package
    - [ ] Implement distance type hierarchy
    - [ ] Unify all methods under a `dist` function (for both tree and node distances)
  - [ ] Tree dissimilarities / differences (semi-metric)
    - [ ] Contradiction dissimilarity (Bapst et a. 2018)
    - [ ] CI contradiction dissimilarity?
- [ ] Implement tree consensus functions
  - [ ] Majority-rule (Jensen et al. 2016)
  - [ ] Strict (Day 1985?)
- [ ] Implement branch swapping functions
  - [x] Swap nodes
  - [ ] NNI
  - [ ] SPR
  - [ ] TBR
- [ ] Newick writing and reading
  - [x] Write (`Newick.write`, `Newick.string`; don't implement print methods)
  - [x] Read (`Newick.parse`, `Newick.read`)
  - [ ] Read metacomment data
    - [ ] Custom readers (RAxML-NG, BEAST2, MrBayes, etc.)
    - [ ] Write metacomment data
- [ ] Implement parsimony
  - [ ] Data structure for character/sequence data
  - [ ] Unweighted Fitch kernel (Goloboff 2022; White & Holland 2011)
  - [ ] Weighted Fitch kernel
  - [ ] Sankoff kernel
  - [ ] Tree scores, ASR
  - [ ] Rearrangement scores
  - [ ] Search
  - [ ] Dynamic kernels

## Design

- [ ] Implement taxon deletion from `TaxonSet`?
- [ ] "Edge-matrix" (APE-style) tree data structure
- [ ] Neighbour-vector (PHYML/GoTree-style) tree data structure
- [ ] Rooted-tree specialisations
  - [ ] Time tree
