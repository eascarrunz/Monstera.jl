"Bipartition table grows by adding its current number of entries times 
`BPTABLE_GROWTH_FACTOR`"
const BPTABLE_GROWTH_FACTOR = 1

"""
    TaxonBipartition

Bit-vector that tags taxa in a list as belonging to one of two disjoint groups.

Beware that many functions assume that bipartitions are normalised so that the first bit is
always naught (see Penny & Hendy 1985).
"""
struct TaxonBipartition <: AbstractArray{Bool,1}
    v::BitVector
end

Base.size(bp::TaxonBipartition) = (length(bp.v), )
Base.getindex(bp::TaxonBipartition, i::Int) = getindex(bp.v, i)
Base.IndexStyle(::Type{<:TaxonBipartition}) = IndexLinear()
Base.setindex!(bp::TaxonBipartition, v, i::Int) = setindex!(bp.v, v, i)
Base.length(bp::TaxonBipartition) = length(bp.v)
Base.similar(bp::TaxonBipartition) = TaxonBipartition(falses(length(bp)))
Base.firstindex(::TaxonBipartition) = firstindex(bp.v)
Base.lastindex(bp::TaxonBipartition) = lastindex(bp.v)

Base.:~(bp::TaxonBipartition) = TaxonBipartition(.~ bp.v.chunks)
Base.:|(bp1::TaxonBipartition, bp2::TaxonBipartition) =
    TaxonBipartition(bp1.v.chunks .| bp2.v.chunks)
Base.:&(bp1::TaxonBipartition, bp2::TaxonBipartition) =
    TaxonBipartition(bp1.v.chunks .& bp2.v.chunks)

function Base.count_ones(bp::TaxonBipartition)
    s = 0
    for chunk in bp.v.chunks
        s += count_ones(chunk)
    end

    return s
end

function Base.count_zeros(bp::TaxonBipartition)
    s = 0
    for chunk in bp.v.chunks
        s += count_zeros(chunk)
    end
    
    return s
end


"""
Print taxon bipartition grouping taxa by 5
"""
function Base.print(io::IO, bp::TaxonBipartition)
    l = lastindex(bp)
    for (i, taxon_bit) in enumerate(bp)
        print(io, taxon_bit ? '■' : "□")
        (i % 5 == 0 && i ≠ l) && print(io, ' ')
    end

    return nothing
end

Base.show(io::IO, bp::TaxonBipartition) = print(io, string(bp))

function Base.show(io::IO, ::MIME"text/plain", bp::TaxonBipartition)
    print(io, length(bp), "-taxon TaxonBipartition:\n ")
    show(io, bp)

    return nothing
end


function compute_bipartition(branchnode, bpvec, k)
    branch, node = branchnode
    bp = bpvec[branch.id]

    if node.taxon > 0
        bp.v[node.taxon] = true
    end

    for cbranchnode in children(branchnode)
        bp.v.chunks .|= compute_bipartition(cbranchnode, bpvec, k).v.chunks
    end

    return bp
end


"""
    normalise!(bipartition)

Normalise a bipartition so that the first bit is always naught (Penny & Hendy 1985)
"""
function normalise!(bp)
    if first(bp)
        bp.v.chunks .= .~ bp.v.chunks
    end
end


struct TaxonBipartitionList
    bpvec::Vector{TaxonBipartition}
    dict::OrderedDict{Int,Int}
    nbp::Int
end


Base.getindex(x::TaxonBipartitionList, inds...) = getindex(x.bpvec, inds...)
Base.length(x::TaxonBipartitionList) = x.nbp
Base.size(x::TaxonBipartitionList) = (x.nbp, )
Iterators.eltype(::TaxonBipartitionList) = TaxonBipartition
Base.first(x::TaxonBipartitionList) = x.bpvec[begin]
Base.last(x::TaxonBipartitionList) = x.bpvec[end]
Base.firstindex(x::TaxonBipartitionList) = firstindex(x.bpvec)
Base.lastindex(x::TaxonBipartitionList) = lastindex(x.bpvec)
Base.iterate(iter::TaxonBipartitionList, state=1) =
    state > iter.nbp ? nothing : (iter.bpvec[state], state + 1)


struct TaxonBipartitionListBranchIndexer
    bplist::TaxonBipartitionList
end


function Base.getproperty(bplist::TaxonBipartitionList, key::Symbol)
    if key == :bybranch
        return TaxonBipartitionListBranchIndexer(bplist)
    else
        getfield(bplist, key)
    end
end


Base.getindex(x::TaxonBipartitionListBranchIndexer, ind...) =
    x.bplist.bpvec[x.bplist.dict[ind...]]
Base.firstindex(x::TaxonBipartitionListBranchIndexer) = 1
Base.lastindex(x::TaxonBipartitionListBranchIndexer) = x.bplist.nbp


"""
    is_singleton(bipartition)

Return true is the bipartition contains a singleton, i.e. it separates a single taxon from
all the rest.

In other words, a singleton represents a taxon at a leaf node of a tree. This kind of
bipartition is commonly called "trivial", because they provide no additional information
about tree topology if one knows that every taxon is a leaf node. But singleton bipartitions
    are *not* "trivial" in trees with sampled ancestors.
"""
function is_singleton(bp::TaxonBipartition)
    s = count_ones(bp)

    return s == 1 || s == length(bp) - 1
end


"""
    raw_bipartitions(tree)

Return a the raw (=unnormalised) bipartitions of a tree, singletons included.
"""
function raw_bipartitions(tree::AbstractTree)
    ntax::Int=length(tree.taxonset)
    bpvec = [TaxonBipartition(falses(ntax)) for _ in 1:length(tree.branches)]

    for cbranchnode in children(nothing => tree.root)
        compute_bipartition(cbranchnode, bpvec, ntax)
    end

    return bpvec
end


"""
    bipartitions(tree, singletons=true)

Return a list of taxon bipartitions corresponding to the branches of a tree.

The tree and its nodes must be associated to a taxon set. Singleton bipartitions are 
excluded by default.

The returned list is of type `TaxonBipartitionList`. It is possible to retrieve the 
bipartition corresponding to a branch with the `bybranch` property.
"""
function bipartitions(tree::AbstractTree; singletons=false)
    bpvec = raw_bipartitions(tree)
    dict = OrderedDict{Int,Int}()

    if ! singletons
        n_kept_bp = 0
        keep_mask = falses(length(bpvec))
        for (i, bp) in enumerate(bpvec)
            this_is_singleton = is_singleton(bp)
            keep_mask[i] = ! this_is_singleton
            if ! this_is_singleton
                n_kept_bp += 1
                dict[i] = n_kept_bp
            end
        end

        bpvec = bpvec[keep_mask]
    else
        n_kept_bp = length(bpvec)
    end

    for bp in bpvec
        normalise!(bp)
    end
    
    return TaxonBipartitionList(bpvec, dict, n_kept_bp)
end


"""
    are_compatible(bipartition1, bipartition2)

Check compatibility between two partitions.

Two bipartitions B1 and B2 are compatible if it is logically possible for them to exist in 
the same tree.

This can be verified with one of the following conditions being true:
- B1 = B1 ∪ B2
- B1 = B1 ∪ ¬B2
- ¬B1 = ¬B1 ∪ B2
- ¬B1 = ¬B1 ∪ ¬B2
"""
function are_compatible(bp1::TaxonBipartition, bp2::TaxonBipartition)
    for (chunk1, chunk2) in zip(bp1.v.chunks, bp2.v.chunks)
        (chunk1 == chunk1 | chunk2) ||
        (chunk1 == chunk1 | ~chunk2) ||
        (~chunk1 == ~chunk1 | chunk2) ||
        (~chunk1 == ~chunk1 | ~chunk2) ||
        return false
    end

    return true
end


mutable struct BipartitionTable
    "Mapping of bipartitions to rows in the table"
    dict::OrderedDict{TaxonBipartition,Int}
    "Table of occurrence of bipartitions (rows = bipartitions, columns = trees)"
    data::BitMatrix
    "Number of trees"
    n::Int
    "Index of the last unique bipartition added to the table"
    m::Int
    "Current capacity for unique bipartitions"
    cap::Int

    """
        BipartitionTable(n, k)

    Initalise a taxon bipartition incidence table with the following parameters for the size
    hint:
    - `n`: Number of trees
    - `k`: Number of taxa
    """
    function BipartitionTable(n, k)
        # Expected number of unique bipartitions in the tree set
        cap = (k - 3) * n ÷ 5

        data = falses(cap, n)

        return new(OrderedDict{TaxonBipartition,Int}(), data, n, 0, cap)
    end
end


Base.show(io::IO, bpt::BipartitionTable) =
    print(io, "BipartitionTable: $(bpt.m) bipartitions from $(bpt.n) trees")

Base.getindex(bpt::BipartitionTable, bp::TaxonBipartition) =
    getfield(bpt, :data)[getfield(bpt, :dict)[bp], :]

Base.getindex(bpt::BipartitionTable, i::Int) = getfield(bpt, :data)[1:getfield(bpt, :m), i]

Base.size(bpt::BipartitionTable) = bpt.m, bpt.n


function Base.getproperty(bpt::BipartitionTable, key::Symbol)
    if key == :bipartitions
        return getfield(getfield(bpt, :dict), :keys)
    elseif key == :frequencies
        return sum(getfield(bpt, :data), dims=2)
    else
        return getfield(bpt, key)
    end
end


function record_incidence!(bptable::BipartitionTable, bp::TaxonBipartition, i)
    bpindex = get(bptable.dict, bp, 0)
    if bpindex == 0
        bptable.m += 1
        bpindex = bptable.m
        bptable.dict[bp] = bpindex
        if bptable.m > bptable.cap
            new_entry = falses(BPTABLE_GROWTH_FACTOR * bptable.m, bptable.n)
            bptable.cap += BPTABLE_GROWTH_FACTOR * bptable.m
            bptable.data = vcat(bptable.data, new_entry)
        end
    end

    @inbounds bptable.data[bpindex, i] = true
end


"""
    bipartition_table(trees, trim=true)

Return a record of the incidences of bipartitions in a collection of trees.

The bipartition table contains a matrix where the rows represent unique bipartitions and the
columns represent trees. The entries of the matrix are `true` when a bipartition occurs in
a tree. The matrix can be accessed by indexing the table either by bipartition or by tree 
index (multiple indices and ranges are not supported yet).

A number of entries are pre-allocated based on a guess of how many unique biparitions may be
in the collection. There are removed by default (`trim=true`).

# Properties
(Accessable through the dot-syntax)

- `bipartitions`: The list of unique bipartitions
- `frequencies`: The number of times each bipartition occurs in a tree
- `n`: The number of trees
- `m`: The number of unique bipartitions
"""
function bipartition_table(trees, trim=true; singletons = false)
    n = length(trees)
    k = length(first(trees).taxonset)
    bptable = BipartitionTable(n, k)
    skip_singletons = ! singletons

    for (i, tree) in enumerate(trees)
        bpvec = raw_bipartitions(tree)
        for bp in bpvec
            is_singleton(bp) & skip_singletons && continue
            normalise!(bp)
            record_incidence!(bptable, bp, i)
        end
    end

    if trim
        bptable.data = bptable.data[1:bptable.m, :]
    end

    return bptable
end


#### Alternative implementations


mutable struct BipartitionTable2
    "Mapping of bipartitions to rows in the table"
    dict::OrderedDict{TaxonBipartition,Int}
    "Table of occurrence of bipartitions (rows = bipartitions, columns = trees)"
    data::Vector{Vector{UInt32}}
    "Number of trees"
    n::Int
    "Index of the last unique bipartition added to the table"
    m::Int

    function BipartitionTable2(n, k, singletons=false)
        npb = singletons ? (2 * k - 3) : (k - 3)
        data = [zeros(UInt32, npb) for _ in 1:n]

        return new(OrderedDict{TaxonBipartition,Int}(), data, n, 0)
    end
end


function bipartition_table2(trees; singletons=false)
    ntree = length(trees)
    k = length(first(trees).taxonset)
    bptable = BipartitionTable2(ntree, k, singletons)
    skip_singletons = ! singletons

    for (i, tree) in enumerate(trees)
        bpvec = raw_bipartitions(tree)
        j = 0
        for bp in bpvec
            is_singleton(bp) & skip_singletons && continue
            j += 1
            normalise!(bp)
            record_incidence!(bptable, bp, i, j)
        end
        sort!(bptable.data[i], rev=true)
    end

    return bptable
end




function record_incidence!(bptable::BipartitionTable2, bp::TaxonBipartition, i, j)
    bpindex = get(bptable.dict, bp, 0)
    if bpindex == 0
        bptable.m += 1
        bpindex = bptable.m
        bptable.dict[bp] = bpindex
    end
    
    bptable.data[i][j] = bpindex

    return nothing
end
    
