"""
    TaxonBipartition

Bit-vector that tags taxa in a list as belonging to one of two disjoint groups.

Beware that many functions assume that bipartitions are normalised so that the first bit is
always naught (see Penny & Hendy 1985).
"""
struct TaxonBipartition <: AbstractArray{Bool,1}
    v::BitVector

    TaxonBipartition(v::BitVector) = new(v)
    TaxonBipartition(ntax) = new(falses(ntax))
end

Base.size(bp::TaxonBipartition) = (length(bp.v), )
Base.getindex(bp::TaxonBipartition, i::Int) = getindex(bp.v, i)
Base.IndexStyle(::Type{<:TaxonBipartition}) = IndexLinear()
Base.firstindex(bp::TaxonBipartition) = firstindex(bp.v)
Base.lastindex(bp::TaxonBipartition) = lastindex(bp.v)
Base.setindex!(bp::TaxonBipartition, v, i::Int) = setindex!(bp.v, v, i)
Base.similar(bp::TaxonBipartition) = TaxonBipartition(falses(length(bp)))

Base.:~(bp::TaxonBipartition) = TaxonBipartition(.~ bp.v)
Base.:|(bp1::TaxonBipartition, bp2::TaxonBipartition) = TaxonBipartition(bp1.v .| bp2.v)
Base.:&(bp1::TaxonBipartition, bp2::TaxonBipartition) = TaxonBipartition(bp1.v .& bp2.v)

function Base.isequal(a::TaxonBipartition, b::TaxonBipartition)
    for (chunk1, chunk2) in zip(a.v.chunks, b.v.chunks)
        chunk1 == chunk2 || return false
    end
    
    return true
end

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


"""
    normalise!(bipartition)

Normalise a bipartition so that the first bit is always naught (Penny & Hendy 1985)
"""
function normalise!(bp)
    if first(bp.v)
        bp.v .= .~ bp.v
    end
end



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
    s = 0
    for chunk in bp.v.chunks
        s += count_ones(chunk)
    end

    return s == 1 || s == length(bp) - 1
end


"""
    are_compatible(bipartition1, bipartition2)

Check compatibility between two partitions.

Two bipartitions B1 and B2 are compatible if it is logically possible for them to exist in 
the same tree.

This can be verified with any one of the following conditions being true:
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
    "Map from bp to bp index in the table"
    dict::OrderedDict{TaxonBipartition,UInt}
    "Table with the indices of bps in each tree, column-wise"
    records::Matrix{UInt}
    "Number of bps encountered in the table"
    nbp::UInt

    function BipartitionTable(ntax, ntrees)
        bps_in_tree = ntax - 3
        dict = OrderedDict{TaxonBipartition,UInt}()
        records = zeros(UInt, bps_in_tree, ntrees)

        return new(dict, records, 0)
    end
end


Base.show(io::IO, bpt::BipartitionTable) =
    print(io, "BipartitionTable: $(bpt.nbp) bipartitions from $(size(bpt.records, 2)) trees")


function Base.getproperty(bpt::BipartitionTable, key::Symbol)
    if key == :bipartitions
        return getfield(getfield(bpt, :dict), :keys)
    else
        return getfield(bpt, key)
    end
end


function add_record!(bp_table::BipartitionTable, bp, treecol, treerow)
    normalise!(bp)
    bp_index = get(bp_table.dict, bp, 0)

    if bp_index == 0
        bp_table.nbp += 1
        bp_index = bp_table.nbp
        bp_table.dict[bp] = bp_index
    end

    @inbounds bp_table.records[treerow, treecol] = bp_index

    return treerow + 1
end


function compute_bipartition!(
    bp_table::BipartitionTable,
    branchnode,
    ntax,
    treecol,
    treerow
    )
    _, node = branchnode
    bp = TaxonBipartition(ntax)

    for cbranchnode in children(branchnode)
        _, cnode = cbranchnode
        if isouter(cnode)
            cnode.taxon > 0 || continue
            @inbounds bp.v[cnode.taxon] = true
        else
            treerow, cbp =
                compute_bipartition!(bp_table, cbranchnode, ntax, treecol, treerow)
            bp.v.chunks .|= cbp.v.chunks
            treerow = add_record!(bp_table, cbp, treecol, treerow)
        end
    end

    if node.taxon > 0
        @inbounds bp.v[node.taxon] = true
    end

    return treerow, bp
end


"""
    bipartition_table(trees, trim=true)

Return a record of the incidences of bipartitions in a collection of trees.

The table contains:

- A list with all the unique bipartitions encountered in collection of trees
- A record of the bipartitions that occur in each tree, identified by their index in the list.

The record consists of a matrix where the trees are arranged by columns. Trees that are not
fully resolved will contain zeros representing the absence of bipartitions.

# Properties
(Accessable through the dot-syntax)

- `bipartitions`: The list of unique bipartitions
"""
function bipartition_table(trees)
    NType = nodetype(eltype(trees))
    ntrees = length(trees)
    ntax = length(first(trees).taxonset)
    bp_table = BipartitionTable(ntax, ntrees)

    progmeter = Progress(
        ntrees;
        output=stderr,
        dt=2,
        desc="Enumerating bipartitions",
        showspeed=true,
        enabled=!is_logging(stderr)
        )

    for (treecol, tree) in enumerate(trees)
        compute_root::NType = tree.root
        if isouter(compute_root)
            compute_root = last(only(children(compute_root)))
        end
        
        treerow = 1
        for cbranchnode in children(nothing => compute_root)
            _, cnode::NType = cbranchnode
            isouter(cnode) && continue
            treerow, bp =
                compute_bipartition!(bp_table, cbranchnode, ntax, treecol, treerow)
            treerow = add_record!(bp_table, bp, treecol, treerow)
        end

        next!(progmeter)
    end

    @showprogress desc="Sorting bipartition indices" dt=2 showspeed=true foreach(
        sort!, eachcol(bp_table.records)
        )
    
    return bp_table
end


function compute_bipartition!(bplist::Vector{TaxonBipartition}, branchnode, i)
    _, node = branchnode
    @inbounds bp = bplist[i]

    for cbranchnode in children(branchnode)
        _, cnode = cbranchnode
        if isouter(cnode)
            cnode.taxon > 0 || continue
            @inbounds bp.v[cnode.taxon] = true
        else
            i, cbp = compute_bipartition!(bplist, cbranchnode, i + 1)
            bp.v.chunks .|= cbp.v.chunks
            normalise!(cbp)
        end
    end

    if node.taxon > 0
        @inbounds bp.v[node.taxon] = true
    end

    return i, bp
end


"""
    bipartitions(tree)

Return a vector of `TaxonBipartition`s that represent the bipartitions in the tree.

Note that inner nodes of degree 2 (e.g. a bifurcation at the root) might cause some 
bipartitions to be duplicated.
"""
function bipartitions(tree::T) where T <: AbstractTree
    NType = nodetype(T)
    ntax = length(tree.taxonset)
    n_all_bp = length(tree.branches)
    n_inner_bp = n_all_bp - ntax
    bplist = [TaxonBipartition(ntax) for _ in 1:n_inner_bp]

    i = 1
    compute_root::NType = tree.root
    if isouter(compute_root)
        compute_root = last(only(children(compute_root)))
    end

    for cbranchnode in children(nothing => compute_root)
        _, cnode::NType = cbranchnode
        isouter(cnode) && continue
        i, bp = compute_bipartition!(bplist, cbranchnode, i)
        normalise!(bp)
        i += 1
    end

    return bplist
end
