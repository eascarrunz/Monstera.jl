function rf_dist(trees; singletons = false)
    ntree = length(trees)
    d = zeros(Int, ntree, ntree)
    bptable = bipartition_table(trees, false; singletons = singletons)
    bpincidence_buffer = falses(bptable.m)
    for j in 1:ntree, i in 1:j
        @inbounds bpincidence_buffer .= bptable.data[1:bptable.m, i]
        @inbounds bpincidence_buffer .⊻= bptable.data[1:bptable.m, j]
        @inbounds d[i, j] = sum(bpincidence_buffer)
        bpincidence_buffer .= false
    end

    return d
end

function rf_dist2(trees; singletons = false)
    ntree = length(trees)
    d = zeros(Int, ntree, ntree)
    bptable = bipartition_table2(trees; singletons = singletons)
    for j in 1:ntree, i in 1:j
        @inbounds d[i, j] = count_symdiff(bptable.data[i], bptable.data[j])
    end

    return d
end

function rf_dist(tree1::AbstractTree, tree2::AbstractTree; singletons = false)
    bptable = bipartition_table((tree1, tree2); singletons = singletons)

    return @inbounds sum(bptable.data[:, 1] .⊻ bptable.data[:, 2])
end


"""
    distance(trees, mode = :ci; <keyword arguments>)
    distance(tree1, tree2, mode = :ci; <keyword arguments>)

Return the distance between two trees or a the upper triangle of the matrix of pairwise 
distances among trees in a collection.

Two distances are available through the `mode` argument: **Clustering Information Distance** 
(`:ci`) (Smith 2020) and **Robinson-Foulds** (`:rf`) (Robinson & Foulds 1981).

It is assumed that all trees belong to a same taxon set and have the same taxa among them, 
that the trees are unrooted, and that polytomies are hard.

# Keyword arguments

- `singletons=false`: if `true`, singleton bipartitions (corresponding to outer branches) \
are considered when computing the distances. This deviates from the original definitions of\
 the distances, but might be useful for trees with taxa in inner nodes ("sampled \
ancestors").

## For Clustering Information Distance (`:ci` mode)

- `diffonly=true`: if `true`, the distance is computed only from the bipartitions that are \
not shared between the two trees. This makes small distances more accurate and faster to \
compute.
""" 
function distance(trees, mode = :ci; kwargs...)
    if mode == :ci
        return clustering_information_dist(trees; kwargs...)
    elseif mode == :rf
        return rf_dist(trees)
    end
end

function distance(tree1::AbstractTree, tree2::AbstractTree, mode = :ci; kwargs...)
    if mode == :ci
        return clustering_information_dist(tree1, tree2; kwargs...)
    elseif mode == :rf
        return rf_dist(tree1, tree2)
    end
end

