function rf_dist(trees; progress::Bool=true)
    ntree = length(trees)
    D = zeros(Int, ntree, ntree)
    bptable = bipartition_table(trees)
    Dsize = ntree * (ntree - 1) ÷ 2 
    if progress
        progmeter = Progress(
            Dsize;
            output = stderr,
            dt=2,
            desc="Computing pairwise distances",
            showspeed=true,
            enabled=! is_logging(stderr)
            )
    end
    for j in 1:ntree, i in 1:j
        @inbounds D[i, j] = special_countdiff(bptable.records[:, i], bptable.records[:, j])
        progress && next!(progmeter)
    end
    progress && finish!(progmeter)

    return D
end

function rf_dist(tree1::AbstractTree, tree2::AbstractTree)
    bptable = bipartition_table((tree1, tree2))

    return @inbounds special_countdiff(bptable.records[:, 1], bptable.records[:, 2])
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

## Keyword arguments

- `diffonly=true`: if `true`, only bipartitions that are not shared between the two trees 
  will be considered. Only applicable to the `:ci` distance.
- `progress=true`: whether to print a progress bar and ETA for long runs. Disable for Pluto
compatibility.
""" 
function distance(trees, mode = :ci; kw...)
    if mode == :ci
        return clustering_information_dist(trees; kw...)
    elseif mode == :rf
        return rf_dist(trees; kw...)
    end
end

function distance(tree1::AbstractTree, tree2::AbstractTree, mode = :ci; kw...)
    if mode == :ci
        return clustering_information_dist(tree1, tree2; kw...)
    elseif mode == :rf
        return rf_dist(tree1, tree2; kw...)
    end
end

