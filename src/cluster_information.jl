"""
Return H * n
"""
function _clustering_entropy(bp::TaxonBipartition, n=length(bp))
    A = 0
    for chunk in bp.v.chunks    # Faster than using `sum(bp)`
        A += count_ones(chunk)
    end
    B = n - A
    
    return - A * log2(A / n) - B * log2(B / n)
end

"""
    clustering_entropy(bipartition, n=length(bipartition))

Return the entropy associated to a taxon bipartition as given by the formula from Smith 2020:

``
H = (-A/n) * log₂(A/n) - (B/n) * log₂(B/n)
``

Where A and B represent the sizes of each bipartition and n represents the total number of taxa.
"""
clustering_entropy(bp::TaxonBipartition, n=length(bp)) = _clustering_entropy(bp, n) / n


function _ic_element(x, y, ixy, n)
    ixy > 0 && x > 0 && y > 0 || return 0.0

    return ixy * (log2(n * ixy) - log2(x * y))
end


"""
Return s * n
"""
function _mutual_clustering_information(bp1::TaxonBipartition, bp2::TaxonBipartition, n=length(bp1))
    A1, A2 = sum(bp1), sum(bp2)
    B1, B2 = n - A1, n - A2

    # Hacky bit: Directly accessing BitVector chunks 
    iA1A2 = 0
    for (x1, x2) in zip(bp1.v.chunks, bp2.v.chunks)
        iA1A2 += count_ones(x1 & x2)
    end
    iA1B2 = A1 - iA1A2
    iB1A2 = A2 - iA1A2
    iB1B2 = B1 - iB1A2

    s = _ic_element(A1, A2, iA1A2, n) +
        _ic_element(A1, B2, iA1B2, n) +
        _ic_element(B1, A2, iB1A2, n) +
        _ic_element(B1, B2, iB1B2, n)

    return s
end


"""
    mutual_clustering_information(bipartition1, bipartition2)

Return the mutual clustering information between two bipartitions.
"""
function mutual_clustering_information(bp1::TaxonBipartition, bp2::TaxonBipartition)
    n = length(bp1)

    return _mutual_clustering_information(bp1, bp2, n) / n
end


"""
Return matrix with s * n value of each bipartition matching.
"""
function mci_cost_matrix(bplist1, bplist2, n)
    lrows = length(bplist1)
    lcols = length(bplist2)
    M = zeros(lrows, lcols)

    for j in 1:lcols
        for i in 1:lrows
            s = _mutual_clustering_information(bplist1[i], bplist2[j], n)
            @inbounds M[i, j] = s
        end
    end

    return M
end

function clustering_information(tree1::AbstractTree, tree2::AbstractTree; singletons=false)
    bplist1 = bipartitions(tree1, singletons=singletons)
    bplist2 = bipartitions(tree2, singletons=singletons)

    return clustering_information(bplist1, bplist2)
end

function clustering_information(bplist1::Vector{TaxonBipartition}, bplist2::Vector{TaxonBipartition}, n = length(bplist1[1]))
    M = mci_cost_matrix(bplist1, bplist2, n)
    optimal_matching = find_best_assignment(M, true)
    
    return optimal_matching.cost / n
end


function clustering_information_dist(
    bplist1::Vector{TaxonBipartition},
    bplist2::Vector{TaxonBipartition},
    k = length(bplist1[1]))
    
    empty_bplist1 = isempty(bplist1)
    empty_bplist2 = isempty(bplist2)

    #TODO: Create a specialised method with the assumption that trees are fully resolved.
    # With that assumption we can eliminate two branch instructions
    if empty_bplist1
        if empty_bplist2
            return 0.0
        else
            return 0.5 * sum(_clustering_entropy(bp) for bp in bplist2) / k
        end
    end

    isempty(bplist2) &&
        return 0.5 * sum(_clustering_entropy(bp) for bp in bplist1) / k

    s = clustering_information(bplist1, bplist2, k)
    H = sum(_clustering_entropy(bp) for bp in bplist1) +
        sum(_clustering_entropy(bp) for bp in bplist2)

    return 0.5 * H / k - s
end


function clustering_information_dist(trees; diffonly=true, singletons=false)
    n = length(trees)
    k = length(first(trees).taxonset)
    bptable = bipartition_table(trees; singletons=singletons)
    bpincidence1 = falses(bptable.m)
    bpincidence2 = falses(bptable.m)
    
    d = zeros(Float64, n, n)
    
    if diffonly
        bpincidence_different = falses(bptable.m)

        for j in 1:n, i in 1:(j - 1)
            # Select only the bipartitions that are not shared between the two trees
            bpincidence1 .= bptable.data[1:bptable.m, i]
            bpincidence2 .= bptable.data[1:bptable.m, j]
            bpincidence_different .= bpincidence1 .⊻ bpincidence2
            bplist1 = bptable.bipartitions[bpincidence1 .& bpincidence_different]
            bplist2 = bptable.bipartitions[bpincidence2 .& bpincidence_different]
    
            @inbounds d[i, j] = clustering_information_dist(bplist1, bplist2, k)
        end
    else
        for j in 1:n, i in 1:(j - 1)
            bpincidence1 .= bptable.data[:, i]
            bpincidence2 .= bptable.data[:, j]
            bplist1 = bptable.bipartitions[bpincidence1]
            bplist2 = bptable.bipartitions[bpincidence2]
    
            @inbounds d[i, j] = clustering_information_dist(bplist1, bplist2, k)
        end
    end

    return d
end


function clustering_information_dist(
    tree1::AbstractTree,
    tree2::AbstractTree;
    diffonly=true,
    singletons=false
    )
    k = length(tree1.taxonset)
    bptable = bipartition_table((tree1, tree2); singletons=singletons)

    bpincidence1 = bptable.data[:, 1]
    bpincidence2 = bptable.data[:, 2]

    if diffonly
        # Select only the bipartitions that are not shared between the two trees
        bpincidence_different = bpincidence1 .⊻ bpincidence2
        bpincidence1 .&= bpincidence_different
        bpincidence2 .&= bpincidence_different
    end

    bplist1 = bptable.bipartitions[bpincidence1]
    bplist2 = bptable.bipartitions[bpincidence2]

    d = clustering_information_dist(bplist1, bplist2, k)
        
    return d
end


##### Alternative


function clustering_information_dist2(trees; diffonly=true, singletons=false)
    n = length(trees)
    k = length(first(trees).taxonset)
    bptable = bipartition_table2(trees; singletons=singletons)
    d = zeros(Float64, n, n)
    bpinds1 = zeros(Int, bptable.m)
    bpinds2 = zeros(Int, bptable.m)
    
    if diffonly
        for j in 1:n, i in 1:(j - 1)
            # Select only the bipartitions that are not shared between the two trees
            m1, m2 = two_way_setdiff(bptable.data[i], bptable.data[j], bpinds1, bpinds2)
            bplist1, bplist2 = bptable.dict.keys[bpinds1[1:m1]], bptable.dict.keys[bpinds2[1:m2]]
            # bpincidence1 .= bptable.data[1:bptable.m, i]
            # bpincidence2 .= bptable.data[1:bptable.m, j]
            # bpincidence_different .= bpincidence1 .⊻ bpincidence2
            # bplist1 = bptable.bipartitions[bpincidence1 .& bpincidence_different]
            # bplist2 = bptable.bipartitions[bpincidence2 .& bpincidence_different]
            
            @inbounds d[i, j] = clustering_information_dist(bplist1, bplist2, k)
            bpinds1 .= bpinds2 .= 0
        end
    else
        for j in 1:n, i in 1:(j - 1)
            bpinds1 = bptable.data[i]
            bpinds2 = bptable.data[j]
            bplist1 = bptable.bipartitions[bpinds1]
            bplist2 = bptable.bipartitions[bpinds2]
    
            @inbounds d[i, j] = clustering_information_dist(bplist1, bplist2, k)
        end
    end

    return d
end


