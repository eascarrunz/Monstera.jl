"""
Return H * ntax
"""
function _clustering_entropy(bp::TaxonBipartition, ntax=length(bp))
    A = sum(bp)
    B = ntax - A
    
    return - A * log2(A / ntax) - B * log2(B / ntax)
end

"""
    clustering_entropy(bipartition, n=length(bipartition))
    clustering_entropy(bipartitions)

Return the entropy associated to a taxon bipartition as given by the formula from Smith 2020:

``
H = (-A/n) * log₂(A/n) - (B/n) * log₂(B/n)
``

Where A and B represent the sizes of each bipartition and n represents the total number of taxa.
"""
clustering_entropy(bp::TaxonBipartition, ntax=length(bp)) =
    _clustering_entropy(bp, ntax) / ntax
clustering_entropy(bplist) = sum(clustering_entropy.(bplist))
clustering_entropy(tree::AbstractTree) = clustering_entropy(bipartitions(tree))


# function _ic_element(x, y, ixy, ntax)
#     ixy > 0 && x > 0 && y > 0 || return 0.0

#     return ixy * (log2(ntax * ixy) - log2(x * y))
# end


# """
# Return s * ntax
# """
# function _mutual_clustering_information(
#     bp1::TaxonBipartition,
#     bp2::TaxonBipartition,
#     ntax=length(bp1)
#     )
#     A1, A2 = sum(bp1), sum(bp2)
#     B1, B2 = ntax - A1, ntax - A2

#     iA1A2 = 0
#     for (x1, x2) in zip(bp1.v.chunks, bp2.v.chunks)
#         iA1A2 += count_ones(x1 & x2)
#     end
#     iA1B2 = A1 - iA1A2
#     iB1A2 = A2 - iA1A2
#     iB1B2 = B1 - iB1A2

#     s = _ic_element(A1, A2, iA1A2, ntax) +
#         _ic_element(A1, B2, iA1B2, ntax) +
#         _ic_element(B1, A2, iB1A2, ntax) +
#         _ic_element(B1, B2, iB1B2, ntax)

#     return s
# end


"""
    mutual_clustering_information(bipartition1, bipartition2)

Return the mutual clustering information between two bipartitions.
"""
function mutual_clustering_information(bp1::TaxonBipartition, bp2::TaxonBipartition)
    n = length(bp1)

    return _mutual_clustering_information(bp1, bp2, n) / n
end


function mutual_clustering_information(tree1::AbstractTree, tree2::AbstractTree)
    bplist1 = unique(bipartitions(tree1))
    bplist2 = unique(bipartitions(tree2))

    return mutual_clustering_information(bplist1, bplist2)
end

function _mutual_clustering_information!(M::AbstractMatrix{Float64}, bplist1, bplist2, lrows, lcols, ntax)
    mci_cost_matrix!(M, bplist1, bplist2, lrows, lcols, ntax)
    optimal_matching = find_best_assignment(@view(M[1:lrows, 1:lcols]), true)
    
    return optimal_matching.cost / ntax
end


function mutual_clustering_information(bplist1, bplist2, ntax = length(bplist1[1]))
    lrows = length(bplist1)
    lcols = length(bplist2)
    M = zeros(lrows, lcols)
    
    return _mutual_clustering_information!(M, bplist1, bplist2, lrows, lcols, ntax)
end


function clustering_information_dist!(
    M::AbstractMatrix{Float64},
    bplist1,
    bplist2,
    ntax = length(first(bplist1))
    )
    
    lrows = length(bplist1)
    lcols = length(bplist2)

    if lrows == 0
        if lcols == 0
            return 0.0
        else
            return 0.5 * sum(_clustering_entropy(bp) for bp in bplist2) / ntax
        end
    end

    lcols == 0 &&
        return 0.5 * sum(_clustering_entropy(bp) for bp in bplist1) / ntax

    s = _mutual_clustering_information!(M, bplist1, bplist2, lrows, lcols, ntax)
    H = sum(_clustering_entropy(bp) for bp in bplist1) +
        sum(_clustering_entropy(bp) for bp in bplist2)

    return 0.5 * H / ntax - s
end


function clustering_information_dist(trees; diffonly=true)
    ntrees = length(trees)
    ntax = length(first(trees).taxonset)

    bptable = bipartition_table(trees)

    bpindices1 = zeros(UInt32, size(bptable.records, 1))
    bpindices2 = zeros(UInt32, size(bptable.records, 1))

    M = zeros(ntax - 3, ntax - 3)          # Cost matrix for bipartition matching
    D = zeros(Float64, ntrees, ntrees)     # Distance matrix
    Dsize = ntrees * (ntrees - 1) ÷ 2      # Number of actual non-zero entries
    progmeter = Progress(Dsize;
        output=stderr,
        dt=2,
        desc="Computing pairwise distances",
        showspeed=true,
        enabled=!is_logging(stderr)
        )
    
    if diffonly
        # Select only the bipartitions that are not shared between the two trees
        unique_indices1 = zeros(UInt32, size(bptable.records, 1))
        unique_indices2 = zeros(UInt32, size(bptable.records, 1))

        @inbounds for j in 2:ntrees, i in 1:(j - 1)
            bpindices1 .= bptable.records[:, i]
            bpindices2 .= bptable.records[:, j]
            nbp1, nbp2 =
                special_setdiff!(unique_indices1, unique_indices2, bpindices1, bpindices2)
            bplist1 = bptable.bipartitions[unique_indices1[1:nbp1]]
            bplist2 = bptable.bipartitions[unique_indices2[1:nbp2]]
    
            D[i, j] = clustering_information_dist!(M, bplist1, bplist2, ntax)

            unique_indices1 .= zero(UInt32)
            unique_indices2 .= zero(UInt32)
            next!(progmeter)
        end
    else
        # Use all bipartitions
        @inbounds for j in 1:ntrees, i in 1:(j - 1)
            bpindices1 .= bptable.records[:, i]
            bpindices2 .= bptable.records[:, j]
            firstbp1 = findfirst(!iszero, bpindices1)
            firstbp2 = findfirst(!iszero, bpindices2)
            bplist1 = @view bptable.bipartitions[bpindices1[firstbp1:end]]
            bplist2 = @view bptable.bipartitions[bpindices2[firstbp2:end]]
    
            D[i, j] = clustering_information_dist!(M, bplist1, bplist2, ntax)
            next!(progmeter)
        end
    end

    finish!(progmeter)

    return D
end


function clustering_information_dist(tree1::AbstractTree, tree2::AbstractTree; diffonly=true)
    ntax = length(tree1.taxonset)

    bptable = bipartition_table((tree1, tree2))

    if diffonly
        # Select only the bipartitions that are not shared between the two trees
        unique_indices1 = zeros(UInt32, size(bptable.records, 1))
        unique_indices2 = zeros(UInt32, size(bptable.records, 1))

        bpindices1 = @view bptable.records[:, 1]
        bpindices2 = @view bptable.records[:, 2]

        nbp1, nbp2 =
            special_setdiff!(unique_indices1, unique_indices2, bpindices1, bpindices2)

        M = zeros(nbp1, nbp2)          # Cost matrix for bipartition matching
        bplist1 = @view bptable.bipartitions[unique_indices1[1:nbp1]]
        bplist2 = @view bptable.bipartitions[unique_indices2[1:nbp2]]

        d = clustering_information_dist!(M, bplist1, bplist2, ntax)
    else
        # Use all bipartitions
        bpindices1 = @view bptable.records[:, 1]
        bpindices2 = @view bptable.records[:, 2]
        nbp1 = findfirst(!iszero, bpindices1)
        nbp2 = findfirst(!iszero, bpindices2)
        bplist1 = @view bptable.bipartitions[bpindices1[1:nbp1]]
        bplist2 = @view bptable.bipartitions[bpindices2[1:nbp2]]

        d = clustering_information_dist!(M, bplist1, bplist2, ntax)
    end

    return d
end


"""
Return matrix with ``s * n`` value of each bipartition matching.
"""
function mci_cost_matrix!(M, bplist1, bplist2, lrows, lcols, ntax)
    @inbounds for colnum in 1:lcols
        bp1 = bplist2[colnum]
        A1 = mapreduce(count_ones, +, bp1.v.chunks)
        B1 = ntax - A1
        log2A1 = log2(A1)
        log2B1 = log2(B1)
        log2ntax = log2(ntax)
        
        @inbounds for rownum in 1:lrows
            bp2 = bplist1[rownum]
            A2 = mapreduce(count_ones, +, bp2.v.chunks)
            B2 = ntax - A2

            iA1A2 = 0
            for (chunk1, chunk2) in zip(bp1.v.chunks, bp2.v.chunks)
                iA1A2 += count_ones(chunk1 & chunk2)
            end
            iA1B2 = A1 - iA1A2
            iB1A2 = A2 - iA1A2
            iB1B2 = B1 - iB1A2

            log2A2 = log2(A2)
            log2B2 = log2(B2)
    
            s  =
                ifelse(iA1A2 > 0, iA1A2 * (log2ntax + log2(iA1A2) - log2A1 - log2A2), 0.0)
            s +=
                ifelse(iA1B2 > 0, iA1B2 * (log2ntax + log2(iA1B2) - log2A1 - log2B2), 0.0)
            s +=
                ifelse(iB1A2 > 0, iB1A2 * (log2ntax + log2(iB1A2) - log2B1 - log2A2), 0.0)
            s +=
                ifelse(iB1B2 > 0, iB1B2 * (log2ntax + log2(iB1B2) - log2B1 - log2B2), 0.0)
    
            M[rownum, colnum] = s
        end
    end

    return M
end
