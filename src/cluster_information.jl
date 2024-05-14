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

function _mutual_clustering_information!(C::AbstractMatrix{Float64}, bplist1, bplist2, lrows, lcols, ntax)
    mci_cost_matrix!(C, bplist1, bplist2, lrows, lcols, ntax)
    optimal_matching = find_best_assignment(@view(C[1:lrows, 1:lcols]), true)
    
    return optimal_matching.cost / ntax
end


function mutual_clustering_information(bplist1, bplist2, ntax = length(bplist1[1]))
    lrows = length(bplist1)
    lcols = length(bplist2)
    C = zeros(lrows, lcols)
    
    return _mutual_clustering_information!(C, bplist1, bplist2, lrows, lcols, ntax)
end


function clustering_information_dist!(
    C::AbstractMatrix{Float64},
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
            return sum(_clustering_entropy(bp) for bp in bplist2) / ntax
        end
    end

    lcols == 0 &&
        return sum(_clustering_entropy(bp) for bp in bplist1) / ntax

    s = _mutual_clustering_information!(C, bplist1, bplist2, lrows, lcols, ntax)
    H = sum(_clustering_entropy(bp) for bp in bplist1) +
        sum(_clustering_entropy(bp) for bp in bplist2)

    # return 0.5 * H / ntax - s
    return H / ntax - s - s
end


function clustering_information_dist(trees; diffonly=true)
    ntrees = length(trees)
    ntax = length(first(trees).taxonset)

    bprecord = bipartition_record(trees)

    bpindices1 = zeros(UInt32, size(bprecord.occurrences, 1))
    bpindices2 = zeros(UInt32, size(bprecord.occurrences, 1))

    C = zeros(ntax - 3, ntax - 3)          # Cost matrix for bipartition matching
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
        unique_indices1 = zeros(UInt32, size(bprecord.occurrences, 1))
        unique_indices2 = zeros(UInt32, size(bprecord.occurrences, 1))

        @inbounds for j in 2:ntrees, i in 1:(j - 1)
            bpindices1 .= bprecord.occurrences[:, i]
            bpindices2 .= bprecord.occurrences[:, j]
            nbp1, nbp2 =
                special_setdiff!(unique_indices1, unique_indices2, bpindices1, bpindices2)
            bplist1 = bprecord.bipartitions[unique_indices1[1:nbp1]]
            bplist2 = bprecord.bipartitions[unique_indices2[1:nbp2]]
    
            D[i, j] = clustering_information_dist!(C, bplist1, bplist2, ntax)

            unique_indices1 .= zero(UInt32)
            unique_indices2 .= zero(UInt32)
            next!(progmeter)
        end
    else
        # Use all bipartitions
        @inbounds for j in 1:ntrees, i in 1:(j - 1)
            bpindices1 .= bprecord.occurrences[:, i]
            bpindices2 .= bprecord.occurrences[:, j]
            firstbp1 = findfirst(!iszero, bpindices1)
            firstbp2 = findfirst(!iszero, bpindices2)
            bplist1 = isnothing(firstbp1) ? TaxonBipartition[] : @view bprecord.bipartitions[bpindices1[firstbp1:end]]
            bplist2 = isnothing(firstbp2) ? TaxonBipartition[] : @view bprecord.bipartitions[bpindices2[firstbp2:end]]
            # bplist1 = @view bprecord.bipartitions[bpindices1[firstbp1:end]]
            # bplist2 = @view bprecord.bipartitions[bpindices2[firstbp2:end]]
    
            D[i, j] = clustering_information_dist!(C, bplist1, bplist2, ntax)
            next!(progmeter)
        end
    end

    finish!(progmeter)

    return D
end


function clustering_information_dist(tree1::AbstractTree, tree2::AbstractTree; diffonly=true)
    ntax = length(tree1.taxonset)

    bprecord = bipartition_record((tree1, tree2))
    C = zeros(ntax - 3, ntax - 3)          # Cost matrix for bipartition matching

    if diffonly
        # Select only the bipartitions that are not shared between the two trees
        unique_indices1 = zeros(UInt32, size(bprecord.occurrences, 1))
        unique_indices2 = zeros(UInt32, size(bprecord.occurrences, 1))

        bpindices1 = @view bprecord.occurrences[:, 1]
        bpindices2 = @view bprecord.occurrences[:, 2]

        nbp1, nbp2 =
            special_setdiff!(unique_indices1, unique_indices2, bpindices1, bpindices2)

        bplist1 = @view bprecord.bipartitions[unique_indices1[1:nbp1]]
        bplist2 = @view bprecord.bipartitions[unique_indices2[1:nbp2]]

        d = clustering_information_dist!(C, bplist1, bplist2, ntax)
    else
        # Use all bipartitions
        bpindices1 = @view bprecord.occurrences[:, 1]
        bpindices2 = @view bprecord.occurrences[:, 2]
        firstbp1 = findfirst(!iszero, bpindices1)
        firstbp2 = findfirst(!iszero, bpindices2)
        bplist1 = isnothing(firstbp1) ? TaxonBipartition[] : @view bprecord.bipartitions[bpindices1[firstbp1:end]]
        bplist2 = isnothing(firstbp2) ? TaxonBipartition[] : @view bprecord.bipartitions[bpindices2[firstbp2:end]]

        d = clustering_information_dist!(C, bplist1, bplist2, ntax)
    end

    return d
end


"""
Return matrix with ``s * n`` value of each bipartition matching.
"""
function mci_cost_matrix!(C, bplist1, bplist2, lrows, lcols, ntax)
    log2ntax = log2(ntax)

    nA2list = count_ones.(bplist1)
    log2nA2list = log2.(nA2list)
    log2nB2list = @. log2(ntax - nA2list)

    @inbounds for colnum in 1:lcols
        bp1 = bplist2[colnum]
        nA1 = mapreduce(count_ones, +, bp1.v.chunks)
        nB1 = ntax - nA1
        log2nA1 = log2(nA1)
        log2nB1 = log2(nB1)
        
        @inbounds for rownum in 1:lrows
            bp2 = bplist1[rownum]
            nA2 = nA2list[rownum]
            log2nA2 = log2nA2list[rownum]
            log2nB2 = log2nB2list[rownum]

            nUA1A2 = 0
            for (chunk1, chunk2) in zip(bp1.v.chunks, bp2.v.chunks)
                nUA1A2 += count_ones(chunk1 & chunk2)
            end
            nUA1B2 = nA1 - nUA1A2
            nUB1A2 = nA2 - nUA1A2
            nUB1B2 = nB1 - nUB1A2

    
            s  =
                ifelse(nUA1A2 > 0, nUA1A2 * (log2ntax + log2(nUA1A2) - log2nA1 - log2nA2), 0.0)
            s +=
                ifelse(nUA1B2 > 0, nUA1B2 * (log2ntax + log2(nUA1B2) - log2nA1 - log2nB2), 0.0)
            s +=
                ifelse(nUB1A2 > 0, nUB1A2 * (log2ntax + log2(nUB1A2) - log2nB1 - log2nA2), 0.0)
            s +=
                ifelse(nUB1B2 > 0, nUB1B2 * (log2ntax + log2(nUB1B2) - log2nB1 - log2nB2), 0.0)
    
            C[rownum, colnum] = s
        end
    end

    return C
end
