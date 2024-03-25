using Monstera

include("patristic.jl")

function nj(T::Type{<:AbstractTree}, D₀)
    ntax = size(D₀, 1)
    nnode = 2ntax - 2
    D = deepcopy(D₀)
    # Dlist = typeof(D)[]
    # Qlist = Matrix{Float64}[]

    tree = T(TaxonSet(ntax), nnode)
    for (i, node) in enumerate(tree.nodes[1:ntax])
        node.taxon = Int32(i)
    end

    candidates = tree.nodes[1:ntax]    # Nodes ready to be added to the tree
    n = ntax                           # Number of nodes ready to be added to the tree

    Q = zeros(size(D)...)
    br = 0             # Next branch to add to the tree

    δnew(a, b) = D[a, b] / 2 + (sum(D[a, 1:n]) - sum(D[1:n, b])) / 2(n - 2)

    for id_u in (ntax + 1):nnode
        for j in 1:n, i in 1:n
            i == j && continue
            @inbounds Q[i, j] = (n - 2) * D[i, j] - sum(D[i, 1:n]) - sum(D[1:n, j])
        end

        # push!(Qlist, deepcopy(Q))
        ind_f, ind_g = findmin(Q[1:n, 1:n]) |> last .|> Tuple |> extrema
        f = candidates[ind_f]
        g = candidates[ind_g]
        splice!(candidates, (ind_f, ind_g))
        u = tree.nodes[id_u]
        push!(candidates, u)

        δfu = δnew(ind_f, ind_g)
        dfg = D[ind_f, ind_g]
        δgu = dfg - δfu

        br += 1
        tree.branches[br].length = δfu
        link!(u, tree.branches[br], f)

        br += 1
        tree.branches[br].length = δgu
        link!(u, tree.branches[br], g)

        n -= 1

        Dnew = zeros(eltype(D), n, n)

        jnew = 1
        for j in axes(D, 2)
            (j == ind_f || j == ind_g) && continue
            Dnew[n, jnew] = Dnew[jnew, n] = (D[ind_f, j] + D[ind_g, j] - dfg) / 2

            inew = 1
            for i in axes(D, 1)
                (i == ind_f || i == ind_g) && continue

                @inbounds Dnew[inew, jnew] = D[i, j]
                inew += 1
            end
            jnew += 1
        end

        D = Dnew
        # push!(Dlist, deepcopy(D))
    end

    f, g = candidates
    dfg = D[1, 2]

    br += 1
    tree.branches[br].length = Float64(dfg)
    link!(f, tree.branches[br], g)

    return tree
end
