node_symbol(node::AbstractNode) = node.taxon > 0 ? '●' : '○'
node_symbol(::Nothing) = '×'


function Base.show(io::IO, x::AbstractNode)
    labelled = ! isempty(x.label)
    withtaxon = x.taxon > 0

    print(io, string(typeof(x)) * " #" * string(x.id))
    print(io, 
        (withtaxon ? " (taxon " : " (no taxon") *
        (withtaxon ? string(x.taxon) : "") *
        ")"
    )
    print(io, " - \"" ^ labelled, x.label, "\"" ^ labelled)

    return nothing
end


Base.show(io::IO, x::AbstractBranch) =
    print(io, typeof(x), " #", x.id)


function Base.show(io::IO, ::MIME"text/plain", x::AbstractBranch)
    l, r = nodes_flanking(x)
    lid = isnothing(l) ? "" : string(l.id)
    rid = isnothing(r) ? "" : string(r.id)
    show(io, x)
    print(io, ": ")
    print(io, string(lid), ' ', node_symbol(l))
    print(io, "⎯⎯⎯")
    print(io, node_symbol(r), ' ', string(rid))
    print(io, " (length ", string(x.length), ")")

    return nothing
end


function Base.show(io::IO, tree::AbstractTree)
    b, n = size(tree)
    print(io, typeof(tree), ": ", string(n), " nodes, ", string(b), " branches")
end
