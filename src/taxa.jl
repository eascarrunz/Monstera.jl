"""
    TaxonSet([taxon_names])
    TaxonSet(n::Int)

An ordered collection of taxa with unique names.

TaxonSets can be created empty, or from a vector of name strings, or by giving a number `n`,
which creates taxa successively named "1" through "n".
"""
struct TaxonSet
    list::Vector{String}
    dict::Dict{String,Int32}

    function TaxonSet(list::Vector{String}=String[])
        isempty(list) && return new(list, Dict{String,Int32}())

        list = unique(list)
        any(isempty(list)) && error("empty strings cannot be taxon names")
        dict = Dict{String,Int32}(key => id for (id, key) in enumerate(list))

        return new(list, dict)
    end

    TaxonSet(n::Int) = TaxonSet(string.(1:n))
end




"""
    get(taxonset, name)

Get the taxon ID that matches a name in a taxon set, or return 0 if there is no match.
"""
Base.get(ts::TaxonSet, key::String) = get(ts.dict, key, Int32(0))


"""
    key(taxonset, ids...)

Get the names of the taxa with the given `ids` in the taxon set.
"""
key(ts::TaxonSet, inds...) = ts.list[inds...]


"""
    hasname(taxonset, name)

Check whether a taxon set contains a name.
"""
Base.haskey(ts::TaxonSet, key::String) = haskey(ts.dict, key)

Base.length(ts::TaxonSet) = length(ts.list)


"""
    names(taxonset)

Return the list of names of taxa in a taxon set.
"""
Base.keys(ts::TaxonSet) = ts.list



"""
    add!(taxonset, names...)

Add name(s) to a taxon set and return its id.

Names already in the taxon set will be ignored (skipped or return 0).
"""
function add!(ts::TaxonSet, x::String...)
    n = n_old = length(ts)
    for key in x
        isempty(key) && error("a taxon name must be a non-empty `String`")
        haskey(ts, key) && continue
        n += 1
        ts.list[key] = n
    end

    return (n_old+1):n
end

function add!(ts::TaxonSet, key::String)
    isempty(key) && error("a taxon name must be a non-empty `String`")
    haskey(ts, key) && return 0
    id = length(ts) + 1
    push!(ts.list, key)
    ts.dict[key] = id
    
    return id
end

function Base.show(io::IO, x::TaxonSet)
    print(io, "TaxonSet: ", length(x), " taxa")
    show(io, keys(x))
end
