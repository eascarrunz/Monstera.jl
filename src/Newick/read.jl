const SIZE_NODE_BUFFER = Int32(250)    # Greater than zero


"""
    Newick.ReaderSettings(taxonset = TaxonSet(); taxafrom=:outer)

Settings for reading Newick strings.

Node labels can be interpreted as taxon names. Taxon names will be mapped to the 
corresponding taxon in the taxon set, or added to the taxon set if they were not already 
present.

The `taxafrom` keyword argument is used to customise which node labels will be interpeted as
taxon names, with any of the following values:

- `:outer` (default) only outer nodes
- `:inner` only inner nodes
- `:all` all nodes
- `:none` no nodes
"""
mutable struct ReaderSettings
    taxonset::Union{TaxonSet,Nothing}
    taxafrom::Symbol

    ReaderSettings(taxonset=nothing; taxafrom=:outer) = new(taxonset, taxafrom)
end


struct TreeParseState{T<:AbstractTree}
    tree::T
    last_node_ind::Int32
    last_branch_ind::Int32
    nnode::Int32
    nbranch::Int32
    settings::ReaderSettings

    TreeParseState(tree::T, i, j, nnode, nbranch, settings) where T =
        new{T}(tree, i, j, nnode, nbranch, settings)
end


function assign_taxon!(node, taxonset::TaxonSet)
    taxon_id = get(taxonset, node.label)
    taxon_id = taxon_id == 0 ? add!(taxonset, node.label) : taxon_id
    node.taxon = Int32(taxon_id)

    return nothing
end


assign_taxon!(node, settings::ReaderSettings) =
    assign_taxon!(node, settings.taxonset, Val(settings.taxafrom))


function assign_taxon!(node, taxonset::TaxonSet, ::Val{:outer})
    isouter(node) || return nothing
    assign_taxon!(node, taxonset)

    return nothing
end

function assign_taxon!(node, taxonset::TaxonSet, ::Val{:inner})
    isinner(node) || return nothing
    assign_taxon!(node, taxonset)

    return nothing
end

assign_taxon!(node, taxonset::TaxonSet, ::Val{:all}) = assign_taxon!(node, taxonset)

assign_taxon!(node, ::TaxonSet, ::Val{:none}) = nothing



function _get_free_node(treestate::TreeParseState)
    NType = nodetype(typeof(treestate.tree))
    curr_node_ind = treestate.last_node_ind + 1
    new_nnode = treestate.nnode
    if curr_node_ind > treestate.nnode
        new_nnode += SIZE_NODE_BUFFER
        resize!(treestate.tree.nodes, new_nnode)
        for i in (treestate.nnode + 1):new_nnode
            treestate.tree.nodes[i] = NType(i)
        end
    end

    new_state = TreeParseState(
        treestate.tree,
        curr_node_ind,
        treestate.last_branch_ind,
        new_nnode,
        treestate.nbranch,
        treestate.settings
    )

    return (new_state.tree.nodes[new_state.last_node_ind], new_state)
end

function _get_free_branch(treestate::TreeParseState)
    curr_branch_ind = treestate.last_branch_ind + 1
    new_nbranch = treestate.nbranch
    if curr_branch_ind > treestate.nbranch
        new_nbranch += SIZE_NODE_BUFFER
        resize!(treestate.tree.branches, new_nbranch)
        for i in (treestate.nbranch + 1):new_nbranch
            treestate.tree.branches[i] = branchtype(typeof(treestate.tree))(i)
        end
    end

    new_state = TreeParseState(
        treestate.tree,
        treestate.last_node_ind,
        curr_branch_ind,
        treestate.nnode,
        new_nbranch,
        treestate.settings
    )

    return (new_state.tree.branches[new_state.last_branch_ind], new_state)
end

_flush!(treestate::TreeParseState) =
    resize!(treestate.tree, treestate.last_node_ind, treestate.last_branch_ind)


"""
Parse a quoted string after Base.reading the first quote (`qchar` to specify '\'' vs '\"'), with escape character `escchar`
"""
function parse_quoted(io::IO, qchar::Char, escchar::Char)
    value_io = IOBuffer()
    c = Base.read(io, Char)

    while c ≠ qchar
        eof(io) && error("reached end-of-stream in unterminated quoted string")
        if c == escchar
            Base.write(value_io, c)
            c = Base.read(io, Char)
            Base.write(value_io, c)
            c = Base.read(io, Char)
        end
        Base.write(value_io, c)
        c = Base.read(io, Char)
    end

    return String(take!(value_io))
end

"""
After reading ':'
"""
function parse_branch!(io::IO, branch::AbstractBranch)
    value_io = IOBuffer()
    c = peek(io, Char)
    while c ∉ (',', ')', ';')
        if c == '['
            branch.data = parse_comment(io, branch.data)
        end
        Base.read(io, Char)
        Base.write(value_io, c)
        c = peek(io, Char)
    end
    branch.length = Base.parse(Float64, String(take!(value_io)))

    return nothing
end

function parse_numeric(str)
    value = try
        Base.parse(Float64, str)
    catch
        tryparse(Int, str)
    end

    return value
end


"""
After reading '['. Doesn't work with arrays yet!
"""
function parse_comment(io, target)
    data = target.data
    @assert Base.read(io, Char) == '&'    # NEXUS quirk
    key_io = IOBuffer()
    value_io = IOBuffer()
    key_str = ""
    value_str = ""
    state = :key    # or :value

    c = Base.read(io, Char)
    while ! eof(io)
        if c == ']'
            break
        elseif isspace(c)
            continue
        elseif c == '\'' || c == '\"'
            if state == :key
                key_str = parse_quoted(io, c, '\\')
            elseif state == :value
                value_str = parse_quote(io, c, '\\')
            end
        elseif c == '='
            key_str = String(take!(key_io))
            state = :value
        elseif c == ','
            value_str = String(take!(value_io))
            value_as_numeric = parse_numeric(value_str)
            value = isnothing(value_as_numeric) ?
                value_str : value_as_numeric
            data[key] = value

            key_str = value_str = ""
            state = :key
        else
            if state == :key
                Base.write(key_io, c)
            elseif state == :value
                Base.write(value_io, c)
            end
        end
    end

    return data
end

function _label_and_link!(pnode, branch, cnode, node_label, node_label_io)
    node_label = isempty(node_label) ? String(take!(node_label_io)) : node_label
    cnode.label = node_label

    link!(pnode, branch, cnode)
end

function parse_node(io, pnode, treestate)
    node, treestate = _get_free_node(treestate)

    # Fix this: Not type stable!
    if ! isnothing(pnode)
        branch, treestate = _get_free_branch(treestate)
    else
        branch = nothing
    end

    node_label_io = IOBuffer()
    node_label = ""

    while ! eof(io)
        c = peek(io, Char)

        if isspace(c)
            c = Base.read(io, Char)
            continue
        elseif c == '['
            parse_comment(io, node)
        elseif c == ')'
            Base.read(io, Char)
            break
        elseif c == '('
            Base.read(io, Char)
            treestate = parse_node(io, node, treestate)
            c = peek(io, Char)
        elseif c == ','
            _label_and_link!(pnode, branch, node, node_label, node_label_io)
            Base.read(io, Char)
            isnothing(treestate.settings.taxonset) || assign_taxon!(node, treestate.settings)
            treestate = parse_node(io, pnode, treestate)
            return treestate
        elseif c == ':'
            if isnothing(branch)
                branch, treestate = _get_free_branch(treestate)
            end
            Base.read(io, Char)
            parse_branch!(io, branch)
        elseif c == ';'
            node_label = isempty(node_label) ? String(take!(node_label_io)) : node_label
            node.label = node_label
            
            if isnothing(pnode) && ! isnothing(branch)
                error("the Newick string contains a branch attached at the root. This is not supported.")
            end
            
            isnothing(treestate.settings.taxonset) || assign_taxon!(node, treestate.settings)

            return treestate
        elseif c == '\'' || c == '\"'
            if isempty(node_label)
                node_label = parse_quoted(io, c, '\\')
            else
                error("quote character `$(c)` found in the middle of unquoted string at position $(io.ptr)")
            end
        else
            Base.write(node_label_io, c)
            Base.read(io, Char)
        end
    end

    node_label = isempty(node_label) ? String(take!(node_label_io)) : node_label
    node.label = node_label

    link!(pnode, branch, node)

    isnothing(treestate.settings.taxonset) || assign_taxon!(node, treestate.settings)

    return treestate
end


function read(io::IO, Ttree::Type, settings::ReaderSettings)
    while ! eof(io) && isspace(peek(io, Char))
        Base.read(io, Char)
    end

    eof(io) && error("there is no Newick expression in the string")

    tree = Ttree(SIZE_NODE_BUFFER)
    treestate = TreeParseState(tree, 0, 0, SIZE_NODE_BUFFER, SIZE_NODE_BUFFER - 1, settings)

    treestate = parse_node(io, nothing, treestate)
    Base.read(io, Char) == ';' || error("Newick string not terminated in ';' (position $(io.ptr))")

    treestate.tree.root = treestate.tree.nodes[1]
    _flush!(treestate)

    treestate.tree.taxonset = settings.taxonset

    return treestate.tree
end

function read(io::IO, ::Type{Vector{T}}, settings::ReaderSettings) where T <: AbstractTree
    trees = T[]

    while ! eof(io)
        tree = read(io, T, settings)
        push!(trees, tree)
        while ! eof(io) && isspace(peek(io, Char))
            Base.read(io, Char)
        end
    end

    return trees
end


"""
    Newick.parse(T, string, settings)

Parse a string in Newick format as a tree of type `T`.
"""
function parse(Ttree::Type, str::String, settings::ReaderSettings)
    io = IOBuffer(str)

    return read(io, Ttree, settings)
end

parse(Ttree::Type, str::String, taxonset=nothing; kwargs...) =
    parse(Ttree, str, ReaderSettings(taxonset; kwargs...))


"""
    Newick.read(file, sink)
    Newick.read(io, sink)

Read trees from a file in Newick format.

The `sink` argument specifies the concrete type that will be used to represent the trees. It
can correspond to an `AbstractTree` subtype or a vector thereof. If a tree type is given,
only the first Newick string in the file will be read. If a vector of a tree type is given, 
all the Newick strings in the file will be read.
"""
function read(file::AbstractString, sink::Type, settings::ReaderSettings)
    result = open(file, "r") do f
        read(f, sink, settings)
    end

    return result
end

read(file::AbstractString, sink::Type, taxonset=nothing; args...) =
    read(file, sink, ReaderSettings(taxonset; args...))
