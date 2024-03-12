const DEFAULT_BRLENGTH_DIGITS = floor(Int, log10(floatmax()))

"""
    Newick.WriterSettings(args)

Settings for writing Newick strings.

# Arguments
- `brlength::Bool=true`: whether to write branch lengths
- `label_quote::Bool=false`: wheteher to write double quote marks around node labels
- `brlength_digits::Int=$(DEFAULT_BRLENGTH_DIGITS)`: branch lengths will be rounded to this number of digits
- `skip_nanlength::Bool=true`: whether to write the length of a branch when it is `NaN`
"""
mutable struct WriterSettings
    brlength::Bool
    label_quote::Bool
    brlength_digits::UInt16
    skip_nanlength::Bool

    WriterSettings(;
        brlength=true,
        label_quote=false,
        brlength_digits=DEFAULT_BRLENGTH_DIGITS,
        skip_nanlength=true
    ) = new(brlength, label_quote, brlength_digits, skip_nanlength)
end

function _write_newick_brlength(io::IO, branch::AbstractBranch, settings)
    settings.brlength || return nothing
    isnan(branch.length) && settings.skip_nanlength && return nothing
    print(io,':' , string(round(branch.length, digits=settings.brlength_digits)))

    return nothing
end

_write_newick_brlength(::IO, ::Nothing, ::WriterSettings) = nothing


function _write_newick(io::IO, branchnode, is_first_child, settings)
    branch, node = branchnode

    is_first_child || print(io, ",")

    is_first_child = true
    for (cbranch, cnode) in children(branchnode)
        cbranch === branch && continue
         is_first_child && print(io, "(")
        _write_newick(
            io,
            cbranch => cnode,
            is_first_child,
            settings
            )
        is_first_child = false
    end

    is_first_child || print(io, ")")

    print(io, settings.label_quote ? "\"" : "", node.label, settings.label_quote ? "\"" : "")

    _write_newick_brlength(io, branch, settings)

    return nothing
end

_write_newick(io::IO, branchnothing::Pair{<:AbstractBranch,Nothing}, _, settings) =
    @warn "skipped branch connected only on one side"

_write_newick_root(io::IO, node::AbstractNode, settings) = 
    _write_newick_root(io, nothing => node, settings)

_write_newick_root(io::IO, tree::AbstractTree, settings) =
    _write_newick_root(io, tree.root, settings)

function _write_newick_root(
    io::IO,
    x::Union{Pair{<:AbstractBranch, <:AbstractNode},Pair{Nothing, <:AbstractNode}},
    settings
    )
    _write_newick(io, x, true, settings)
    Base.write(io, ';')

    return nothing
end

function _write_newick_root(io::IO, collection, settings)
    for item in collection
        _write_newick_root(io, item, settings)
        Base.write(io, '\n')
    end

    return nothing
end


"""
    Newick.write(<file name or io>, tree, <keyword arguments or writer settings object>)
    Newick.write(<file name or io>, node, <keyword arguments or writer settings object>)
    Newick.write(<file name or io>, branch => node, <keyword arguments or writer settings object>)

Write the Newick representation of a tree to an IOBuffer or a file.

The Newick string format can be optionally customised by passing a `NewickWriterSettings` 
object, or with the keyword arguments below.

See also [`Newick.string`](@ref), [`Newick.WriterSettings`](@ref)

# Arguments
- `brlength::Bool=true`: whether to write branch lengths
- `label_quote::Bool=false`: wheteher to write double quote marks around node labels
- `brlength_digits::Int=$(DEFAULT_BRLENGTH_DIGITS)`: branch lengths will be rounded to this\
 number
- `skip_nanlength::Bool=false`: whether to write the length of a branch when it is `NaN`

## Writing to file
- `append=false`: whether to append to the file
- `lock=true`: lock operations for safe multi-threaded access.
"""
function write(io::IO, x, settings::WriterSettings)
    _write_newick_root(io, x, settings)

    return nothing
end

function write(io, x; kwargs...)
    settings = WriterSettings(; kwargs...)
    write(io, x, settings)

    return nothing
end

function write(file::AbstractString, x; append=false, lock=true, kwargs...)
    open(file, append ? "a" : "w"; lock=lock) do io
        write(io, x; kwargs...)
    end
end


"""
    Newick.string(tree, <keyword arguments or writer settings object>)
    Newick.string(node, <keyword arguments or writer settings object>)
    Newick.string(branch => node, <keyword arguments or writer settings object>)

Return a string with the Newick representation of a tree.

The string format can be optionally customised by passing a `WriterSettings` object, or with
 the keyword arguments below.

See also [`Newick.read`](@ref), [`Newick.WriterSettings`](@ref)

# Arguments
- `brlength::Bool=true`: whether to write branch lengths
- `label_quote::Bool=false`: wheteher to write double quote marks around node labels
- `brlength_digits::Int=$(DEFAULT_BRLENGTH_DIGITS)`: branch lengths will be rounded to this number
- `skip_nanlength::Bool=false`: whether to write the length of a branch when it is `NaN`
"""
function string(x, settings::WriterSettings)
    io = IOBuffer()
    write(io, x, settings)

    return String(take!(io))
end

string(x; args...) = string(x, WriterSettings(; args...))
