abstract type RoundaboutSide end
struct RoundaboutLeft <: RoundaboutSide end
struct RoundaboutRight <: RoundaboutSide end
struct RoundaboutNoSide <: RoundaboutSide end


mutable struct RoundaboutBranch{T} <: AbstractBranch
    id::Int32
    length::Float64
    _left_node::Union{Nothing,T}
    _left_next_branch::Union{Nothing,RoundaboutBranch{T}}
    _left_next_side::RoundaboutSide
    _right_node::Union{Nothing,T}
    _right_next_branch::Union{Nothing,RoundaboutBranch{T}}
    _right_next_side::RoundaboutSide
end


Base.show(io::IO, x::RoundaboutBranch) =
    print(io, "RoundaboutBranch #", x.id)


const ROUNDABOUTBRANCH_PRIVATE_PROPERTIES = (
    :_left_node,
    :_left_next_branch,
    :_left_next_side,
    :_right_node,
    :_right_next_branch,
    :_right_next_side,
)


function Base.propertynames(::RoundaboutBranch, private::Bool=false)
    namelist = fieldnames(RoundaboutBranch)
    namelist = private ? intersect(namelist, ROUNDABOUTBRANCH_PRIVATE_PROPERTIES) : namelist

    return namelist
end


function Base.getproperty(branch::RoundaboutBranch, key::Symbol)
    key == :id && return getfield(branch, :id)
    key == :length && return getfield(branch, :length)

    pubpropertyerror(Type{RoundaboutBranch}, key)
end


function Base.setproperty!(branch::RoundaboutBranch, key::Symbol, val)
    key == :id && return setfield!(branch, :id, val)
    key == :length && return setfield!(branch, :length, val) 

    pubpropertyerror(Type{RoundaboutBranch}, key)
end


nodes_flanking(branch::RoundaboutBranch) = (
    getfield(branch, :_left_node),
    getfield(branch, :_right_node)
)


mutable struct RoundaboutNode <: AbstractNode
    id::Int32
    taxon::Int32
    label::String

    # Access points to the *last* branch in the circular list of neighbours
    # The last branch is connect to the first branch. This allows quick head and tail
    # insertions.
    _last_branch::Union{Nothing, RoundaboutBranch}
    _last_branch_side::RoundaboutSide

    RoundaboutNode(i) = new(i, 0, "", nothing, RoundaboutNoSide())
end


const RoundaboutBranchNodePair = Pair{RoundaboutBranch{RoundaboutNode},RoundaboutNode}


RoundaboutBranch(i) =
    RoundaboutBranch{RoundaboutNode}(
        i, NaN, nothing, nothing, RoundaboutNoSide(), nothing, nothing, RoundaboutNoSide()
        )


const ROUNDABOUTNODE_PRIVATE_PROPERTIES = (
    :_last_branch,
    :_last_branch_side,
)


function Base.propertynames(::RoundaboutNode, private::Bool=false)
    namelist = fieldnames(RoundaboutNode)
    private || intersect(namelist, ROUNDABOUTNODE_PRIVATE_PROPERTIES)

    return namelist
end


function Base.getproperty(node::RoundaboutNode, key::Symbol)
    key == :id && return getfield(node, :id)
    key == :taxon && return getfield(node, :taxon) 
    key == :label && return getfield(node, :label) 

    pubpropertyerror(Type{RoundaboutNode}, key)
end


function Base.setproperty!(node::RoundaboutNode, key::Symbol, val)
    key == :id && return setfield!(node, :id, val)
    key == :taxon && return setfield!(node, :taxon, val) 
    key == :label && return setfield!(node, :label, val) 

    pubpropertyerror(Type{RoundaboutNode}, key)
end


mutable struct RoundaboutTree <: AbstractTree
    taxonset::Union{Nothing,TaxonSet}
    root::Union{Nothing,RoundaboutNode}
    nodes::Vector{RoundaboutNode}
    branches::Vector{RoundaboutBranch}
    _ind_last_node::Int
    _ind_last_branch::Int

    function RoundaboutTree(taxonset, n)
        n > 0 || error("trees must be created with at least one node")
        b = n - 1
        nodes = [RoundaboutNode(i) for i in 1:n]
        branches = [RoundaboutBranch(i) for i in 1:b]

        new(taxonset, first(nodes), nodes, branches, n, b)
    end

    RoundaboutTree(n) = RoundaboutTree(nothing, n)
end


treetype(::Type{RoundaboutNode}) = RoundaboutTree
treetype(::Type{RoundaboutBranch}) = RoundaboutTree
nodetype(::Type{RoundaboutTree}) = RoundaboutNode
nodetype(::Type{RoundaboutBranch}) = RoundaboutNode
branchtype(::Type{RoundaboutTree}) = RoundaboutBranch
branchtype(::Type{RoundaboutNode}) = RoundaboutBranch


Base.size(tree::RoundaboutTree) = (
    branches = length(getfield(tree, :branches)),
    nodes = length(getfield(tree, :nodes)),
)


function newelem!(
    T::Type,
    tree::RoundaboutTree,
    n,
    elemvec_field_name::Symbol,
    ind_last_elem_field_name::Symbol
    )

    ind_last_elem_current = getfield(tree, ind_last_elem_field_name)
    ind_last_elem_new = ind_last_elem_current + n
    elem_capacity_current = length(getfield(tree, elemvec_field_name))
    if ind_last_elem_new > elem_capacity_current
        ind_elems_to_create = (elem_capacity_current + 1):ind_last_elem_new
        push!(
            getfield(tree, elemvec_field_name),
            (T(i) for i in ind_elems_to_create)...
        )
    end
    ind_new_elems = (ind_last_elem_current + 1):ind_last_elem_new
    setfield!(tree, ind_last_elem_field_name, ind_last_elem_new)

    return (getfield(tree, elemvec_field_name)[i] for i in ind_new_elems)
end


newnode!(tree::RoundaboutTree, n) =
    newelem!(RoundaboutNode, tree, n, :nodes, :_ind_last_node)
newbranch!(tree::RoundaboutTree, n) =
    newelem!(RoundaboutBranch, tree, n, :branches, :_ind_last_branch)


_opposite_side(::RoundaboutLeft) = RoundaboutRight()
_opposite_side(::RoundaboutRight) = RoundaboutLeft()


"""
Return the node that is connected to a given side of the branch
"""
getnode(branch::RoundaboutBranch, ::RoundaboutLeft) = getfield(branch, :_left_node)
getnode(branch::RoundaboutBranch, ::RoundaboutRight) = getfield(branch, :_right_node)
getnode(::RoundaboutBranch, ::RoundaboutNoSide) = error("\"no side\" given")


"""
Return the side of the branch (left, right or no side) that is connected to the node.
"""
getside(branch::RoundaboutBranch, node::RoundaboutNode) =
    getnode(branch, RoundaboutLeft()) ≡ node ? RoundaboutLeft() : RoundaboutRight()


"""
Return the next branch to move on to on the left or right side.
"""
@inline get_next_branch(branch::RoundaboutBranch, ::RoundaboutLeft) =
    getfield(branch, :_left_next_branch)

@inline get_next_branch(branch::RoundaboutBranch, ::RoundaboutRight) =
    getfield(branch, :_right_next_branch)

@inline get_next_branch(::Nothing, ::RoundaboutSide) = nothing


"""
Return the next side to move on to on the left or right side.
"""
@inline get_next_side(branch::RoundaboutBranch, ::RoundaboutLeft) =
    getfield(branch, :_left_next_side)

@inline get_next_side(branch::RoundaboutBranch, ::RoundaboutRight) =
    getfield(branch, :_right_next_side)

@inline get_next_side(::Nothing, ::RoundaboutNoSide) = :_no_side

struct RoundaboutBranchIterator{T<:Union{Nothing,RoundaboutBranch}}
    init_branch::T
    init_side::RoundaboutSide 
    stopbefore::RoundaboutBranch

    RoundaboutBranchIterator(init_branch::RoundaboutBranch, init_side, stopbefore) =
        new{RoundaboutBranch}(init_branch, init_side, stopbefore)

    RoundaboutBranchIterator(::Nothing, _, __) = new{Nothing}()
end

Base.iterate(::RoundaboutBranchIterator{Nothing}) = nothing

function Base.iterate(iter::RoundaboutBranchIterator{RoundaboutBranch})
    next_branch = get_next_branch(iter.init_branch, iter.init_side)
    next_side = get_next_side(iter.init_branch, iter.init_side)
    
    node = getnode(iter.init_branch, _opposite_side(iter.init_side))

    return iter.init_branch => node, (next_branch, next_side)
end

function Base.iterate(
    iter::RoundaboutBranchIterator,
    state::Tuple{RoundaboutBranch,RoundaboutSide}
    )
    branch, side = state
    branch ≡ iter.stopbefore && return nothing
    next_branch = get_next_branch(branch, side)
    next_side = get_next_side(branch, side)
    node = getnode(branch, _opposite_side(side))

    return branch => node, (next_branch, next_side)
end


Iterators.eltype(::RoundaboutBranchIterator) = RoundaboutBranchNodePair

Base.length(::RoundaboutBranchIterator{Nothing}) = 0

function Base.length(iter::RoundaboutBranchIterator{RoundaboutBranch})
    l = 0
    branch = iter.init_branch
    side = iter.init_side
    while true
        prev_side = side
        side = get_next_side(branch, side)
        branch = get_next_branch(branch, prev_side)
        l += 1
        branch ≡ iter.stopbefore && break
    end

    return l
end

function Base.last(iter::RoundaboutBranchIterator{RoundaboutBranch})
    branch = iter.init_branch
    side = iter.init_side

    while true
        next_branch = get_next_branch(branch, side)
        next_side = get_next_side(branch, side)
        next_branch ≡ iter.stopbefore && break
        branch = next_branch
        side = next_side
    end

    node = getnode(branch, _opposite_side(side))

    return branch => node
end


function neighbours(node::RoundaboutNode)
    last_branch = getfield(node, :_last_branch)
    last_side = getfield(node, :_last_branch_side)
    first_branch = get_next_branch(last_branch, last_side)
    first_side = get_next_side(last_branch, last_side)

    return RoundaboutBranchIterator(first_branch, first_side, first_branch)
end


function _children(branch::RoundaboutBranch, side::RoundaboutSide)
    next_branch = get_next_branch(branch, side)
    next_side = get_next_side(branch, side)
    
    iter = branch ≡ next_branch ? 
        RoundaboutBranchIterator(nothing, RoundaboutNoSide(), nothing) :
            RoundaboutBranchIterator(next_branch, next_side, branch)

    return iter
end


function children(branchnode::RoundaboutBranchNodePair)
    branch, node = branchnode
    branch_side = getside(branch, node)

    return _children(branch, branch_side)
end

function Base.resize!(tree::RoundaboutTree, n, b)
    n_old, b_old = length(tree.nodes), length(tree.branches)
    resize!(tree.nodes, n)
    resize!(tree.branches, b)

    @simd for i in (n_old + 1):n
        tree.nodes[i] = RoundaboutNode[i]
    end

    @simd for i in (b_old + 1):b
        tree.branches[i] = RoundaboutBranch[i]
    end

    return nothing
end



function _reindex!(branch_list, node_list, n, b, branchnode)
    branch, node = branchnode
    n += 1
    b += 1

    if n > length(node_list)
        push!(node_list, node)
    else
        node_list[n] = node
    end
    node.id = n

    if b > length(branch_list)
        push!(branch_list, branch)
        branch_list[b] = branch
    end
    branch.id = b

    for cbranchnode in children(branchnode)
        n, b = _reindex!(branch_list, node_list, n, b, cbranchnode)
    end

    return n, b
end


function reindex!(tree::RoundaboutTree)
    new_node_list = similar(tree.nodes)
    new_branch_list = similar(tree.branches)

    _reindex!(new_branch_list, new_node_list, 0, 0, nothing => tree.root)

    tree.nodes = new_node_list
    tree.branches = new_branch_list

    return nothing
end


"""
Return the next branch and side, and the previous branch and side in a branch cycle
"""
function find_prev_and_next_branches(branch::RoundaboutBranch, side::RoundaboutSide)
    # current_branch = branch
    current_branch = get_next_branch(branch, side)
    side = get_next_side(branch, side)
    prev_branch = current_branch
    prev_side = side
    next_side = side
    next_branch = current_branch
    while current_branch ≢ branch
        prev_branch = current_branch
        prev_side = side
        current_branch = get_next_branch(prev_branch, prev_side)
        side = get_next_side(prev_branch, prev_side)
    end

    return (prev_branch, prev_side), (next_branch, next_side)
end


sidestring(::RoundaboutLeft) = "left"
sidestring(::RoundaboutRight) = "right"
sidestring(::RoundaboutNoSide) = "noside"


function isouter(node::RoundaboutNode)
    last_branch = getfield(node, :_last_branch)
    isnothing(last_branch) && return true
    last_side = getfield(node, :_last_branch_side)
    first_branch = get_next_branch(last_branch, last_side)

    return last_branch ≡ first_branch
end


function details(branch::RoundaboutBranch)
    rnode = getfield(branch, :_right_node)
    lnode = getfield(branch, :_left_node)
    rnode_id = isnothing(rnode) ? 0 : rnode.id
    lnode_id = isnothing(lnode) ? 0 : lnode.id
    rside = sidestring(getfield(branch, :_right_next_side))
    lside = sidestring(getfield(branch, :_left_next_side))
    print("Branch #$(getfield(branch, :id)): ")
    print("N$(lnode_id) NB$(getfield(branch, :_left_next_branch).id), S[$(lside)]")
    print("⎯⎯⎯")
    print("N$(rnode_id) NB$(getfield(branch, :_right_next_branch).id), S[$(rside)]")

    return nothing
end

function cycle_details(branch::RoundaboutBranch, side)
    visited_branches = RoundaboutBranch[]
    current_branch = branch
    current_side = side
    while current_branch ∉ visited_branches
        details(current_branch)
        print('\n')
        next_branch = get_next_branch(current_branch, current_side)
        next_side = get_next_side(current_branch, current_side)
        next_branch ≡ branch && break
        push!(visited_branches, current_branch)
        current_branch, current_side = next_branch, next_side
    end

    return nothing
end
