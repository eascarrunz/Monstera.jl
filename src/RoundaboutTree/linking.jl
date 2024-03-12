_setnode!(branch::RoundaboutBranch, node, ::RoundaboutLeft) = setfield!(branch, :_left_node, node)
_setnode!(branch::RoundaboutBranch, node, ::RoundaboutRight) = setfield!(branch, :_right_node, node)

function _set_next_branch!(
    branch::RoundaboutBranch,
    next_branch::RoundaboutBranch,
    ::RoundaboutLeft,
    next_side::RoundaboutSide
    )

    setfield!(branch, :_left_next_branch, next_branch)
    setfield!(branch, :_left_next_side, next_side)
end

function _set_next_branch!(branch::RoundaboutBranch,
    next_branch::RoundaboutBranch,
    ::RoundaboutRight,
    next_side::RoundaboutSide
    )

    setfield!(branch, :_right_next_branch, next_branch)
    setfield!(branch, :_right_next_side, next_side)
end


function _link!(branch::RoundaboutBranch, node::AbstractNode, side::RoundaboutSide)
    last_branch = getfield(node, :_last_branch)

    if isnothing(last_branch)
        _set_next_branch!(branch, branch, side, side)
    else
        last_branch_side = getside(last_branch, node)
        first_branch, first_branch_side = get_link_fields(last_branch, last_branch_side)
        # first_branch = _next_branch(last_branch, last_branch_side)
        # first_branch_side = _next_side(last_branch, last_branch_side)

        _set_next_branch!(branch, first_branch, side, first_branch_side)
        
        _set_next_branch!(last_branch, branch, last_branch_side, side)
    end

    setfield!(node, :_last_branch, branch)
    setfield!(node, :_last_branch_side, side)

    _setnode!(branch, node, side)

    return nothing
end


function link!(node1, branch, node2)
    _link!(branch, node1, RoundaboutLeft())
    _link!(branch, node2, RoundaboutRight())

    return nothing
end


function unlink!(branch::RoundaboutBranch, node::RoundaboutNode)
    last_branch = getfield(node, :_last_branch)

    branch_side = getside(branch, node)
    (prev_branch, prev_branch_side), (next_branch, next_branch_side) =
        find_prev_and_next_branches(branch, branch_side)

    _setnode!(branch, nothing, branch_side)
    
    # The unlinked branch is left pointing to itself
    _set_next_branch!(branch, branch, branch_side, branch_side)

    # Heal the branch circle of the node
    _set_next_branch!(prev_branch, next_branch, prev_branch_side, next_branch_side)

    if last_branch ≡ branch
        if prev_branch ≡ branch
            setfield!(node, :_last_branch, nothing)
            setfield!(node, :_last_branch_side, RoundaboutNoSide())
        else
            setfield!(node, :_last_branch, prev_branch)
            setfield!(node, :_last_branch_side, prev_branch_side)
        end
    end

    return nothing
end
