const _PLOT_SYMBOLS = (
    hline = '─',
    vline = '│',
    top_elbow = '┌',
    bottom_elbow = '└',
    threeway = '├'
)

function get_largest_id_and_number_of_outers(branchnode, id, n)
    _, node = branchnode
    n += isouter(node)
    id = max(id, node.id)
    for cbranchnode in children(branchnode)
        _, cnode = cbranchnode
        id = max(id, cnode.id)
        id, n = get_largest_id_and_number_of_outers(cbranchnode, id, n)
    end

    return id, n
end


"""
    CharCanvas

A canvas for drawing text plots upon.

Drawing coordinates are in "X,Y" format with the origin at the top right corner.

Canvases are created filled with `\n`, and the first occurence of `\n` in a line causes 
printing to skip the rest of the line and move on to the next. The `paint_white_background`
function can be used to change '\n' to ' ' in selected areas. This awkward system exists to
prevent printing excessive right-padding whitespaces, which makes tree plots more prone to
deforming when resizing the text area.
"""
struct CharCanvas
    m::Matrix{Char}
    h::Int
    w::Int
    hot::BitMatrix

    CharCanvas(h, w) = new(fill('\n', h, w), h, w, falses(h, w))
end


function Base.print(io::IO, canvas::CharCanvas)
    for (char_row, hot_row) in zip(eachrow(canvas.m), eachrow(canvas.hot))
        for (c, h) in zip(char_row, hot_row)
            c == '\n' && continue
            Base.printstyled(io, c; reverse=h, blink=h)
        end
        Base.print(io, '\n')
    end

    return nothing
end


"""
Draw on a canvas. Points outside the canvas will be ignored.
"""
function Base.setindex!(canvas::CharCanvas, v::Char, x, y)
    for i in y
        i > canvas.h && break
        for j in x
            j > canvas.w && break
            setindex!(canvas.m, v, i, j)
        end
    end

    return nothing
end

function Base.setindex!(canvas::CharCanvas, v::AbstractString, x, y)
    viter = eachindex(v)
    nextiter = iterate(viter)
    isnothing(nextiter) && return nothing
    charindex, state = nextiter

    charindex, state = iterate(viter)
    for i in y
        i > canvas.h && break
        for j in x
            j > canvas.w && break
            v[charindex] == '\n' && continue
            setindex!(canvas.m, v[charindex], i, j)
            nextiter = iterate(viter, state)
            isnothing(nextiter) && break
            charindex, state = nextiter
        end
    end

    return nothing
end


"""
Set coordinates in a canvas to appear highlighted (blinking + reversed).
"""
function highight!(canvas::CharCanvas, x, y)
    for i in y
        i > canvas.h && break
        for j in x
            j > canvas.w && break
            setindex!(canvas.hot, true, i, j)
        end
    end

    return nothing
end


Base.getindex(canvas::CharCanvas, x, y) = getindex(canvas.m, y, x)
Base.size(canvas::CharCanvas) = canvas.w, canvas.h


get_highlight_ids(x) = collect(x)
get_highlight_ids(x::Vector{Int}) = x
get_highlight_ids(x::AbstractNode) = Int[x.id]
get_highlight_ids(x::Vector{<:AbstractNode}) = map(node -> node.id, x)

struct TreeDrawing{T}
    x_coords_nodes::Vector{T}
    y_coords_nodes::Vector{T}
    ids::Bool
    innerlabels::Bool
    outerlabels::Bool
    highlight::Vector{Int}
    
    function TreeDrawing{T}(
        n_nodes;
        ids=true,
        outerlabels=true,
        innerlabels=true,
        highlight=Int[]
        ) where T
        x_coords_nodes = zeros(T, n_nodes)
        y_coords_nodes = zeros(T, n_nodes)
        
        return new(
            x_coords_nodes,
            y_coords_nodes,
            ids, innerlabels,
            outerlabels,
            get_highlight_ids(highlight)
            )
    end
end

Base.size(x::TreeDrawing) = (length(x.x_coords_nodes), length(x.y_coords_nodes))

fork_centre(ytop, ybottom) = (ytop + ybottom) * 0.5
fork_centre(ytop::Integer, ybottom::Integer) = (ytop + ybottom) ÷ 2


"""
Replace non-line break characters with from the left up to the x-th character at height y
"""
function paint_white_background!(canvas, x, y)
    y > canvas.h && return nothing
    x = min(x, canvas.w)
    for i in 1:x
        canvas[i, y] = canvas[i, y] == '\n' ? ' ' : canvas[i, y]
    end

    return nothing
end

function compute_y_coords_nodes!(branchnode, drawing::TreeDrawing{T}, y_next_outer) where T
    branch, node = branchnode
    y_coords_nodes = drawing.y_coords_nodes

    # Use `! isnothing(branch)` in case we have a root of degree 1
    if ! isnothing(branch) && isouter(node)
        y_coords_nodes[node.id] = y_next_outer
        
        return y_next_outer + T(2)
    end

    nchildren = 0
    top_child_id = 0
    bottom_child_id = 0
    for cbranchnode in children(branchnode)
        y_next_outer = compute_y_coords_nodes!(cbranchnode, drawing, y_next_outer)
        nchildren += 1
        _, cnode = cbranchnode
        top_child_id += (nchildren < 2) * cnode.id
        bottom_child_id = cnode.id
    end

    y_coords_nodes[node.id] =
        fork_centre(y_coords_nodes[top_child_id], y_coords_nodes[bottom_child_id])

    return y_next_outer
end

function compute_x_coords_nodes!(branchnode, drawing::TreeDrawing, nan_length=0.0, xparent=0.0)
    branch, node = branchnode
    brlen = 0.0
    if isnothing(branch)
        brlen = 0.0
    else
        brlen = branch.length
        brlen = isnan(brlen) ? nan_length : brlen
    end

    drawing.x_coords_nodes[node.id] = xparent + brlen
    
    for cbranchnode in children(branchnode)
        compute_x_coords_nodes!(cbranchnode, drawing, nan_length, drawing.x_coords_nodes[node.id])
    end

    return nothing
end


function draw_textplot(branchnode, drawing, canvas, x, l)
    _, node = branchnode
    y = drawing.y_coords_nodes[node.id]
    
    node_symbol = node.taxon > 0 ? '●' : '○'
    node_draw_width = 1
    
    id_text = drawing.ids ? string('(', node.id, ')') : ""
    node_draw_width += length(id_text)
    
    this_is_outer = isouter(node)
    draw_label = (this_is_outer & drawing.outerlabels) | (! this_is_outer & drawing.innerlabels)
    label_text = draw_label ? node.label : ""
    node_draw_width += length(label_text)
    
    node_text = string(node_symbol, id_text, label_text)
    x_textend = x + node_draw_width

    xchildren = x_textend + 1 + l
    
    nchildren = 0
    y_first_child = 0
    y_last_child = 0
    y_this_child = 0
    for cbranchnode in children(branchnode)
        _, cnode = cbranchnode
        nchildren += 1
        y_this_child = drawing.y_coords_nodes[cnode.id]
        if nchildren < 2
            y_first_child = y_last_child = y_this_child
        end
        
        if y == y_this_child
            canvas[(x_textend + 1):xchildren, y_this_child] = _PLOT_SYMBOLS.hline
        else
            canvas[(x + 1):xchildren, y_this_child] = _PLOT_SYMBOLS.hline
        end
        canvas[x + 1, y_last_child+1:y_this_child] = _PLOT_SYMBOLS.vline
        canvas[x + 1, y_this_child] = _PLOT_SYMBOLS.threeway
        for y_in_fork in y_last_child:y_this_child
            paint_white_background!(canvas, x, y_in_fork)
        end
        draw_textplot(cbranchnode, drawing, canvas, xchildren, 3)
        y_last_child = y_this_child
    end

    if nchildren > 1
        canvas[x + 1, y_first_child] = _PLOT_SYMBOLS.top_elbow
        canvas[x + 1, y_last_child] = _PLOT_SYMBOLS.bottom_elbow
    end
    canvas[(x + 1):x_textend, y] = node_text
    node.id ∈ drawing.highlight && highight!(canvas, (x+1):x_textend, y)

    return nothing
end


"""
    textplot([io::IO=stdout], tree; <keyword arguments>)
    textplot([io::IO=stdout], node; <keyword arguments>)
    textplot([io::IO=stdout], branch => node; <keyword arguments>)

Use characters to plot a tree or subtree.

Warning: Branch lengths are not respected!

# Arguments
- `ids=true`: Show the ids of the nodes (between parentheses).
- `outerlabels=true`: Show the labels of outer nodes.
- `innerlabels=false`: Show the labels of inner nodes.
- `highlight=Int[]`: Vector of nodes or node IDs that will appear highlighted in the plot.

# Examples
```jldoctest
newick = "((A,B)H,((C,D,E)I,(F)K)J)G;"
tree = Newick.parse(RoundaboutTree, newick, TaxonSet())
textplot(tree)

# output

     ┌───────●(3)A
┌────○(2)
│    └───────●(4)B
│
○(1)         ┌───────●(7)C
│            │
│    ┌───────○(6)────●(8)D
│    │       │
└────○(5)    └───────●(9)E
     │
     └───────○(10)────●(11)F
```
"""
function textplot(io::IO, branchnode; kwargs...)
    _, node = branchnode
    max_id, n_outers = get_largest_id_and_number_of_outers(branchnode, 0, 0)
    drawing = TreeDrawing{Int}(max_id; kwargs...)
    if length(neighbours(node)) == 0
        drawing.y_coords_nodes[1] = 1
    else
        compute_y_coords_nodes!(branchnode, drawing, 1)
    end
    canvas = CharCanvas(2 * n_outers - 1, last(displaysize(io)))
    draw_textplot(branchnode, drawing, canvas, 0, 0)

    print(io, canvas)

    return nothing
end

function textplot(io::IO, node::AbstractNode; kwargs...)
    if length(neighbours(node)) == 0
        textplot(io, nothing => node; kwargs...)

        return nothing
    end

    if hasparent(node)
        branch = first(Monstera.parent(node))
        textplot(io, branch => node; kwargs...)
    else
        textplot(io, nothing => node; kwargs...)
    end
    
    return nothing
end
textplot(io::IO, tree::AbstractTree; kwargs...) = textplot(io, nothing => tree.root; kwargs...)
textplot(x; kwargs...) = textplot(stdout, x; kwargs...)
