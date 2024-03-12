#=
MIT License

Copyright (c) 2021 Matthew Griffiths <matthewghgriffiths@gmail.com> and contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

https://github.com/matthewghgriffiths/Assignment.jl
=#

"""
Stores the solution of an assigment problem

Can be directly used to index the appropriate values from the cost matrix 
"""
struct AssignmentSolution{T, S}
	row4col::Vector{Int}
	col4row::Vector{Int}
	cost::Union{T, Nothing}
	v::Vector{S}
	u::Vector{S}
	use_row::Bool

    # AssignmentSolution(row4col, col4row, cost::Nothing, u::Vector{S}, v::Vector{S}, use_row) where S = new{S, S}(row4col, col4row, cost, u, v, use_row)
    AssignmentSolution{T, S}(row4col, col4row, cost::T, u::Vector{S}, v::Vector{S}, use_row) where {T, S} = new{T, S}(row4col, col4row, cost, u, v, use_row)
    AssignmentSolution{T, S}(row4col, col4row, cost::Nothing, u::Vector{S}, v::Vector{S}, use_row) where {T, S} = new{T, S}(row4col, col4row, cost, u, v, use_row)
    AssignmentSolution(row4col, col4row, cost::T, u::Vector{S}, v::Vector{S}, use_row) where {T, S} = new{T, S}(row4col, col4row, cost, u, v, use_row)
end

function AssignmentSolution{T, S}(row4col, col4row, cost::Union{T, Nothing}, u::Vector{S}, v::Vector{S}) where {T, S}
	n = length(row4col)
	m = length(col4row)
	@assert n == length(v) 
	@assert m == length(u)
	AssignmentSolution{T, S}(row4col, col4row, cost, u, v, n <=m)
end

AssignmentSolution(row4col::Vector{Int}, col4row::Vector{Int}, cost::T, u::Vector{S}, v::Vector{S}) where {T, S} = 
    AssignmentSolution{isnothing(cost) ? S : T, S}(row4col, col4row, cost, u, v)

Base.adjoint(sol::AssignmentSolution{T, S}) where {T, S} = AssignmentSolution{T, S}(
	sol.col4row, sol.row4col, sol.cost, sol.u, sol.v, length(sol.col4row) <= length(sol.row4col))

Base.getindex(sol::AssignmentSolution, i::Integer) = sol.use_row ? CartesianIndex(i, sol.row4col[i]) : CartesianIndex(sol.col4row[i], i)

Base.eltype(sol::AssignmentSolution) = CartesianIndex{2}
Base.length(sol::AssignmentSolution) = sol.use_row ? length(sol.row4col) : length(sol.col4row)

Base.iterate(sol::AssignmentSolution) = Base.iterate(sol, 0)

function Base.iterate(sol::AssignmentSolution, i)
	i += 1
	i > length(sol) ? nothing : (sol[i], i)
end

Base.getindex(A::AbstractMatrix, sol::AssignmentSolution) = A[collect(sol)]

Base.show(io::IO, sol::AssignmentSolution) = if sol.use_row
    print(io, "AssignmentSolution(CartesianIndex.(1:$(length(sol.u)), $(sol.row4col)), $(sol.cost))")
else
    print(io, "AssignmentSolution(CartesianIndex.($(sol.col4row), 1:$(length(sol.v))), $(sol.cost))")
end

"""
find_best_assignment(C, maximize=false) = solution

Solve the two-dimensional assignment problem with a rectangular cost matrix C, 
scanning row-wise.

Note that cost returned can overflow if using smaller integer types

# Example 
```julia
julia> M=rand(1:100,3,4)
3×4 Matrix{Int64}:
 77  51  42  67
 72  53  47   4
 24  50  77  96

 julia> sol = find_best_assignment(M)
 AssignmentSolution(CartesianIndex.(1:3, [3, 4, 1]), 70)
 
 julia> sum(M[sol])
 70
 
 julia> max_sol = find_best_assignment(M', true)
 AssignmentSolution(CartesianIndex.([1, 2, 4], 1:3), 226)
```

This code is a port of Matlab code released [1] in the public domain by the
US Naval Research Laboratory. 

This work is not affliated with or endorsed by the US Naval Research Laboratory.

The algorithm is described in detail in [2].

REFERENCES:
[1] D. F. Crouse, "The Tracker Component Library: Free Routines for Rapid 
   Prototyping," IEEE Aerospace and Electronic Systems Magazine, vol. 32, 
   no. 5, pp. 18-27, May. 2017
[2] D. F. Crouse, "Advances in displaying uncertain estimates of multiple
   targets," in Proceedings of SPIE: Signal Processing, Sensor Fusion, and
   Target Recognition XXII, vol. 8745, Baltimore, MD, Apr. 2013
"""
function find_best_assignment(cost_matrix::AbstractMatrix{T}, maximize::Bool=false) where T
	   
    num_row, num_col = size(cost_matrix)
	gain = nothing
    
    num_col > num_row && return find_best_assignment(cost_matrix', maximize)'
	
	# The cost matrix must have all non-negative elements for the assignment
	# algorithm to work. This forces all of the elements to be positive. The
	# delta is added back in when computing the gain in the end.
    if maximize
        ΔC = maximum(cost_matrix)
        C = ΔC .- cost_matrix 
    else
        ΔC = minimum(cost_matrix)
        C = cost_matrix .- ΔC
    end

    # These store the assignment as it is made.
    col4row=zeros(Int, num_row)
    row4col=zeros(Int, num_col)
    S = _dual_eltype(T)
    u=zeros(S, num_col) # The dual variable for the columns
    v=zeros(S, num_row) # The dual variable for the rows.
	
	pred = zeros(Int, num_row)
	scanned_cols = falses(num_col)
	scanned_rows = falses(num_row)
	row2scan = zeros(Int, num_row)
	# Need to be able to fill with infinity
	shortest_path_cost = fill(typemax(float(T)), num_row)
	
	# Initially, none of the columns are assigned.
    for curr_unassigned_col = 1:num_col       
        # This finds the shortest augmenting path starting at k and returns
        # the last node in the path.
        sink, pred, u, v = ShortestPath!(
			curr_unassigned_col,u,v,C,col4row,row4col,
			pred,scanned_cols,scanned_rows,row2scan,shortest_path_cost
		)
        
        # If the problem is infeasible, mark it as such and return.
		sink == 0 && return AssignmentSolution{T, S}(row4col, col4row, nothing, v, u, true)
        
        # We have to remove node k from those that must be assigned.
        j = sink
        @inbounds while true
            i = col4row[j] = pred[j]
            j, row4col[i] = row4col[i], j   			
			
            if i == curr_unassigned_col 
				break
			end
        end
    end
	

	gain = T(0)
	@inbounds for curCol=1:num_col
		gain += cost_matrix[row4col[curCol],curCol]
	end
	return AssignmentSolution{T, S}(col4row, row4col, gain, u, v)
end


_dual_eltype(T) = T <: Integer ? signed(promote_type(T, Int32)) : T


function ShortestPath!(
		curr_unassigned_col::Int,
		u::AbstractVector{S},
		v::AbstractVector{S},
		C::AbstractMatrix{T},
		col4row::AbstractVector{<:Integer},
		row4col::AbstractVector{<:Integer}, 
		pred::AbstractVector{<:Integer} = zeros(Int, size(C, 1)), 
		scanned_cols::AbstractVector{<:Bool} = falses(size(C,2)), 
		scanned_rows::AbstractVector{<:Bool} = falses(size(C,2)), 
		row2scan::AbstractVector{<:Integer} = Vector(1:size(C, 1)),
    	shortest_path_cost::AbstractVector = fill(Inf, size(C, 1))
	) where {T, S}
    # This assumes that unassigned columns go from 1:numUnassigned
    num_row, num_col = size(C)
	
	pred .= 0
    # Initially, none of the rows and columns have been scanned.
    # This will store a 1 in every column that has been scanned.
	scanned_cols .= false 
    # This will store a 1 in every row that has been scanned.
	scanned_rows .= false
	row2scan .= 1:num_row # Columns left to scan.
    num_row2scan = num_row
    
    sink=0
    delta=T(0)
	closest_row_scan = -1
    current_col=curr_unassigned_col
    shortest_path_cost .= Inf
    
    while sink==0        
        # Mark the current row as having been visited.
        scanned_cols[current_col] = 1
        
        # Scan all of the columns that have not already been scanned.
        minval = Inf
        @inbounds for current_row_scan in 1:num_row2scan
            current_row=row2scan[current_row_scan]
            reduced_cost =
                delta + C[current_row,current_col] - u[current_col] - v[current_row];
            if reduced_cost < shortest_path_cost[current_row]
                pred[current_row] = current_col;
                shortest_path_cost[current_row]=reduced_cost;
            end
            
            #Find the minimum unassigned column that was
            #scanned.
            if shortest_path_cost[current_row]<minval 
                minval=shortest_path_cost[current_row];
                closest_row_scan=current_row_scan;
            end
        end
                
	   #If the minimum cost column is not finite, then the problem is
	   #not feasible.
        isfinite(minval) || return 0, pred, u, v
        
        closest_row = row2scan[closest_row_scan]
        #Add the column to the list of scanned columns and delete it from
        #the list of columns to scan.
        scanned_rows[closest_row] = 1;
        num_row2scan -= 1;
		@inbounds for i in closest_row_scan:num_row2scan
			row2scan[i] = row2scan[i + 1]
		end
        
        delta=shortest_path_cost[closest_row];
        #If we have reached an unassigned row.
        if col4row[closest_row]==0
            sink=closest_row;
        else
            current_col=col4row[closest_row];
        end
    end
    
    #Dual Update Step
    
    #Update the first row in the augmenting path.
    u[curr_unassigned_col] += delta
	
    #Update the rest of the rows in the agumenting path.
	@inbounds for (i, col) in enumerate(scanned_cols)
		if col != 0 && i != curr_unassigned_col
			u[i] += delta - convert(T, shortest_path_cost[row4col[i]])
		end
	end
    
    # Update the scanned columns in the augmenting path.
	@inbounds for (i, row) in enumerate(scanned_rows)
		if row != 0
			v[i] -= delta - convert(T, shortest_path_cost[i])
		end
	end
	
	return sink, pred, u, v
end
