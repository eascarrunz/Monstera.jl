"""
    is_logging(io)

Helper function from https://github.com/timholy/ProgressMeter.jl
"""
is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")


"""
    special_countdiff(v1, v2)

Count the number of elements that two sorted vectors do not have in common, fast.

Zero-valued elements are ignored.

## Assumptions:
- `v1` and `v2` are sorted
- `v1` and `v2` have the same length
"""
function special_countdiff(v1, v2)
    n = length(v1)
    i, j = 1, 1
    d = 0
    reached_end_of_one_vector = false

    while i ≤ n && j ≤ n
        # Skip zeros
        while v1[i] == 0
            if i == n
                reached_end_of_one_vector = true
                break
            end
            i += 1
        end

        while @inbounds v2[j] == 0
            if j == n
                reached_end_of_one_vector = true
                break
            end
            j += 1
        end
        reached_end_of_one_vector && break

        # Add to result if the elements are not the same
        if @inbounds v1[i] < v2[j]
            d += 1
            i += 1
        elseif @inbounds v1[i] > v2[j]
            d += 1
            j += 1
            
        # Same element, move on
        else
            i += 1
            j += 1
        end
    end
    
    while i ≤ n
        @inbounds v1[i] == 0 && break
        d += 1
        i += 1
    end
    
    while j ≤ n
        @inbounds v2[j] == 0 && break
        d += 1
        j += 1
    end

    return d
end


"""
    special_setdiff!(result1, result2, v1, v2)

Find the set difference between two sorted vectors `v1` and `v2`, fast.

The elements present in `v1` but not in `v2` are stored in-place in `result1`. Vice-versa
for `v2` and `result2`.

Zero-valued elements are ignored.

## Assumptions:
- `v1` and `v2` are sorted
- `v1` and `v2` have the same length
"""
function special_setdiff!(res1, res2, v1, v2)
    n = length(v1)
    i, j = 1, 1
    r1, r2 = 0, 0
    reached_end_of_one_vector = false

    while i ≤ n && j ≤ n
        # Skip zeros
        while @inbounds v1[i] == 0
            if i == n
                reached_end_of_one_vector = true
                break
            end
            i += 1
        end

        while @inbounds v2[j] == 0
            if j == n
                reached_end_of_one_vector = true
                break
            end
            j += 1
        end

        reached_end_of_one_vector && break

        # Add to result if the elements are not the same
        if @inbounds v1[i] < v2[j]
            r1 += 1
            @inbounds res1[r1] = v1[i]
            i += 1
        elseif @inbounds v1[i] > v2[j]
            r2 += 1
            @inbounds res2[r2] = v2[j]
            j += 1
            
        # Same element, move on
        else
            i += 1
            j += 1
        end
    end
    
    if i ≤ n && @inbounds v1[i] ≠ 0
        while i ≤ n
            r1 += 1
            @inbounds res1[r1] = v1[i]
            i += 1
        end
    end
    
    if j ≤ n && @inbounds v2[j] ≠ 0
        while j ≤ n
            r2 += 1
            @inbounds res2[r2] = v2[j]
            j += 1
        end
    end

    return r1, r2
end
