"""
    fast_symdiff(v1, v2)

Return the symmetric difference between two vectors of integers.

Zeros are ignored.

This is *not* general-purpose alternative to `symdiff`.

# Assumptions 
- Vectors of equal size
- Vectors ordered from greatest to smallest
- Elements are unique within vectors
"""
function fast_symdiff(v1, v2)
    n = length(v1)
    result = zeros(Int, 2 * n)
    i = j = 1
    k = 0
    while true
        if i > n || @inbounds v1[i] == 0
            while j ≤ n && @inbounds v2[j] ≠ 0
                k += 1
                @inbounds result[k] = v2[j]
                j += 1
            end                
            break
        elseif j > n || @inbounds v2[j] == 0
            while i ≤ n && @inbounds v1[i] ≠ 0
                k += 1
                @inbounds result[k] = v1[i]
                i += 1
            end
            break
        end

        if @inbounds v1[i] < v2[j]
            k += 1
            @inbounds result[k] = v2[j]
            j += 1
        elseif @inbounds v1[i] > v2[j]
            k += 1
            @inbounds result[k] = v1[i]
            i += 1
        else
            i += 1
            j += 1
        end
    end

    resize!(result, k)

    return result
end


"""
    two_way_setdiff(v1, v2, result1, result2)

Compute `setdiff(v1, v2)` and store it on `result1`, and compute `setdiff(v2, v1)` and store
 it in `result2`

Returns the indices of the last element stored in each result vector. Result vectors must be
 of the same length as input vectors. Zeros are ignored.

This is *not* general-purpose alternative to `symdiff`.

# Assumptions 
- Vectors of equal size
- Vectors ordered from greatest to smallest
- Elements are unique within vectors
"""
function two_way_setdiff(v1, v2, result1, result2)
    n = length(v1)
    i = j = 1
    r1 = r2 = 0
    while true
        if i > n || @inbounds v1[i] == 0
            while j ≤ n && @inbounds v2[j] ≠ 0
                r2 += 1
                @inbounds result2[r2] = v2[j]
                j += 1
            end                
            break
        elseif j > n || @inbounds v2[j] == 0
            while i ≤ n && @inbounds v1[i] ≠ 0
                r1 += 1
                @inbounds result1[r1] = v1[i]
                i += 1
            end
            break
        end

        if @inbounds v1[i] < v2[j]
            r2 += 1
            @inbounds result2[r2] = v2[j]
            j += 1
        elseif @inbounds v1[i] > v2[j]
            r1 += 1
            @inbounds result1[r1] = v1[i]
            i += 1
        else
            i += 1
            j += 1
        end
    end

    return r1, r2
end


"""
    count_symdiff(v1, v2)

Count the number of elements that two integer vectors do not have in common.

Zeros are ignored.

# Assumptions 
- Vectors of equal size
- Vectors ordered from greatest to smallest
- Elements are unique within vectors
"""
function count_symdiff(v1, v2)
    n = length(v1)
    i = j = 1
    result = 0
    while true
        if i > n || @inbounds v1[i] == 0
            while j ≤ n && @inbounds v2[j] ≠ 0
                result += 1
                j += 1
            end                
            break
        elseif j > n || @inbounds v2[j] == 0
            while i ≤ n && @inbounds v1[i] ≠ 0
                result += 1
                i += 1
            end
            break
        end

        if @inbounds v1[i] < v2[j]
            result += 1
            j += 1
        elseif @inbounds v1[i] > v2[j]
            result += 1
            i += 1
        else
            i += 1
            j += 1
        end
    end

    return result
end


