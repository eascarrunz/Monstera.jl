using Monstera
using BioSequences
using FASTX
using BioSymbols

records = open("data/pines.fas") do f
    collect(FASTAReader(f))
end

struct FitchModel{S}
    "Maximum number of character states"
    nstates::Int
    "Length of the sequence"
    seglength::Int
    "Number of sequence chunks (64 bits)"
    nchunks::Int
    "Character weights"
    weights::Union{Nothing,Vector{Int}}
    "Down-pass Fitch kernel for binary trees"
    down_pass_kernel::Function
    "Up-pass Fitch kernel for binary trees"
    up_pass_kernel::Function
    "Down-pass state vectors"
    dstates::Vector{S}
    "Up-pass state vectors"
    ustates::Vector{S}
    "Backup state vectors"
    bstates::Vector{S}
    "Total parsimony score"
    treescore::UInt64
    nodescore::Vector{UInt64}
end



lowmask = 0x77777777
highmask =  0x88888888
lowmask =  0x7777
highmask = 0x8888
65535
bitstring(0x8888888888888888)
bitstring(0x8000000000000000)

function create_farris_bitmasks(nstates)
    nstates > 64 && error("number of states cannot exceed 64")
    highmask = UInt64(0x1) << (nstates - 1)
    for _ in 1:(64 ÷ nstates)
        highmask |= highmask << nstates
    end
    lowmask = ~ highmask
    
    return highmask, lowmask
end

bitstring.(create_farris_bitmasks(64))


lstate = 0b1100_0011_0011_1110_1100_0011_0011_0001
rstate = 0b0011_0110_0011_1111_0011_0110_0011_0000

highmask, lowmask = create_farris_bitmasks(4)

# lstate = 0b0001_0001_0001_0001_0001_0001_0001_0001_0001_0001_0001_0001_0001_0001_0001_0001
lstate = 0b1111_1111_1111_1111_1111_1111_1111_1111_0001_0001_0001_0001_0001_0001_0001_0001
rstate = 0b1111_1111_1111_1111_1111_1111_1111_1111_0000_0000_0000_0000_0000_0000_0000_0000
lstate = 0b0001_0001_0001_0001_0001_0001_0001_0000_1111_1111_1111_1111_1111_1111_1111_1111
rstate = 0b0000_0000_0000_0000_0000_0001_0001_0000_1111_1111_1111_1111_1111_1111_1101_1111

function fitchcore_xmp(lstate, rstate)
    NSTATES = 4
    c = lstate & rstate    # c for "common" bits (intersection)
    v = highmask & ~((((c & lowmask) + lowmask) | c)) # >> (NSTATES - 1)

    return v
end

x, y =lstate[1], rstate[1]
c = x & y
v = c & lowmask
println(spaced_bitstring(x));println(spaced_bitstring(y));println(spaced_bitstring(v));
v = v + lowmask
println(spaced_bitstring(x));println(spaced_bitstring(y));println(spaced_bitstring(v));
v = v | c
println(spaced_bitstring(x));println(spaced_bitstring(y));println(spaced_bitstring(v));
v = (~v) & highmask
println(spaced_bitstring(x));println(spaced_bitstring(y));println(spaced_bitstring(v));


function fitchcore(lstate, rstate)
    s = lstate & rstate
    return highmask & (~(s | ((s & lowmask) + lowmask)))
end

lstate = rand(UInt64, 1_000)
rstate = rand(UInt64, 1_000)

y = fitchcore(lstate, rstate)
y = fitchcore.(lstate, rstate)
u = fitchcore_xmp.(lstate, rstate)

println(spaced_bitstring(lstate[2]));println(spaced_bitstring(rstate[2]));println(spaced_bitstring(u[2]));

function spaced_bitstring(x::T, by=4) where T <: Integer
    n = sizeof(T) * 8
    y = bitstring(x)
    join((y[i:i+3] for i in 1:by:(n - by)), '_')
end

function compare_binary(a, b)
    aa = bitstring(a)
    bb = bitstring(b)
    for i in 1:4:60
        print(aa[i:i+3], '_')
    end
    print('\n')
    for i in 1:4:60
        print(bb[i:i+3], '_')
    end
end

compare_binary(lstate[2], rstate[2])


foo = (typemax(UInt32) & (y | y >> 33))
count_ones(foo)
bar = count_ones(65535 & ( y | ( y >> 15 ) )) + count_ones(65535 & ( y >> 32 | ( y >> 47 ) ))
bitstring(65535 & ( y | ( y >> 15 ) ))
bitstring(65535 & ( y >> 32 | ( y >> 47 ) ))

count_ones(foo)

lstate = 0b1100_0011_0011_1110
rstate = 0b0011_0110_0011_1111
u = lstate | rstate
bitstring(u)

#y  = HIGH4 & ~(s | ((s & LOW4 ) + LOW4 ) ) ;

bitstring(s)
bitstring(lowmask)
bitstring((s & lowmask))
bitstring((s & lowmask) + lowmask)
bitstring(s | ((s & lowmask) + lowmask))
bitstring(y)
bitstring(y >> 17)
bitstring(y | (y >> 17))
bar = y | (y >> 17)
bitstring(typemax(UInt16) & bar)

foo = typemax(UInt32) & (y | y >> 33)
foo = typemax(UInt8) & (y | y >> 9)
foo = typemax(UInt16) & (y | y >> 17)
count_ones(foo)
bitstring(foo)



"""
Fitch down-pass kernel. Return score, preliminary states, and backup states
"""
function fitch_pass1_kernel(lstate, rstate, datablock)
    s = lstate & rstate
    u = lstate | rstate
    tmp = highmask & ~(s | ((s & lowmask) + lowmask))
    score = typemax(UInt32) & (tmp | tmp >> 33)
    tmp |= tmp - (tmp >> (datablock.nstates - 1))
    ustate = s | (tmp & u)
    bstate = u | tmp

    return score, ustate, bstate
end


"""
Fitch up-pass kernel. Return the final states of the node.
"""
function fitch_pass2_kernel(astate, pstate, bstate)
    x = astate & ~pstate
    x = (x | ((x & lowmask) + lowmask)) & highmask
    x |= x - (x >> (nstates - 1))

    return (pstate & ~x) | (x & (pstate | astate & bstate))
end


function make_weighted_score_counter(nstates)
    for i in ((64 ÷ nstates) - 1):-1:0
        expr = :((y >> $i * nstates) * weights[$nstates + $i])
    end
end

create_site_masks(nstates) = [UInt64(1) << i for i in (0:nstates:(64 - nstates))]


function create_fitch_kernels(nstates)
    #TODO: Prevent variable boxing with a `let` block
    HIGHMASK, LOWMASK = create_farris_bitmasks(nstates)
    SITEMASKS = create_site_masks(nstates)

    f1 = function up_pass_kernel(c1state, c2state)
        s = c1state & c2state
        u = c1state | c2state
        tmp = HIGHMASK & ~(s | ((s & LOWMASK) + LOWMASK))
        score = count_ones(tmp)
        tmp |= tmp - (tmp >> (nstates - 1))
        upstate = s | (tmp & u)
        backupstate = u | tmp
    
        return score, upstate, backupstate
    end

    f2 = function down_pass_kernel(parent_state, upstate, backupstate)
        x = parent_state & ~upstate
        x = (x | ((x & LOWMASK) + LOWMASK)) & HIGHMASK
        x |= x - (x >> (nstates - 1))
        return (upstate & ~x) | (x & (upstate | parent_state & backupstate))
    end

    return f1, f2
end

function fitch_up_pass(branchnode, m)
    _, node = branchnode
    c1brnode, c2brnode = children(branchnode)
    _, c1node = c1brnode
    _, c2node = c2brnode
    @inbounds for i in 1:m.nchunks
        c1state = m.ustates[c1node.dataid][i]
        c2state = m.ustates[c2node.dataid][i]
        score, ustate, bstate = m.up_pass_kernel(c1state, c2state, weights)
        m.treescore += score
        m.ustate[node.dataid][i] = ustate
        m.bstate[node.dataid][i] = bstate
    end
end




