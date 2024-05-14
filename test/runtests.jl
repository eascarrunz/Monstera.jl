using Monstera
using Test
using Random
using DelimitedFiles
import InteractiveUtils: subtypes

function get_concrete_subtypes!(list, T)
    if isconcretetype(T)
        push!(list, T)
        
        return nothing
    end

    for SubT in subtypes(T)
        get_concrete_subtypes!(list, SubT)
    end

    return nothing
end

CONCRETE_TREE_TYPES = Type[]
get_concrete_subtypes!(CONCRETE_TREE_TYPES, AbstractTree)

@testset "Utility functions" begin
    include("util.jl")
end

@testset "Core tree interface" begin
    include("core_interface.jl")
end

@testset "Newick reading and writing" begin
    include("newick.jl")
end

@testset "Tree distances" begin
    include("distances.jl")
end

