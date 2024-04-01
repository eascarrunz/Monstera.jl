@testset "Custom setdiff functions" begin
    v1 = [1, 2, 3, 4, 5, 6]
    v2 = [1, 2, 3, 7, 8, 9]
    res1 = zeros(Int, length(v1))
    res2 = zeros(Int, length(v2))
    r1, r2 = Monstera.special_setdiff!(res1, res2, v1, v2)
    @test r1 == 3
    @test r2 == 3
    @test res1 == [4, 5, 6, 0, 0, 0]
    @test res2 == [7, 8, 9, 0, 0, 0]
    res1 = zeros(Int, length(v1))
    res2 = zeros(Int, length(v2))
    r1, r2 = Monstera.special_setdiff!(res1, res2, v2, v1)
    @test r1 == 3
    @test r2 == 3
    @test res1 == [7, 8, 9, 0, 0, 0]
    @test res2 == [4, 5, 6, 0, 0, 0]
    @test Monstera.special_countdiff(v1, v2) == 6
    @test Monstera.special_countdiff(v2, v1) == 6


    v1 = [1, 2, 3, 4, 5, 6]
    v2 = [4, 5, 6, 7, 8, 9]
    res1 = zeros(Int, length(v1))
    res2 = zeros(Int, length(v2))
    r1, r2 = Monstera.special_setdiff!(res1, res2, v1, v2)
    @test r1 == 3
    @test r2 == 3
    @test res1 == [1, 2, 3, 0, 0, 0]
    @test res2 == [7, 8, 9, 0, 0, 0]
    res1 = zeros(Int, length(v1))
    res2 = zeros(Int, length(v2))
    r1, r2 = Monstera.special_setdiff!(res1, res2, v2, v1)
    @test r1 == 3
    @test r2 == 3
    @test res1 == [7, 8, 9, 0, 0, 0]
    @test res2 == [1, 2, 3, 0, 0, 0]
    @test Monstera.special_countdiff(v1, v2) == 6
    @test Monstera.special_countdiff(v2, v1) == 6


    v1 = [0, 0, 3, 4, 5, 6]
    v2 = [0, 0, 0, 7, 8, 9]
    res1 = zeros(Int, length(v1))
    res2 = zeros(Int, length(v2))
    r1, r2 = Monstera.special_setdiff!(res1, res2, v1, v2)
    @test r1 == 4
    @test r2 == 3
    @test res1 == [3, 4, 5, 6, 0, 0]
    @test res2 == [7, 8, 9, 0, 0, 0]
    res1 = zeros(Int, length(v1))
    res2 = zeros(Int, length(v2))
    r1, r2 = Monstera.special_setdiff!(res1, res2, v2, v1)
    @test r1 == 3
    @test r2 == 4
    @test res1 == [7, 8, 9, 0, 0, 0]
    @test res2 == [3, 4, 5, 6, 0, 0]
    @test Monstera.special_countdiff(v1, v2) == 7
    @test Monstera.special_countdiff(v2, v1) == 7


    v1 = [0, 0, 3, 4, 5, 6]
    v2 = [0, 0, 0, 4, 8, 9]
    res1 = zeros(Int, length(v1))
    res2 = zeros(Int, length(v2))
    r1, r2 = Monstera.special_setdiff!(res1, res2, v1, v2)
    @test r1 == 3
    @test r2 == 2
    @test res1 == [3, 5, 6, 0, 0, 0]
    @test res2 == [8, 9, 0, 0, 0, 0]
    res1 = zeros(Int, length(v1))
    res2 = zeros(Int, length(v2))
    r1, r2 = Monstera.special_setdiff!(res1, res2, v2, v1)
    @test r1 == 2
    @test r2 == 3
    @test res1 == [8, 9, 0, 0, 0, 0]
    @test res2 == [3, 5, 6, 0, 0, 0]
    @test Monstera.special_countdiff(v1, v2) == 5
    @test Monstera.special_countdiff(v2, v1) == 5


    v1 = [1, 2, 3, 4, 5, 6]
    v2 = [0, 0, 0, 0, 0, 0]
    res1 = zeros(Int, length(v1))
    res2 = zeros(Int, length(v2))
    r1, r2 = Monstera.special_setdiff!(res1, res2, v1, v2)
    @test r1 == 6
    @test r2 == 0
    @test res1 == [1, 2, 3, 4, 5, 6]
    @test res2 == [0, 0, 0, 0, 0, 0]
    res1 = zeros(Int, length(v1))
    res2 = zeros(Int, length(v2))
    r1, r2 = Monstera.special_setdiff!(res1, res2, v2, v1)
    @test r1 == 0
    @test r2 == 6
    @test res1 == [0, 0, 0, 0, 0, 0]
    @test res2 == [1, 2, 3, 4, 5, 6]
    @test Monstera.special_countdiff(v1, v2) == 6
    @test Monstera.special_countdiff(v2, v1) == 6


    v1 = [0, 0, 0, 0, 0, 0]
    v2 = [0, 0, 0, 0, 0, 0]
    res1 = zeros(Int, length(v1))
    res2 = zeros(Int, length(v2))
    r1, r2 = Monstera.special_setdiff!(res1, res2, v1, v2)
    @test r1 == 0
    @test r2 == 0
    @test res1 == [0, 0, 0, 0, 0, 0]
    @test res2 == [0, 0, 0, 0, 0, 0]
    @test Monstera.special_countdiff(v1, v2) == 0
end
