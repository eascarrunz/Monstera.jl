@testset "`fast_symdiff`" begin
    N = 1000
    a = [sort(Random.shuffle(1:(2*N))[1:N], rev=true) for _ in 1:20] 
    b = [sort(Random.shuffle(1:(2*N))[1:N], rev=true) for _ in 1:20]

    for (aa, bb) in zip(a, b)
        aa[rand((N ÷ 2):N):end] .= 0
        bb[rand((N ÷ 2):N):end] .= 0
    end

    for (v1, v2) in zip(a, b)
        @test Monstera.fast_symdiff(v1, v2) == sort(symdiff(v1, v2), rev=true)
    end
end
