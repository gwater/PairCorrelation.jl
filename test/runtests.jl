using Test
using PairCorrelation
using FFTW

import PairCorrelation: rfftkspace2

@testset "k-space" begin
    a = zeros(10, 7)
    @test size(collect(rfftkspace2(size(a)..., 2, 3))) == size(rfft(a))
end

@testset "normalization" begin
    p = PairAutoCorrelation(ones(10, 10), 1, 1)
    @test p(0) == 1
    p = PairAutoCorrelation(ones(10, 10), 4, 5)
    @test p(0) == 1
    p = PairAutoCorrelation(2 * ones(10, 10), 4, 5)
    @test 2^2 == p(0)
end
