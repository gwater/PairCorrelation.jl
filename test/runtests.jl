using Test
using PairCorrelation
using FFTW

import PairCorrelation: rfftkspace2

@testset "k-space" begin
    a = zeros(10, 7)
    @test size(collect(rfftkspace2(size(a)..., 2, 3))) == size(rfft(a))
end

@testset "normalization" begin
    a = ones(10, 10)
    p = PairAutoCorrelation(a, 1, 1)
    @test p(0) == 1
    p = PairAutoCorrelation(a, 4, 5)
    @test p(0) == 1
    p = PairAutoCorrelation(2a, 4, 5)
    @test 2^2 == p(0)
    p = PairRadialCorrelation(a, a, 1, 1)
    @test p(0) == 1
    p = PairRadialCorrelation(a, a, 4, 5)
    @test p(0) == 1
    p = PairRadialCorrelation(2a, a, 4, 5)
    @test 2 == p(0)
end
