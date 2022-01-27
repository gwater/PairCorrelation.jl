module PairCorrelation

using StaticArrays
using LinearAlgebra
using FFTW
using SpecialFunctions

export PairAutoCorrelation

function rfftkspace2(nx, ny, lx, ly)
    return (SVector(x, y)
        for x in rfftfreq(nx, 2pi * nx / lx),
            y in fftfreq(ny, 2pi * ny / ly))
end

struct PairAutoCorrelation{T, S}
    ks::T
    ak::S
    n::Int
    function PairAutoCorrelation(a, lx, ly)
        ak = rfft(a)
        ks = rfftkspace2(size(a)..., lx, ly)
        return new{typeof(ks), typeof(ak)}(ks, ak, length(a))
    end
end

(p::PairAutoCorrelation)(r) =
    sum(
        abs2.(p.ak) .*
        besselj0.(r .* norm.(p.ks)) .*
        ifelse.(iszero.(getindex.(p.ks, 1)), 1, 2)
    ) / abs2(p.n)

end # module
