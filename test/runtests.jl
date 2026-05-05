using MehremLonderganMacfarlaneFactors
using Test
using Jacobi: legendre

@testset "MehremLonderganMacfarlaneFactors.jl" begin

    #ENV["JULIA_DEBUG"] = MehremLonderganMacfarlaneFactors

    # l1 = l2 = l3 = 0
    krange = 0.1:0.1:0.3

    for k1 in krange, k2 in krange, k3 in krange

        delta = (k1^2 + k2^2 - k3^2) / (2 * k1 * k2)

        if MehremLonderganMacfarlaneFactors.beta(delta) == 1//2
            @test mehremlonderganmacfarlanefactor3j(0, 0, 0, k1, k2, k3) ≈ π / (8 * k1 * k2 * k3)
        elseif abs(k1 - k2) < k3 < k1 + k2
            @test mehremlonderganmacfarlanefactor3j(0, 0, 0, k1, k2, k3) ≈ π / (4 * k1 * k2 * k3)
        else
            @test mehremlonderganmacfarlanefactor3j(0, 0, 0, k1, k2, k3) ≈ 0  atol=eps(π / (8 * k1 * k2 * k3))
        end

        # I hate floating point with discontinuities:
        if abs(k3 - (k1 + k2)) < eps(1.0) || abs(k3 - abs(k1 - k2)) < eps(1.0)
            broken = false
            if delta == -0.9999999999999994
                broken=true
            end
            @test MehremLonderganMacfarlaneFactors.beta(delta) == 1//2  broken=broken
        elseif abs(k1 - k2) < k3 < k1 + k2
            @test MehremLonderganMacfarlaneFactors.beta(delta) == 1
        else
            @test MehremLonderganMacfarlaneFactors.beta(delta) == 0
        end
    end


    # scaling formula
    krange = 0.1:0.1:0.3
    lrange = 0:5
    for k1 in krange, k2 in krange, k3 in krange
        for l1 in lrange, l2 in lrange, l3 in abs(l1 - l2):1:(l1 + l2)
            isodd(l1 + l2 + l3) && continue
            if abs(k3 - (k1 + k2)) < eps(1.0) || abs(k3 - abs(k1 - k2)) < eps(1.0)
                # I hate floating point.
                continue
            end

            @test mehremlonderganmacfarlanefactor3j(l1, l2, l3, k1, k2, k3) ≈ k1^-3 * mehremlonderganmacfarlanefactor3j(l1, l2, l3, 1, k2/k1, k3/k1)
        end
    end


    # Table 1 in Mehrem etal. 1991
    @test mehremlonderganmacfarlanefactor(0, 0, 0, 1.0, 2.0, 1.5) ≈  0.2617993877992  rtol=1e-12 atol=abs(0.2617993877991 - 0.2617993877992)
    @test mehremlonderganmacfarlanefactor(0, 1, 1, 1.0, 2.0, 1.5) ≈  0.2290744643243  rtol=1e-12 atol=abs(0.2290744643243 - 0.2290744643243)
    @test mehremlonderganmacfarlanefactor(1, 1, 0, 1.0, 2.0, 1.5) ≈  0.1799870791119  rtol=1e-12 atol=abs(0.1799870791119 - 0.1799870791119)
    @test mehremlonderganmacfarlanefactor(3, 2, 1, 1.0, 2.0, 1.5) ≈  0.1128754196419  rtol=1e-12 atol=abs(0.1128754196419 - 0.1128754196419)
    @test mehremlonderganmacfarlanefactor(4, 4, 4, 1.0, 2.0, 1.5) ≈ -0.0456849264425  rtol=1e-12 atol=abs(-0.0456849264425 - -0.0456849264425)
    @test mehremlonderganmacfarlanefactor(0, 0, 0, 1.0, 5.0, 5.5) ≈  0.0285599332145  rtol=1e-11 atol=abs(0.0285599332145 - 0.0285599332145)  # insufficient comparison
    @test mehremlonderganmacfarlanefactor(4, 4, 4, 1.0, 5.0, 5.5) ≈ -0.0080238752369  rtol=1e-12 atol=abs(-0.0080238752371 - -0.0080238752369)
    @test mehremlonderganmacfarlanefactor(0, 0, 0, 1.0, 1.05, 0.06) ≈ 12.46663751425  rtol=1e-12 atol=abs(12.46663751425 - 12.46663751425)
    @test mehremlonderganmacfarlanefactor(4, 4, 4, 1.0, 1.05, 0.06) ≈ -1.511653808097  rtol=1e-12 atol=3*abs(-1.511653806929 - -1.511653808097)  # catastrophic cancellation
    @test mehremlonderganmacfarlanefactor(0, 0, 0, 1.0, 1.02, 0.03) ≈  25.6666066469758  rtol=1e-12 atol=abs(25.6666066469754 - 25.6666066469758)
    @test mehremlonderganmacfarlanefactor(4, 4, 4, 1.0, 1.02, 0.03) ≈ -10.9755226591813  rtol=1e-12 atol=abs(-10.9755226231663 - -10.9755226591813)


    # BigFloat
    @test mehremlonderganmacfarlanefactor(0, 0, 0, 1.0, 1.05, big(0.06)) ≈ 12.46663751425  rtol=1e-12
    @test mehremlonderganmacfarlanefactor(4, 4, 4, 1.0, 1.05, big(0.06)) ≈ -1.511653808097  rtol=1e-12
    @test mehremlonderganmacfarlanefactor(0, 0, 0, 1.0, 1.02, big(0.03)) ≈  25.6666066469758  rtol=1e-12
    @test mehremlonderganmacfarlanefactor(4, 4, 4, 1.0, 1.02, big(0.03)) ≈ -10.9755226591813  rtol=1e-12


    # Special case: (ℓ₁, ℓ₂, ℓ₃) = (ℓ, ℓ, 0)
    krange = 0.1:0.1:0.3
    for k1 in krange, k2 in krange, k3 in krange, ell=0:1000
        delta = (k1^2 + k2^2 - k3^2) / (2 * k1 * k2)
        if MehremLonderganMacfarlaneFactors.beta(delta) == 0
            simplified = 0.0  # protect the legendre polynomial
        else
            simplified = π * MehremLonderganMacfarlaneFactors.beta(delta) / (4 * k1 * k2 * k3) * legendre(delta, ell)
        end
        @test mehremlonderganmacfarlanefactor(ell, ell, 0, k1, k2, k3) ≈ simplified
    end


    # types
    @inferred Float32 mehremlonderganmacfarlanefactor3j(0, 0, 0, 1, 2, 1.5f0)
    @inferred Float32 mehremlonderganmacfarlanefactor(0, 0, 0, 1, 2, 1.5f0)
    @inferred BigFloat mehremlonderganmacfarlanefactor3j(0, 0, 0, 1, big(2), 1.5f0)
    @inferred BigFloat mehremlonderganmacfarlanefactor(0, 0, 0, 1, big(2), 1.5f0)
end
