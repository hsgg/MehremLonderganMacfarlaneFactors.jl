module MehremLonderganMacfarlaneFactors


export mehremlonderganmacfarlanefactor3j, mehremlonderganmacfarlanefactor


using SpecialFunctions: loggamma
using WignerSymbols: wigner6j
using Jacobi: legendre



function mehremlonderganmacfarlanefactor(l1, l2, l3, k1, k2, k3)
    mlm3j = mehremlonderganmacfarlanefactor3j(l1, l2, l3, k1, k2, k3)
    w3j000 = wigner3j000(eltype(promote(k1, k2, k3)), l1, l2, l3)
    mlm = mlm3j / w3j000
    @debug "mlm" l1,l2,l3 k1,k2,k3 w3j000 mlm
    return mlm
end



function mehremlonderganmacfarlanefactor3j(l1, l2, l3, k1::T, k2::T, k3::T) where {T<:Real}
    delta = (k1^2 + k2^2 - k3^2) / (2 * k1 * k2)
    beta_of_delta = beta(delta)

    if beta_of_delta == 0
        # protect the Legendre polynomial
        return T(0)
    end

    mlm3j = T(0)

    for L=0:l3

        lmin = max(abs(l1 - l3 + L), abs(l2 - L))
        lmax = min(l1 + l3 - L, l2 + L)
        lsum = T(0)

        odd1 = isodd(l1 + l3 - L + lmin)
        odd2 = isodd(l2 + L + lmin)
        if odd1 != odd2
            continue
        elseif odd1 || odd2
            lmin += 1
        end

        for l=lmin:2:lmax
            lterm = T(2 * l + 1)
            lterm *= wigner3j000(T, l1, l3 - L, l)
            lterm *= wigner3j000(T, l2, L, l)
            lterm *= wigner6j(T, l1, l2, l3, L, l3 - L, l)
            lterm *= legendre(delta, l)
            lsum += lterm
            # @show l,lterm,lsum
        end

        mlm3jterm = sqrt(binom(T, 2 * l3, 2 * L)) * (k2 / k1)^L * lsum
        mlm3j += mlm3jterm
        # @show L,mlm3jterm,mlm3j
    end

    mlm3j *= π * beta_of_delta / (4 * k1 * k2 * k3)

    mlm3j *= sqrt(T(2 * l3 + 1)) * (k1 / k3)^l3

    mlm3j *= (-1)^(((l1 + l2 + l3) ÷ 2) + l3)

    @debug "mlm3j" l1,l2,l3 k1,k2,k3 delta mlm3j

    return mlm3j
end

mehremlonderganmacfarlanefactor3j(l1, l2, l3, k1, k2, k3) = begin
    mehremlonderganmacfarlanefactor3j(l1, l2, l3, promote(k1, k2, k3)...)
end



function beta(delta::T) where {T<:Real}
    if abs(delta) > 1
        return T(0)
    elseif abs(delta) == 1
        return T(1//2)
    else
        return T(1)
    end
end


function binom(::Type{T}, n, k) where {T<:Real}
    return exp(loggamma(T(1 + n))
               - loggamma(T(1 + k))
               - loggamma(T(1 + n - k)))
end

binom(n, k) = binom(Float64, n, k)


function wigner3j000(::Type{T}, l, l′, L) where {T<:Real}
    (abs(l-l′) <= L <= l+l′) || return T(0)
    J = l + l′ + L
    iseven(J) || return T(0)
    wig3j = (-1)^(J÷2) * exp(loggamma(T(1+J-2l))/2 + loggamma(T(1+J-2l′))/2
                             + loggamma(T(1+J-2L))/2 - loggamma(T(1+J+1))/2
                             + loggamma(T(1+J÷2))
                             - loggamma(T(1+J÷2-l)) - loggamma(T(1+J÷2-l′))
                             - loggamma(T(1+J÷2-L)))
    return wig3j
end


wigner3j000(l, l′, L) = wigner3j000(Float64, l, l′, L)



end
