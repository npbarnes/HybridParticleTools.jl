module Sensors
export density, bulkvelocity, pressuretensor, pressure, thermalenergy, flux,
       energy, energies, energypercharge, fluxes

using LinearAlgebra
using Unitful
using ..SphericalShapes
using ..Distributions

# Get moments from a distribution
density(f::Distribution) = sum(f.n)
bulkvelocity(f::Distribution) = sum(f.v .* f.n)/density(f)
function pressuretensor(f::Distribution)
    uf = bulkvelocity(f)
    sum(zip(f.v, f.n, f.m)) do (v,n,m)
        w = v-uf
        (m * w)*(w' * n)
    end
end
pressure(f::Distribution) = tr(pressuretensor(f))/3
"Also called kinetic temperature or kT = P/n"
thermalenergy(f::Distribution) = pressure(f)/density(f)

flux(e::DistElement) = norm(e.n .* e.v)
flux(d::Distribution{U,V,W,X}) where {U,V,W,X} = isempty(d) ? zero(U)*zero(X) : sum(flux, d)
fluxes(d::Distribution) = flux.(d)

energy(e::DistElement) = (1/2 * e.m * norm(e.v)^2)
energies(d::Distribution) = energy.(d)
energypercharge(e::DistElement) = energy(e) / e.q

function _differential_intensity(d::Distribution, Elo, Ehi, Ω)
    ΔE = Ehi - Elo
    dd = filter(e->Elo<=energy(e)<=Ehi, d)
    flux(dd)/(Ω*ΔE)
end
function differential_intensity(fov::SphericalShape, bin_edges, d::Distribution)
    Ω = area(fov)u"sr"
    dd = filter(fov, d)

    first = _differential_intensity(dd, bin_edges[1], bin_edges[2], Ω)

    ret = Vector{typeof(first)}(undef, length(bin_edges)-1)
    ret[1] = first
    for (i,(Elo,Ehi)) in enumerate(zip(@view(bin_edges[2:end]), @view(bin_edges[3:end])))
        @inbounds ret[i+1] = _differential_intensity(dd, Elo, Ehi, Ω)
    end
    return ret
end
end # module
