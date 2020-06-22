module Sensors
export density, bulkvelocity, pressuretensor, pressure, thermalenergy, flux,
       energy, energies, energypercharge, fluxes

using LinearAlgebra
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

function differential_intensity(d, fov, Emin, Emax)
    dd = filter(fov, d)
    dd = filter(e->Emin <= energy(e) <= Emax, dd)
    return flux(dd)/(area(fov)*(Emax-Emin))
end
end # module
