module Sensors
export density, bulkvelocity, pressuretensor, pressure, thermalenergy, flux,
    fovfilter, spectrum_data, espec

using StaticArrays
using LinearAlgebra
using StatsBase
using Unitful
using PyCall
using StatsBase

using ..Utility
using ..Simulations
using ..Distributions
using ..Spacecraft
using ..SphericalShapes

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
"Also called kinetic temperature"
thermalenergy(f::Distribution) = pressure(f)/density(f)

flux(e::DistElement) = norm(e.n .* e.v)
flux(d::Distribution{U,V,W,X}) where {U,V,W,X} = isempty(d) ? zero(U)*zero(X) : sum(flux, d)

"Compute E/q spectrum for the given distribution"
function energypercharge(e::DistElement)
    (1/2 .* e.m .* norm(e.v).^2) ./ e.q
end

const pepssifov = let
    st = pyimport("spice_tools")
    shape, frame, bsight, n, bounds = st.sp.getfov(-98402, 25, 255, 255)
    SphericalPolygon(bsight, listofvectors(bounds))
end

function fovfilter(d::Distribution, fov=SphericalShapes.FullSphere())
    filter(e->contains(fov, -e.v), d)
end

function espec(s::Simulation, t::Trajectory, J=1.5u"cm^2*sr*s")
    full_dists = Distribution.(s, t.pos)
    rfd = @view full_dists[length.(full_dists) .!= 0]
    dists = [fovfilter(d, Rotated(pepssifov, cmat)) for (d,cmat) in zip(full_dists, t.cmat)]
    counts = ustrip.(Unitful.NoUnits, J.*flux.(dists))
    Eq = [energypercharge.(d) for d in dists]
    (counts=counts, Eq=Eq)
end

function spectrum_data(d::Distribution, J=1.5u"cm^2*sr*s")
    (count = ustrip.(Unitful.NoUnits, J.*flux.(d)),
     Eq = ustrip.(u"J/C", energypercharge.(d))
    )
end

function prepare_espec(counts, Eq, bins=geomspace(100,10000,200))
    histlist = [fit(Histogram, ustrip.(Unitful.NoUnits, e), weights(c), bins) for (c,e) in zip(counts, Eq)]
end

end # module
