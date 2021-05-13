using LinearAlgebra
using Unitful
using StaticArrays

struct Dipole{T}
    m::SVector{3,T}
    function Dipole(m::SVector{3})
        new{eltype(m)}(m)
    end
end

function _dipole(m⃗, r⃗, μ)
    r = norm(r⃗)
    r̂ = r⃗/r
    μ/(4π*r^3) * (3(m⃗⋅r̂)r̂ - m⃗)
end
(d::Dipole)(x,y,z) = d(SA[x,y,z])
(d::Dipole{<:Quantity})(r⃗::AbstractVector{<:Quantity}) = _dipole(d.m, r⃗, μ_0)
(d::Dipole{<:Number})(r⃗::AbstractVector{<:Number}) = _dipole(d.m, r⃗, ustrip(u"kg*km*s^-2*A^-2", μ_0)) # Hybrid code units
