module Distributions
export Distribution, DistElement

using StaticArrays
using StructArrays
using Unitful
using ..Simulations

struct DistElement{U,V,W,X}
    v::SVector{3,U}
    m::V
    q::W
    n::X
end
Base.show(io::IO, e::DistElement) = print(io, "DistElement(v=SA[$(e.v[1]),$(e.v[2]),$(e.v[3])], m=$(e.m), q=$(e.q), n=$(e.n))")
const Distribution = StructArray{DistElement{U,V,W,X}} where {U,V,W,X}
Base.show(io::IO, ::Type{<:Distribution{U,V,W,X}}) where {U,V,W,X} = print(io, "Distribution{$U,$V,$W,$X}")
Base.show(io::IO, d::Distribution) = print(io, "$(length(d))-element Distribution...")
function Base.show(io::IO, ::MIME"text/plain", d::Distribution)
    print(io, "$(length(d))-element Distribution:")
    if length(d) <= 10
        for e in d
            print(io, "\n   ")
            print(io, e)
        end
    else
        for e in @view d[1:3]
            print(io, "\n   ")
            print(io, e)
        end
        print(io,"\n   ⋮")
        for e in @view d[end-2:end]
            print(io, "\n   ")
            print(io, e)
        end
    end
end
function Distribution(v::AbstractArray,m::AbstractArray,q::AbstractArray,n::AbstractArray)
    Distribution{
        eltype(eltype(v)),
        eltype(m),
        eltype(q),
        eltype(n)
    }((v,m,q,n))
end

function weight(x,xp,dx)
    (abs(1 - abs(x[1]-xp[1])/dx) * abs(1 - abs(x[2]-xp[2])/dx) * abs(1 - abs(x[3]-xp[3])/dx)) / dx^3
end

"Assume units are kilometers, unless given"
Distribution(s::Simulation, x::AbstractVector{<:Number}, dx=s.dx) = Distribution(s, x*u"km", dx*u"km")
function Distribution(s::Simulation, x::AbstractVector{<:Quantity}, dx=s.dx*u"km")
    t = touching(s, x, dx)
    n = similar(s.N[t])*u"km^-3"
    for (i,(xp,N)) in enumerate(zip(s.x[t],s.N[t]))
        n[i] = N*weight(x,xp,dx)
    end
    Distribution(s.v[t], s.m[t], s.q[t], n)
end
Base.:*(cmat::AbstractMatrix, d::Distribution) = Distribution(
    [cmat * v for v in d.vs],
    d.ms,
    d.qs,
    d.ns
)

end # module
