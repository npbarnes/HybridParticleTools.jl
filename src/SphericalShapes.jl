module SphericalShapes
export SphericalShape, SphericalPolygon, Rotated, area, contains

using StaticArrays
using LinearAlgebra

function circular_index(a, i)
    N = length(a)
    if i<1
        while i < 1
            i += N
        end
    elseif i>N
        while i > N
            i -= N
        end
    end
    return a[i]
end

abstract type SphericalShape end
function area(::SphericalShape) end
function contains(::SphericalShape, P) end
(s::SphericalShape)(P) = contains(s, P)
# Broadcast as a scalar
Base.broadcastable(ss::SphericalShape) = Ref(ss)

struct FullSphere{T<:Number} <: SphericalShape end
FullSphere() = FullSphere{Float64}()
area(::FullSphere) = 4T(π)
contains(::FullSphere, P) = true


struct NoSphere{T<:Number} <: SphericalShape end
NoSphere() = NoSphere{Float64}()
area(::NoSphere) = zero(T)
contains(::NoSphere) = false

"Abstract Spherical Polygons have edges that are great circles"
abstract type AbstractSPolygon <: SphericalShape end
function inside(::AbstractSPolygon) end
function vertices(::AbstractSPolygon) end
function edges(ss::AbstractSPolygon)
    v = vertices(ss)
    ((v[i], circular_index(v,i+1)) for i in eachindex(v))
end
function corners(ss::AbstractSPolygon)
    v = vertices(ss)
    ((circular_index(v,i-1), v[i], circular_index(v,i+1))
        for i in eachindex(v))
end

"""
    contains(sp::AbstractSPolygon, P)
Check if the point P is contained in the SphericalPolygon.
It works by considering the arc between the known inside point and P. If this
arc intersects edges of the polygon and even number of times (including zero),
then P is inside. If the arc intersects the edges of the polygon an odd number
of times, then P is outside.
"""
function contains(sp::AbstractSPolygon, P)
    P = P ./ norm(P)
    I = (intersection(A,B,P,inside(sp)) for (A,B) in edges(sp))
    return count(!isnothing, I) % 2 == 0
end

"""
    area(sp::SphericalPolygon)
Compute the area of a convex spherical polygon.
If sp is concave the result is meaningless.
"""
function area(sp::AbstractSPolygon)
    s = sum(x->angle(x...), corners(sp))
    s - π * (length(sp.vertices) - 2)
end

"""
A polygon on a sphere
    Assumptions:
        1. Vertices are listed in order so that there is an edge between v[n] and v[n+1] for each n
        and from v[end] to v[1].
        2. No edges cover exactly pi radians. i.e. vertices connected by an edge are not antipodal.
        3. Edges are always defined to be the shorter arc between vertices (less than pi radians).
        4. The polygon is simply connected and not self intersecting.
        5. The given interior point is not on an edge.
    Of course, if you would like your polygon to have an edge equal to or longer than pi radians
    just use an extra vertex to make it into two edges on the same great circle.

    Some of the methods defined for SphericalPolygons may make additional assumptions. Be sure
    to read their respective documentation.
"""
struct SphericalPolygon{T<:Number} <: AbstractSPolygon
    inside::SVector{3,T}
    vertices::Vector{SVector{3,T}}
    function SphericalPolygon(inside::AbstractVector{<:Number}, vertices::AbstractVector{<:AbstractVector{<:Number}})
        new{eltype(inside)}(inside./norm(inside), vertices./norm.(vertices))
    end
end
inside(sp::SphericalPolygon) = sp.inside
vertices(sp::SphericalPolygon) = sp.vertices

"""
Compute angle between the spherical arc AB and the spherical arc BC.
Method is based on the spherical law of cosines.
"""
function angle(A,B,C)
    cb = A⋅C
    cc = B⋅A
    ca = C⋅B

    sc = norm(A×B)
    sa = norm(B×C)

    return acos((cb - cc*ca)/(sc*sa))
end

"""
Computes the intesection point of spherical arcs AB and CD.
Returns nothing if they do not intersect.
"""
function intersection(A,B,C,D)
    ABX = A × B
    CDX = C × D
    T = ABX × CDX
    T = T ./ norm(T)

    s1 = (ABX × A) ⋅ T
    s2 = (B × ABX) ⋅ T
    s3 = (CDX × C) ⋅ T
    s4 = (D × CDX) ⋅ T

    if all( (s1,s2,s3,s4) .>= 0 )
        return T
    elseif all( (s1,s2,s3,s4) .<= 0 )
        return -T
    else
        return nothing
    end
end

struct Rotated{SS, M} <: AbstractSPolygon where {SS<:AbstractSPolygon, M<:AbstractMatrix}
    sp::SS
    m::M
end
Rotated(r::Rotated, m) = Rotated(r.sp, m*r.m)
area(r::Rotated) = area(r.sp)
inside(r::Rotated) = r.m * inside(r.sp)
vertices(r::Rotated) = [r.m * v for v in vertices(r.sp)]

SphericalPolygon(r::Rotated) = SphericalPolygon(inside(r), collect(vertices(r)))

end # module
