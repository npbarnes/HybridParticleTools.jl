module SphericalPolygons
export SphericalPolygon, area, contains

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
struct SphericalPolygon{T<:Number}
    inside::SVector{3,T}
    vertices::Vector{SVector{3,T}}
    function SphericalPolygon(inside::AbstractVector{<:Number}, vertices::AbstractVector{<:AbstractVector{<:Number}})
        new{eltype(inside)}(inside./norm(inside), vertices./norm.(vertices))
    end
end
Base.broadcastable(sp::SphericalPolygon) = Ref(sp)

edges(sp::SphericalPolygon) = ((sp.vertices[i], circular_index(sp.vertices,i+1)) for i in 1:length(sp.vertices))
corners(sp::SphericalPolygon) = ((circular_index(sp.vertices,i-1), sp.vertices[i], circular_index(sp.vertices,i+1))
                                        for i in 1:length(sp.vertices))

function contains(sp::SphericalPolygon, P)
    P = P ./ norm(P)
    I = (intersection(A,B,P,sp.inside) for (A,B) in edges(sp))
    return count(!isnothing, I) % 2 == 0
end

"""
    area(sp::SphericalPolygon)
Compute the area of a convex spherical polygon.
If sp is concave the result is meaningless.
"""
function area(sp::SphericalPolygon)
    s = sum(x->angle(x...), corners(sp))
    s - π * (length(sp.vertices) - 2)
end

"""
Compute the spherical angle ABC.
"""
function angle(A,B,C)
    cb = A⋅C
    cc = B⋅A
    ca = C⋅B

    sc = norm(A×B)
    sa = norm(B×C)

    return acos((cb - cc*ca)/(sc*sa))
end

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

end # module
