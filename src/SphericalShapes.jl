module SphericalShapes
export  SphericalShape, SphericalPolygon, SPolygon, SCircle, Rotated, area, contains, vertices, edges, corners, inside,
        r_transform, latlon, SREAG_grid, LLAlignedRectangle

using PyCall
using StaticArrays
using LinearAlgebra

using ..Utility

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
area(::FullSphere{T}) where T = 4T(π)
contains(::FullSphere, P) = true


struct NoSphere{T<:Number} <: SphericalShape end
NoSphere() = NoSphere{Float64}()
area(::NoSphere{T}) where T = zero(T)
contains(::NoSphere, P) = false

"Abstract Spherical Polygons have edges that are great circles"
abstract type SphericalPolygon <: SphericalShape end
function inside(::SphericalPolygon) end
function vertices(::SphericalPolygon) end
function edges(ss::SphericalPolygon)
    v = vertices(ss)
    ((v[i], circular_index(v,i+1)) for i in eachindex(v))
end
function corners(ss::SphericalPolygon)
    v = vertices(ss)
    ((circular_index(v,i-1), v[i], circular_index(v,i+1))
        for i in eachindex(v))
end

"""
    contains(sp::SphericalPolygon, P)
Check if the point P is contained in the SPolygon.
It works by considering the arc between the known inside point and P. If this
arc intersects edges of the polygon and even number of times (including zero),
then P is inside. If the arc intersects the edges of the polygon an odd number
of times, then P is outside.
"""
function contains(sp::SphericalPolygon, P)
    P = P ./ norm(P)
    I = (intersection(A,B,P,inside(sp)) for (A,B) in edges(sp))
    return count(!isnothing, I) % 2 == 0
end

"""
    area(sp::SPolygon)
Compute the area of a convex spherical polygon.
If sp is concave the result is meaningless.
"""
function area(sp::SphericalPolygon)
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

    Some of the methods defined for SPolygons may make additional assumptions. Be sure
    to read their respective documentation.
"""
struct SPolygon{T<:Number} <: SphericalPolygon
    inside::SVector{3,T}
    vertices::Vector{SVector{3,T}}
    function SPolygon(inside::AbstractVector{<:Number}, vertices::AbstractVector{<:AbstractVector{<:Number}})
        new{eltype(inside)}(inside./norm(inside), vertices./norm.(vertices))
    end
end
SPolygon(o::PyObject) = SPolygon(o.inside, listofvectors(o.vertices))
inside(sp::SPolygon) = sp.inside
vertices(sp::SPolygon) = sp.vertices

Utility.r_transform(s::SPolygon) = SPolygon(r_transform(s.inside), r_transform.(s.vertices))

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

struct Rotated{SS, M} <: SphericalPolygon
    sp::SS
    m::M
end
Rotated(r::Rotated, m) = Rotated(r.sp, m*r.m)
area(r::Rotated) = area(r.sp)
inside(r::Rotated) = r.m * inside(r.sp)
vertices(r::Rotated) = [r.m * v for v in vertices(r.sp)]

SPolygon(r::Rotated) = SPolygon(inside(r), collect(vertices(r)))

struct SCircle{T,U} <: SphericalShape
    center::SVector{3,T}
    angle::U # Angle from center to the edge of the circle
    cosangle::U
    function SCircle(center,angle)
        @assert angle >= 0 && angle <= π
        normalized = center./norm(center)
        new{eltype(normalized),typeof(angle)}(normalized,angle,cos(angle))
    end
end
inside(c::SCircle) = c.center
area(c::SCircle) = 2π*(1-c.cosangle)
function contains(c::SCircle, P)
    P = P ./ norm(P)
    cosine = c.center ⋅ P
    return cosine >= c.cosangle
end

function latlon(P)
    lat = atan(P[3], sqrt(P[1]^2 + P[2]^2))
    lon = atan(P[2],P[1])
    lat,lon
end
abstract type LatitudeLongitudePolygon <: SphericalShape end
struct LLAlignedRectangle{T<:Number} <: LatitudeLongitudePolygon
    lat::Tuple{T,T}
    lon::Tuple{T,T}
    function LLAlignedRectangle(lat::Tuple{T,T},lon::Tuple{T,T}) where T<:Number
        @assert -π/2 <= lat[1] <= lat[2] <= π/2
        @assert -π <= lon[1] <= π
        @assert -π <= lon[2] <= π
        new{T}(lat, lon)
    end
end
function LLAlignedRectangle(lat, lon)
    @assert length(lat) == 2
    @assert length(lon) == 2
    LLAlignedRectangle(lat[1],lat[2],lon[1],lon[2])
end
LLAlignedRectangle(lat1,lat2,lon1,lon2) = LLAlignedRectangle(promote(lat1, lat2, lon1, lon2)...)
LLAlignedRectangle(lat1::T, lat2::T, lon1::T, lon2::T) where T = LLAlignedRectangle((lat1,lat2), (lon1,lon2))
vertices(r::LLAlignedRectangle) = [
    r.lon[1] r.lat[1];
    r.lon[2] r.lat[1];
    r.lon[2] r.lat[2];
    r.lon[1] r.lat[2];
    r.lon[1] r.lat[1]
]

_lat_okay(r,lat) = r.lat[1] < lat < r.lat[2]
function _lon_okay(r,lon)
    if r.lon[1] <= r.lon[2]
        lon_okay = r.lon[1] <= lon < r.lon[2]
    else
        lon_okay = lon >= r.lon[1] || lon < r.lon[2]
    end
    lon_okay
end

function contains(r::LLAlignedRectangle, P)
    lat,lon = latlon(P)
    return _lat_okay(r,lat) && _lon_okay(r,lon)
end
function area(r::LLAlignedRectangle)
    lon_diff = abs(r.lon[2] - r.lon[1])
    if r.lon[2] < r.lon[1]
        lon_diff = 2π - lon_diff
    end
    abs(sin(r.lat[1]) - sin(r.lat[2])) * lon_diff
end

_centers_fromedges(b) = [(b1+b2)/2 for (b1,b2) in zip(@view(b[begin:end-1]), @view(b[begin+1:end]))]
"""Generate an equal area grid on the sphere using the SREAG method described in the paper:
Malkin, Zinovy "A New Equal-area Isolatitudinal Grid on a Spherical Surface" The Astronomical Journal (2019)
Two values are returned, `l,b`. `l` is an array of ranges where l[i] gives the longitudes of the cell
boundaries on the ith row of cells. Longitude is taken to range from -π to π. `b` is the array of latitudes
of cell boundaries. Latitude is taken to be in the range -π/2 to π/2.
"""
function _sphericalgrid_nodes(N_ring)
    if N_ring % 2 == 1
        error("Number of rings must be even. Got $(N_ring).")
    end
    dB = π/N_ring
    b_init = -π/2:dB:π/2
    dL_init = dB .* sec.(_centers_fromedges(b_init))
    N_cell_row = round.(Int, 2π ./ dL_init)
    N_cell = sum(N_cell_row)

    A  = 4π/N_cell
    dL = 2π ./ N_cell_row

    b = similar(b_init)
    b[1] = -π/2
    for i in 1:(N_ring ÷ 2)
        b[i+1] = asin(A/dL[i] + sin(b[i]))
    end
    b[N_ring÷2+2:end] = -reverse(b[begin:N_ring÷2])

    l = [-π:dl:π for dl in dL]
    l,b
end
function _gridsquares_latlon(l,b)
    r = LLAlignedRectangle{eltype(b)}[]
    for (i, row_l) in pairs(l)
        for j in 1:(length(row_l)-1)
            push!(r, LLAlignedRectangle((b[i], b[i+1]), (row_l[j], row_l[j+1])) )
        end
    end
    r
end

function SREAG_grid(N)
    l,b = _sphericalgrid_nodes(N)
    _gridsquares_latlon(l,b)
end

end # module
