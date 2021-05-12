module HybridGrids

export loadfields, loadvector, loadscalar, VectorField, Covariant, Contravariant, curl, loadB, loadE, Covariant_fromfunction

using PyCall
using Interpolations
using Interpolations: extrapolate, Flat, Periodic
using StaticArrays
using Unitful
using LinearAlgebra
using PhysicalConstants.CODATA2018: m_p, e, μ_0

using ..ParameterSets

const hr = PyNULL()

function __init__()
    copy!(hr,pyimport("HybridReader2").HybridReader2)
end

struct Grid{T}
    nodes::NTuple{3, Vector{T}}
end
Grid(x,y,z) = Grid((x,y,z))

function _ave(q)
    ret = similar(q)
    @. ret[1:end-1] = @views 0.5 * (q[1:end-1] + q[1+1:end])
    ret[end] = q[end] + (q[end] - ret[end-1])
    return ret
end
function dualgrid(g::Grid)
    (_ave(g.nodes[1]), _ave(g.nodes[2]), _ave(g.nodes[3]))
end

abstract type VectorField{G,D,I} end

struct Covariant{G,D,I} <: VectorField{G,D,I}
    maingrid::Grid{G}
    data::Array{D,4}
    xinterp::I
    yinterp::I
    zinterp::I
    function Covariant(maingrid::Grid{G}, data::AbstractArray{D,4}) where {G,D}
        mx, my, mz = maingrid.nodes
        dx, dy, dz = dualgrid(maingrid)
        @views begin
            xinterp = extrapolate(interpolate((mx, dy, dz), data[:,:,:,1], Gridded(Linear())), Flat())
            yinterp = extrapolate(interpolate((dx, my, dz), data[:,:,:,2], Gridded(Linear())), Flat())
            zinterp = extrapolate(interpolate((dx, dy, mz), data[:,:,:,3], Gridded(Linear())), Flat())
        end
        V = typeof(xinterp)
        new{G,D,V}(maingrid, data, xinterp, yinterp, zinterp)
    end
end
function Base.:+(a::Covariant, b::Covariant)
    @assert a.maingrid == b.maingrid
    Covariant(a.maingrid, a.data .+ b.data)
end

Covariant_fromfunction(maingrid, f) = Covariant_fromfunction(maingrid, (x,y,z)->f(x,y,z)[1], (x,y,z)->f(x,y,z)[2], (x,y,z)->f(x,y,z)[3])
function Covariant_fromfunction(maingrid, fx, fy, fz)
    mx, my, mz = maingrid.nodes
    dx, dy, dz = dualgrid(maingrid)
    data = Array{Float64,4}(undef, length(mx), length(my), length(mz), 3)
    for (i,x) in pairs(mx), (j,y) in pairs(dy), (k,z) in pairs(dz)
        data[i,j,k,1] = fx(x,y,z)
    end
    for (i,x) in pairs(dx), (j,y) in pairs(my), (k,z) in pairs(dz)
        data[i,j,k,2] = fy(x,y,z)
    end
    for (i,x) in pairs(dx), (j,y) in pairs(dy), (k,z) in pairs(mz)
        data[i,j,k,3] = fz(x,y,z)
    end
    Covariant(maingrid, data)
end


struct Contravariant{G,D,I} <: VectorField{G,D,I}
    maingrid::Grid{G}
    data::Array{D,4}
    xinterp::I
    yinterp::I
    zinterp::I
    function Contravariant(maingrid::Grid{G}, data::AbstractArray{D,4}) where {G,D}
        mx, my, mz = maingrid.nodes
        dx, dy, dz = dualgrid(maingrid)
        @views begin
            xinterp = extrapolate(interpolate((dx, my, mz), data[:,:,:,1], Gridded(Linear())), Flat())
            yinterp = extrapolate(interpolate((mx, dy, mz), data[:,:,:,2], Gridded(Linear())), Flat())
            zinterp = extrapolate(interpolate((mx, my, dz), data[:,:,:,3], Gridded(Linear())), Flat())
        end
        V = typeof(xinterp)
        new{G,D,V}(maingrid, data, xinterp, yinterp, zinterp)
    end
end

function (f::VectorField)(x, y, z)
    SA[f.xinterp(x, y, z), f.yinterp(x, y, z), f.zinterp(x, y, z)]
end

function spacings(g)
    dx_cells = similar(g.nodes[1], Union{Missing, eltype(g.nodes[1])})
    dy_cells = similar(g.nodes[2], Union{Missing, eltype(g.nodes[2])})
    dz_cells = similar(g.nodes[3], Union{Missing, eltype(g.nodes[3])})
    ret = (dx_cells, dy_cells, dz_cells)
    for (i,q) in pairs(g.nodes)
        @. ret[i][2:end-1] = @views ((q[3:end] + q[2:end-1])/2) - ((q[2:end-1] + q[1:end-2])/2)
        ret[i][1] = missing
        ret[i][end] = missing
    end
    ret
end

function curl(vf::Covariant)
    mx, my, mz = vf.maingrid.nodes
    dx, dy, dz = dualgrid(vf.maingrid)
    c = similar(vf.data, Union{Missing, eltype(vf.data)})
    c .= missing
    dx,dy,dz = spacings(vf.maingrid)
    for i in 2:size(c,1), j in 2:size(c,2), k in 2:size(c,3)
        c[i,j,k,1] = (vf.data[i,j,k,3] - vf.data[i,j-1,k,3])/dy[j] - (vf.data[i,j,k,2] - vf.data[i,j,k-1,2])/dz[k]
        c[i,j,k,2] = (vf.data[i,j,k,1] - vf.data[i,j,k-1,1])/dz[k] - (vf.data[i,j,k,3] - vf.data[i-1,j,k,3])/dx[i]
        c[i,j,k,3] = (vf.data[i,j,k,2] - vf.data[i-1,j,k,2])/dx[i] - (vf.data[i,j,k,1] - vf.data[i,j-1,k,1])/dy[j]
    end
    Contravariant(vf.maingrid, c)
end

const α = e^2*μ_0/m_p
Ep(Bp, cBp, np, up, α) = -(up - cBp/(α*np)) × Bp
Ep(Bp::AbstractVector{<:Quantity}, cBp::AbstractVector{<:Quantity}, np::AbstractVector{<:Quantity}, up::AbstractVector{<:Quantity}) = Ep(Bp,cBp,np,up,α)
Ep(Bp, cBp, np, up) = Ep(Bp, cBp, np, up, ustrip(u"C*T*km*s/kg", α))

struct E_Field{T,U,V,Q}
    B::T
    cB::U
    u_i::V
    n_i::Q
end
(F::E_Field)(x,y,z) = Ep(F.B(x,y,z), F.cB(x,y,z), F.n_i(x,y,z), F.u_i(x,y,z))

function loadvector(prefix, name, step=-1)
    h = hr(prefix, name)
    _, raw = h.get_timestep(step)
    data = similar(raw, SVector{3,Float64}, size(raw)[1:end-1])
    for i in axes(raw,1), j in axes(raw,2), k in axes(raw,3)
        data[i,j,k] = SVector{3,Float64}(raw[i,j,k,:])
    end
    para  = ParameterSet(prefix)
    nodes = Tuple(convert(Array{Float64}, gp) for gp in para.grid_points)

    interpolate(nodes, data, Gridded(Linear()))
end

function loadscalar(prefix, name, step=-1)
    _, data = hr(prefix, name).get_timestep(step)
    para = ParameterSet(prefix)
    nodes = Tuple(convert(Array{Float64}, gp) for gp in para.grid_points)
    interpolate(nodes, data, Gridded(Linear()))
end

# There was a problem writing bt data for 2020-Jan-23/pluto-2.
# Timestep 27 is the last good step.
function loadB(prefix, step=-1)
    _, B_data = hr(prefix, "bt").get_timestep(step)
    para = ParameterSet(prefix)
    mion = para.ion_amu
    B_data = ustrip.(u"T", B_data.*u"s^-1".*(mion*m_p/e))
    nodes = Tuple(convert(Array{Float64}, gp) for gp in para.grid_points)
    maingrid = Grid(nodes)
    B = Covariant(maingrid, B_data)
    return B
end
function loadE(prefix, step=-1, B=loadB(prefix, step))
    cB = curl(B)
    u_i = loadvector(prefix, "up", step)
    n_i = loadscalar(prefix, "np", step)
    E_Field(B,cB,u_i,n_i)
end
function loadfields(prefix, step=-1)
    B = loadB(prefix, step)
    E = loadE(prefix, step, B)
    E,B
end

end # module
