module HybridGrids

export loadfields, loadvector, loadfields2, loadscalar, VectorField

using PyCall
using Interpolations
using Interpolations: extrapolate, Flat, Periodic
using StaticArrays
using Unitful
using LinearAlgebra
using PhysicalConstants.CODATA2018: m_p, e

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
    ret[1:end-1] = 0.5*(q[1:end-1] + q[1+1:end])
    ret[end] = q[end] + (q[end] - ret[end-1])
    return ret
end
function _dualgrid(g::Grid)
    (_ave(g.nodes[1]), _ave(g.nodes[2]), _ave(g.nodes[3]))
end

struct VectorField{T,U,V}
    xgrid::Grid{T}
    ygrid::Grid{T}
    zgrid::Grid{T}
    F::Array{U, 4}
    xinterp::V
    yinterp::V
    zinterp::V
    function VectorField(gx::Grid{T},gy,gz, F::AbstractArray{U,4}) where {T,U}
        @views begin
            xinterp = extrapolate(interpolate(gx.nodes, F[:,:,:,1], Gridded(Linear())), Flat())
            yinterp = extrapolate(interpolate(gy.nodes, F[:,:,:,2], Gridded(Linear())), Flat())
            zinterp = extrapolate(interpolate(gz.nodes, F[:,:,:,3], Gridded(Linear())), Flat())
        end
        V = typeof(xinterp)
        new{T,U,V}(gx, gy, gz, F, xinterp, yinterp, zinterp)
    end
end

function covariant(maingrid::Grid, field)
    mx, my, mz = maingrid.nodes
    dx, dy, dz = _dualgrid(maingrid)
    VectorField(
        Grid(mx, dy, dz),
        Grid(dx, my, dz),
        Grid(dx, dy, mz),
        field
    )
end
function contravariant(maingrid::Grid{T}, field) where T
    mx, my, mz = maingrid.nodes
    dx, dy, dz = _dualgrid(maingrid)
    VectorField(
        Grid(dx, my, mz),
        Grid(mx, dy, mz),
        Grid(mx, my, dz),
        field
    )
end

function (f::VectorField)(x, y, z)
    SA[f.xinterp(x, y, z), f.yinterp(x, y, z), f.zinterp(x, y, z)]
end

"""Loads a vector field from the simulation and interpolates it.
This method ignores the Yee grid and assumes that the vector components
are colocated.
"""
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

function loadfields(prefix, step=-1)
    _, E_data = hr(prefix, "E").get_timestep(step)
    _, B_data = hr(prefix, "bt").get_timestep(step)
    para = ParameterSet(prefix)
    mion = para.ion_amu
    E_data = ustrip.(u"V/m", E_data.*u"km/s^2".*(mion*m_p/e))
    B_data = ustrip.(u"T", B_data.*u"s^-1".*(mion*m_p/e))
    nodes = Tuple(convert(Array{Float64}, gp) for gp in para.grid_points)
    maingrid = Grid(nodes)
    E = contravariant(maingrid, E_data)
    B = covariant(maingrid, B_data)
    return E,B
end

end # module
