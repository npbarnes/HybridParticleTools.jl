module HybridGrids

export loadfields, loadvector, loadfields2, loadscalar, VectorField

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
struct YeeGrid{T}
    maingrid::Grid{T}
    dualgrid::Grid{T}
    function YeeGrid(maingrid::Grid{T}) where T
        new{T}(maingrid, _dualgrid(maingrid))
    end
end
function _ave(q)
    ret = similar(q)
    @. ret[1:end-1] = @views 0.5 * (q[1:end-1] + q[1+1:end])
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
function (f::VectorField)(x, y, z)
    SA[f.xinterp(x, y, z), f.yinterp(x, y, z), f.zinterp(x, y, z)]
end

function covariant(maingrid::Grid, field) # main cell faces
    mx, my, mz = maingrid.nodes
    dx, dy, dz = _dualgrid(maingrid)
    VectorField(
        Grid(mx, dy, dz),
        Grid(dx, my, dz),
        Grid(dx, dy, mz),
        field
    )
end
function contravariant(maingrid::Grid{T}, field) where T # main cell edges
    mx, my, mz = maingrid.nodes
    dx, dy, dz = _dualgrid(maingrid)
    VectorField(
        Grid(dx, my, mz),
        Grid(mx, dy, mz),
        Grid(mx, my, dz),
        field
    )
end

_dgrid(maingrid) = diff.(maingrid.nodes)
function _dcell(maingrid)
    d = _dgrid(maingrid)
    ret = similar.(maingrid.nodes)
    for (i,q) in pairs(maingrid.nodes)
        ret[i][begin] = d[i][begin]
        ret[i][end] = d[i][end]
        ret[i][begin+1:end-1] .= (
            ((q[begin+2:end] .+ q[begin+1:end-1])./2)
            .- ((q[begin+1:end-1] .+ q[begin:end-2])./2)
        )
    end
    ret
end
function _ratios(maingrid)
    ret = similar.(maingrid.nodes)
    for (i,q) in pairs(maingrid.nodes)
        ret[i] = similar(q)
        ret[i][1] = 0.5
        ret[i][end] = 0.5
        plus = @views (q[3:end] .+ q[2:end-1])./2
        minus = @views (q[2:end-1] .+ q[1:end-2])./2
        ret[i][2:end-1] = @views (q[2:end-1] - minus)/(plus - minus)
    end
    ret
end

const α = e^2*μ_0/m_p
Ep(Bp, cBp, np, up) = -(up - cBp/(α*np)) × Bp

function calculate_Ep_ongrid(prefix, step=-1)
    error("Maybe just delete this function, I never finished writing it")
    # Grid stuff
    para = ParameterSet(prefix)
    nodes = Tuple(convert(Array{Float64}, gp) for gp in para.grid_points)
    maingrid = Grid(nodes)

    # B stuff (B is covariant)
    _, B_data = hr(prefix, "bt").get_timestep(step)
    mion = para.ion_amu
    B = ustrip.(u"T", B_data.*u"s^-1".*(mion*m_p/e))

    # ∇ × B stuff (cB is contravariant)
    cB = similar(B, Union{Missing, eltype(B)})
    cB .= missing
    dx,dy,dz = _spacings(maingrid)
    for i in 2:size(cB,1), j in 2:size(cB,2), k in 2:size(cB,3)
        cB[i,j,k,1] = (B[i,j,k,3] - B[i,j-1,k,3])/dy[j] - (B[i,j,k,2] - B[i,j,k-1,2])/dz[k]
        cB[i,j,k,2] = (B[i,j,k,1] - B[i,j,k-1,1])/dz[k] - (B[i,j,k,3] - B[i-1,j,k,3])/dx[i]
        cB[i,j,k,3] = (B[i,j,k,2] - B[i-1,j,k,2])/dx[i] - (B[i,j,k,1] - B[i,j-1,k,1])/dy[j]
    end

    # density stuff (np_data is on main cell centers, so np_faces is on the main cell faces)
    _, np_data = hr(prefix, "np").get_timestep(step)
    np_faces = similar(B, Union{Missing, eltype(np_data)})
    np_faces .= missing
    np_faces[:,:,:,1] = 0.5*(np_data[1:end-1,:,:] + np_data[2:end,:,:])
    np_faces[:,:,:,2] = 0.5*(np_data[:,1:end-1,:] + np_data[:,2:end,:])
    np_faces[:,:,:,3] = 0.5*(np_data[:,:,1:end-1] + np_data[:,:,2:end])

    # E stuff
    aa = face_to_center(aj)
    a = aa - up
    btc = edge_to_center(B)
    c = cross(a, btc)
    E = c

end

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
function loadfields2(prefix, step=-1)
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

function loadfields(prefix, step=-1)
    E_data = loadvector(prefix, "E", step)u"km/s^2" * (m_p/e)
    B_data = loadvector(prefix, "bt", step)u"s^-1" * (m_p/e)
    E_data = [uconvert.(u"V/m", dat) for dat in E_data]
    B_data = [uconvert.(u"T", dat) for dat in B_data]
    para = ParameterSet(prefix)
    nodes = Tuple(convert(Array{Float64}, gp)*u"km" for gp in para.grid_points)
    E = interpolate(nodes, E_data, Gridded(Linear()))
    B = interpolate(nodes, B_data, Gridded(Linear()))

    # Fixed extrapolations are broken in Interpolations.jl
    #Bsw = SA[0.0, Float64(para.b0_init), 0.0]u"T"
    #Esw = -SA[-Float64(para.vsw), 0.0, 0.0]u"km/s" × Bsw

    E = extrapolate(E, Flat())
    B = extrapolate(B, Flat())
    return E,B
end



end # module
