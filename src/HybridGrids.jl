module HybridGrids

export loadfields, loadvector

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

function loadvector(prefix, name, step=-1)
    h = hr(prefix, name)
    _, raw = h.get_timestep(step)
    ret = similar(raw, SVector{3,Float64}, size(raw)[1:end-1])
    for i in axes(raw,1), j in axes(raw,2), k in axes(raw,3)
        ret[i,j,k] = SVector{3,Float64}(raw[i,j,k,:])
    end
    return ret
end

# There was a problem writing bt data for 2020-Jan-23/pluto-2.
# Timestep 27 is the last good step.
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
