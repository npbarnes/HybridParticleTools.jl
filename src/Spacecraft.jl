module Spacecraft

export Trajectory

using PyCall
using StaticArrays
using StructArrays
using Unitful

using ..Utility

const st = PyNULL()

function __init__()
    copy!(st, pyimport("spice_tools"))
end

struct SpacecraftState{T,U,V}
    pos::SVector{3,T}
    cmat::SMatrix{3,3,U,9}
    time::V
end
const Trajectory = StructArray{SpacecraftState{T,U,V}} where {T,U,V}
function Trajectory(p::AbstractArray, c::AbstractArray, t::AbstractArray)
    Trajectory{
        eltype(eltype(p)),
        eltype(eltype(c)),
        eltype(t)
    }((p,c,t))
end

function Trajectory(start::Number, stop::Number, step::Number)
    o = pycall(st.trajectory, PyObject, start, stop, step)
    pypos = get(o, PyArray, 0)
    pycmat = get(o, PyArray, 1)

    pos = listofvectors(pypos)*u"km"
    cmat = listofmatrices(pycmat)
    times = get(o, 2)*u"s"
    return Trajectory(pos, cmat, times)
end
end # module
