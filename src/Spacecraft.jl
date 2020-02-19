module Spacecraft

export Trajectory

using PyCall
using StaticArrays
using Unitful

using ..Utility


struct Trajectory{T,U,V}
    pos::Vector{SVector{3,T}}
    cmat::Vector{SMatrix{3,3,U,9}}
    times::Vector{V}
    function Trajectory(pos::AbstractVector{<:AbstractVector{<:Number}}, cmat::AbstractVector{<:AbstractMatrix{<:Number}}, times::AbstractVector{<:Number})
        N = length(times)
        if length(pos) != N || length(cmat) != N
            error("The lists must have the same lengths")
        end
        new{eltype(eltype(pos)), eltype(eltype(cmat)), eltype(times)}(pos, cmat, times)
    end
end
function Trajectory(start::Number, stop::Number, step::Number)
    spice_tools = pyimport("spice_tools")
    o = pycall(spice_tools.trajectory, PyObject, start, stop, step)
    pypos = get(o, PyArray, 0)
    pycmat = get(o, PyArray, 1)

    pos = listofvectors(pypos)*u"km"
    cmat = listofmatrices(pycmat)
    times = get(o, 2)*u"s"
    return Trajectory(pos, cmat, times)
end
end # module
