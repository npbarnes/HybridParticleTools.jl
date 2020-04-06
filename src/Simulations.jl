module Simulations

export ParameterSet, Simulation, touching, Distribution, distributions

using PyCall
using NearestNeighbors: KDTree, inrange
using Distances: Chebyshev
using StaticArrays
using LinearAlgebra: norm
using PhysicalConstants.CODATA2018: m_p, e
using Unitful
using StructArrays

using ..Utility
using ..ParameterSets

function usual_mrat_breakdown(mrat)
    ms = Vector{typeof(1m_p)}(undef, size(mrat,1))
    qs = Vector{typeof(1e)}(undef, size(mrat,1))
    for (i, mr) in enumerate(mrat)
        if mr == 0.5
            ms[i] = 4m_p
            qs[i] = 2e
        else
            ms[i] = (1/mr)*m_p
            qs[i] = 1e
        end
    end
    (ms,qs)
end

struct MacroParticle{T,U,V,W,X}
    x::SVector{3,T}
    v::SVector{3,U}
    m::V
    q::W
    N::X
end
const ParticleData = StructArray{MacroParticle{T,U,V,W,X}} where {T,U,V,W,X}
function ParticleData(x::AbstractArray,v::AbstractArray,m::AbstractArray,q::AbstractArray,N::AbstractArray)
    ParticleData{
        eltype(eltype(x)),
        eltype(eltype(v)),
        eltype(m),
        eltype(q),
        eltype(N)
    }((x,v,m,q,N))
end

mutable struct Simulation{T,U,V,W,X}
    para::ParameterSet
    particles::ParticleData{T,U,V,W,X}
    tree::KDTree
    function Simulation(para, xs, vs, ms, qs, Ns)
        # Leave the tree field undefined
        # the tree field will be defined when it is accessed
        # (see getproperty for this type)
        new{eltype(eltype(xs)), eltype(eltype(vs)), eltype(ms), eltype(qs), eltype(Ns)}(para,ParticleData(xs,vs,ms,qs,Ns))
    end
end

function Simulation(o::PyObject, separate_mrat=usual_mrat_breakdown)
    ax = PyArray(o."x")
    av = PyArray(o."v")

    xs = listofvectors(ax)*u"km"
    vs = listofvectors(av)*u"km/s"

    (ms,qs) = separate_mrat(PyArray(o."mrat"))
    Ns = 1 ./ o.beta
    para = ParameterSet(o."para")

    Simulation(para, xs, vs, ms, qs, Ns)
end
function Simulation(path; n=0, buildtree=false)
    hpr = pyimport("HybridParticleReader")
    ls = hpr.LastStep(path, n=n)
    sim = Simulation(ls)
    if buildtree
        # Accessing the tree the first time causes it to be constructed
        sim.tree
    end
    return sim
end
function Base.getproperty(s::Simulation, name::Symbol)
    if name in fieldnames(MacroParticle)
        return getproperty(getfield(s,:particles), name)
    elseif name in fieldnames(ParameterSet)
        return getproperty(getfield(s,:para), name)
    elseif name === :tree
        # Since generating the tree is expensive, we defer until it is requested.
        if isdefined(s, :tree)
            return getfield(s, :tree)
        else
            return setfield!(s, :tree, KDTree([ustrip.(u"km",x) for x in s.x], Chebyshev()))
        end
    else
        error("type $(typeof(s)) has no property $name")
    end
end
Base.propertynames(::Simulation) = (fieldnames(ParticleData)..., fieldnames(ParameterSet)..., :tree)
Broadcast.broadcastable(s::Simulation) = Ref(s)

touching(s::Simulation, x::AbstractVector{<:Quantity}, dx::Quantity) = inrange(s.tree, ustrip.(u"km",x), ustrip(u"km",dx))
touching(s::Simulation, x::AbstractVector{<:Number}, dx=s.dx) = inrange(s.tree, x, dx)

end # module
