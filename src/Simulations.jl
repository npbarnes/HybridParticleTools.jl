module Simulations

export Parameters, Simulation, touching, Distribution, distributions

using PyCall
using NearestNeighbors: KDTree, inrange
using Distances: Chebyshev
using StaticArrays
using LinearAlgebra: norm
using PhysicalConstants.CODATA2018: m_p, e
using Unitful

using ..Utility

## Parameters
struct Parameters
    nx::Int32
    ny::Int32
    nz::Int32
    dx::Float32
    dy::Float32
    delz::Float32
    nt::Int32
    dtsub_init::Float32
    ntsub::Int32
    dt::Float32
    nout::Int32
    out_dir::String
    vtop::Float32
    vbottom::Float32
    Ni_max::Int32
    mproton::Float32
    m_pu::Float32
    m_heavy::Float32
    np_top::Float32
    np_bottom::Float32
    b0_top::Float32
    b0_bottom::Float32
    vth_top::Float32
    vth_bottom::Float32
    alpha::Float64
    beta::Float32
    RPluto::Float32
    RIo::Float32
    b0_init::Float32
    ion_amu::Int32
    mpu::Float32
    nf_init::Float32
    dt_frac::Float32
    vsw::Float32
    vth::Float32
    Ni_tot_frac::Float32
    dx_frac::Float32
    nu_init_frac::Float32
    mrestart::Int32
    pluto_offset::Int64
    ri0::Int32
    part_nout::Int32
    qx::Array{Float32,1}
    qy::Array{Float32,1}
    qz::Array{Float32,1}
    dz_grid::Array{Float32,1}
    dz_cell::Array{Float32,1}
    num_proc::Int64
    zrange::Int64
    saved_steps::Float64
    simulation_height::Float64
    qzrange::Array{Float64,1}
    grid_points::Tuple{Array{Float32,1},Array{Float32,1},Array{Float64,1}}
end
Parameters(o::PyObject) = Parameters(convert(Dict{Symbol,Any}, o))
function Parameters(d::Dict{Symbol,Any})
    Parameters(
        d[:nx],
        d[:ny],
        d[:nz],
        d[:dx],
        d[:dy],
        d[:delz],
        d[:nt],
        d[:dtsub_init],
        d[:ntsub],
        d[:dt],
        d[:nout],
        d[:out_dir],
        d[:vtop],
        d[:vbottom],
        d[:Ni_max],
        d[:mproton],
        d[:m_pu],
        d[:m_heavy],
        d[:np_top],
        d[:np_bottom],
        d[:b0_top],
        d[:b0_bottom],
        d[:vth_top],
        d[:vth_bottom],
        d[:alpha],
        d[:beta],
        d[:RPluto],
        d[:RIo],
        d[:b0_init],
        d[:ion_amu],
        d[:mpu],
        d[:nf_init],
        d[:dt_frac],
        d[:vsw],
        d[:vth],
        d[:Ni_tot_frac],
        d[:dx_frac],
        d[:nu_init_frac],
        d[:mrestart],
        d[:pluto_offset],
        d[:ri0],
        d[:part_nout],
        d[:qx],
        d[:qy],
        d[:qz],
        d[:dz_grid],
        d[:dz_cell],
        d[:num_proc],
        d[:zrange],
        d[:saved_steps],
        d[:simulation_height],
        d[:qzrange],
        d[:grid_points]
    )
end
Base.show(io::IO, hp::Parameters) = print(io, "Parameters(...)")
function Base.show(io::IO, ::MIME"text/plain", hp::Parameters)
    print(io, "HybridParams:")
    for f in fieldnames(typeof(hp))
        print(io, "\n\t$f: $(getfield(hp,f))")
    end
end

function Parameters(path)
    pymod = pyimport("HybridParams")
    pyhp = pymod.HybridParams(path)
    return Parameters(pyhp."para")
end

## Hybrid Simulation Particle Data
struct BasicSimulation{T,U,V,W,X}
    para::Parameters
    xs::Vector{SVector{3,T}}
    vs::Vector{SVector{3,U}}
    ms::Vector{V}
    qs::Vector{W}
    Ns::Vector{X}
    function BasicSimulation(para, xs, vs, ms, qs, Ns)
        len = length(xs)
        if len != length(vs) || len != length(ms) || len != length(qs) || len != length(Ns)
            error("each vector should have the same length")
        end
        new{eltype(eltype(xs)), eltype(eltype(vs)), eltype(ms), eltype(qs), eltype(Ns)}(para,xs,vs,ms,qs,Ns)
    end
end
mutable struct Simulation{T,U,V,W,X}
    sim::BasicSimulation{T,U,V,W,X}
    tree::KDTree
    function Simulation(para, xs, vs, ms, qs, Ns)
        # Leave the tree field undefined
        # the tree field will be defined when it is accessed
        # (see getproperty for this type)
        new{eltype(eltype(xs)), eltype(eltype(vs)), eltype(ms), eltype(qs), eltype(Ns)}(BasicSimulation(para,xs,vs,ms,qs,Ns))
    end
end

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

function Simulation(o::PyObject, separate_mrat=usual_mrat_breakdown)
    ax = PyArray(o."x")
    av = PyArray(o."v")

    xs = listofvectors(ax)*u"km"
    vs = listofvectors(av)*u"km/s"

    (ms,qs) = separate_mrat(PyArray(o."mrat"))
    Ns = 1 ./ o.beta
    para = Parameters(o."para")

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
    if name in fieldnames(BasicSimulation)
        return getproperty(getfield(s,:sim), name)
    elseif name in fieldnames(Parameters)
        return getproperty(getfield(s,:sim).para, name)
    elseif name === :tree
        # Since generating the tree is expensive, we defer until it is requested.
        if isdefined(s, :tree)
            return getfield(s, :tree)
        else
            return setfield!(s, :tree, KDTree([ustrip.(u"km",x) for x in s.xs], Chebyshev()))
        end
    else
        error("type $(typeof(s)) has no property $name")
    end
end
Base.propertynames(::Simulation) = (fieldnames(BasicSimulation)..., fieldnames(Parameters)..., :tree)
Broadcast.broadcastable(s::Simulation) = Ref(s)

"If x doesn't have units, it is assumed to be km"
function touching(s::Simulation, x::AbstractVector{Quantity{T,Unitful.𝐋,Unitful.FreeUnits{(u"km",),Unitful.𝐋,nothing}}}, dx=s.dx) where T
    inrange(s.tree, ustrip(x), dx)
end
function touching(s::Simulation, x::AbstractVector{<:Quantity}, dx=s.dx)
    inrange(s.tree, ustrip.(u"km", x), dx)
end
touching(s::Simulation, x::AbstractVector{<:Real}, dx=s.dx) = inrange(s.tree, x, dx)

struct Distribution{U<:Number,V<:Number,W<:Number,X<:Number}
    vs::Vector{SVector{3,U}}
    ms::Vector{V}
    qs::Vector{W}
    ns::Vector{X}
    function Distribution(vs, ms, qs, ns)
        N = length(vs)
        if length(ms) != N || length(qs) != N  || length(ns) != N
            error("lists must have the same lengths")
        end
        new{eltype(eltype(vs)), eltype(ms), eltype(qs), eltype(ns)}(vs, ms, qs, ns)
    end
end

function weight(x,xp,dx)
    (abs(1 - abs(x[1]-xp[1])/dx) * abs(1 - abs(x[2]-xp[2])/dx) * abs(1 - abs(x[3]-xp[3])/dx)) / dx^3
end

function Distribution(s::Simulation, x::AbstractVector{<:Number}, dx=s.dx*u"km")
    t = touching(s, x)
    ns = similar(s.Ns[t])*u"km^-3"
    for (i,(xp,N)) in enumerate(zip(s.xs[t],s.Ns[t]))
        ns[i] = N*weight(x,xp,dx)
    end
    Distribution(s.vs[t], s.ms[t], s.qs[t], ns)
end

end # module
