using HybridTools
using Dates
using PyCall
using PyPlot
using StaticArrays
using Unitful
using LinearAlgebra
using Statistics
using NearestNeighbors
using Random
using Base.Threads
using IterTools: firstrest
using Distributions: Pareto
import PhysicalConstants.CODATA2018: m_p, e

using HybridTools.SphericalShapes
using HybridTools.Utility
using HybridTools.Sensors
using HybridTools.PlottingTools
using HybridTools.HybridGrids
using HybridTools.Boris
using HybridTools.ParameterSets
using HybridTools.PlutoUnits
using HybridTools.Spacecraft

#folder = "/media/nathan/DATAPART11/2020-Thu-Jan-23/pluto-2"
folder = "/media/nathan/DATAPART11/2020-Mon-Aug-10/pluto-1"
const np = pyimport("numpy")
const plt = pyimport("matplotlib.pyplot")
const st = pyimport("spice_tools")
const s = Simulation(folder, n=6)

const pepssifov = fov_polygon("NH_PEPSSI_S0", st.nh_in_wake)
const swapfov = fov_polygon("NH_SWAP", st.nh_in_wake)

const et_start = utc2et"11:30"
const et_end   =  utc2et"14:01"
ets = collect(et_start:(10*60):et_end)

Hipuidist(d::Distribution) = filter(x->x.t==H_ipui, d)
Hipuidist(pos::AbstractVector) = Hipuidist(Distribution(s,pos))
Hipuidist(et::Number) = Hipuidist(location(et))

ipuidist(d::Distribution) = filter(x->x.t==He_ipui, d)
ipuidist(pos::AbstractVector) = ipuidist(Distribution(s,pos, 2s.dx*u"km"))
ipuidist(et::Number) = ipuidist(location(et))

swdist(d::Distribution) = filter(x->x.t==H_sw || x.t==He_sw, d)
swdist(pos::AbstractVector) = swdist(Distribution(s,pos))
swdist(et::Number) = swdist(location(et))

notipuidist(d::Distribution) = filter(x->x.t!=H_ipui && x.t!=He_ipui && x.t!=dummy, d)
notipuidist(pos::AbstractVector) = notipuidist(Distribution(s,pos))
notipuidist(et::Number) = notipuidist(location(et))

const prefix = joinpath(folder,"data")
const para = ParameterSet(prefix)
const dt = Float64(para.dt)*u"s"
#const E,B = loadfields2(prefix, 27)
const E,B = loadfields2(prefix)

function split(a)
    return (
        getindex.(a,1),
        getindex.(a,2),
        getindex.(a,3)
    )
end

struct UnBufferedGenerator{S,U,T}
    v_sampler::S
    x::U
    dx::U
    y_extent::Tuple{U,U}
    z_extent::Tuple{U,U}
    N::Int64
    dt::T
    function UnBufferedGenerator(v_sampler, x, dx, ys, zs, N, dt)
        new{typeof(v_sampler), typeof(x), typeof(dt)}(v_sampler,x,dx,ys,zs,N,dt)
    end
end
function (g::UnBufferedGenerator)(lst)
    ymin,ymax = g.y_extent
    zmin,zmax = g.z_extent
    v = g.v_sampler()
    x = @SVector zeros(typeof(g.x), 3)
    for i in 1:g.N
        while true # rejection sampling step
            v = g.v_sampler()
            if v[1] > zero(v[1]) # early rejection if the paricle is moving upstream
                continue
            end
            x0 = g.x - rand()*g.dx
            x1 = x0 + v[1]*g.dt
            if x1 <= g.x-g.dx
                x = SA[x1, (ymax-ymin)*rand()+ymin, (zmax-zmin)*rand()+zmin]
                break
            end
        end
        push!(lst, Boris.Particle(x, v, e/(4m_p)))
    end
end
function macroparticle_N(ubg::UnBufferedGenerator, trueflux)
    A = (ubg.y_extent[2] - ubg.y_extent[1])*(ubg.z_extent[2] - ubg.z_extent[1])
    macroflux = ubg.N/(A*ubg.dt)
    uconvert(Unitful.NoUnits, trueflux/macroflux)
end

maxwell() = sqrt(1.5u"eV"/(4m_p)).*@SVector(randn(3)) + SA[-400.0,0.0,0.0]u"km/s"
function sphere_sample()
    u = 2rand() - 1
    v = 2rand() - 1
    s = u^2 + v^2
    while s >= 1
       u = 2rand() - 1
       v = 2rand() - 1
       s = u^2 + v^2
    end
    z = sqrt(1-s)
    SA[2u*z, 2v*z, 1-2s]
end
function superthermal(p=5)
    mag = rand(Pareto(p-1, 400))u"km/s"
    mag*sphere_sample() + SA[-400.0, 0.0, 0.0]u"km/s"
end
function superthermal2()
    mag = (400rand()+400)u"km/s"
    mag*sphere_sample() + SA[-400.0, 0.0, 0.0]u"km/s"
end
function shell()
    mag = 400.0u"km/s"
    mag*sphere_sample() + SA[-400.0, 0.0, 0.0]u"km/s"
end

function _inv_F(y, va, vb)
    f0 = 1/(5/4*vb - va)
    C = f0*(vb-va)
    if y < C
        return y/f0 + va
    else
        return (4/(f0*vb^5)*(C - y) + vb^-4)^(-1/4)
    end
end
"""va is lower cutoff velocity in sw frame, vb is injection velocity"""
function flat_w_superthermal(va=0u"km/s", vb=219u"km/s")
    y = rand()
    mag = _inv_F(y, va, vb)
    mag*sphere_sample() + SA[-400.0, 0.0, 0.0]u"km/s"
end

const f = Boris.Fields(E,B)
#const dom = Boris.Domain(Tuple(t.*u"km" for t in extrema.(E.xgrid.nodes)))


# approx bins for PEPSSI Helium (I assume that means signly ionized):
const pepssi_bin_names = ["L13", "L11", "L09", "L07", "L05", "L03", "L01"]
const pepssi_bin_edges = [2.2, 4.79, 9.485, 18.05, 33.4, 59.8, 106, 183]u"keV"
const pepssi_S0_area = 0.092677

pepssi_view(R, time::String) = pepssi_view(R, utc2et(time))
function pepssi_view(R, time::Float64)
    d = ipuidist(time)
    d.v .= Ref(R) .* d.v
    fig, ax = plot_dist(d)
    #plotshape(ax, fov_polygon("NH_PEPSSI_S0", utc2et"12:18"), color="blue")

    fig, ax
end

nothing
