using HybridTools
using Dates
using PyCall
using PyPlot
using StaticArrays
using Unitful
using LinearAlgebra
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

folder = "/media/nathan/DATAPART11/2020-Thu-Jan-23/pluto-2"
const np = pyimport("numpy")
const plt = pyimport("matplotlib.pyplot")
const st = pyimport("spice_tools")
const s = Simulation(folder, n=12)

const pepssifov = fov_polygon("NH_PEPSSI_S0", st.nh_in_wake)
const swapfov = fov_polygon("NH_SWAP", st.nh_in_wake)

const et_start = utc2et"11:30"
const et_end   =  utc2et"14:01"
ets = collect(et_start:(10*60):et_end)

ipuidist(d::Distribution) = filter(x->x.t==He_ipui, d)
ipuidist(pos::AbstractVector) = ipuidist(Distribution(s,pos))
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
const E,B = loadfields(prefix, 27)

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
    y_extent::Tuple{U,U}
    z_extent::Tuple{U,U}
    N::Int64
    dt::T
    function UnBufferedGenerator(v_sampler, x, ys, zs, N, dt)
        new{typeof(v_sampler), typeof(x), typeof(dt)}(v_sampler,x,ys,zs,N,dt)
    end
end
function (g::UnBufferedGenerator)(lst)
    ymin,ymax = g.y_extent
    zmin,zmax = g.z_extent
    for i in 1:g.N
        v = g.v_sampler()
        x = SA[g.x + rand()*v[1]*g.dt, (ymax-ymin)*rand()+ymin, (zmax-zmin)*rand()+zmin]
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
function superthermal()
    mag = rand(Pareto(4, 400))u"km/s"
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

const f = Boris.Fields(E,B)
const dom = Boris.Domain(extrema.(parent(E).knots))
finaltime = 50000*dt
step = 10dt
N = ceil(Int, finaltime/step)
shell_pg! = UnBufferedGenerator(shell, dom.max_x, (dom.min_y,dom.max_y), (dom.min_z,dom.max_z), 50000, step)
superthermal_pg! = UnBufferedGenerator(superthermal, dom.max_x, (dom.min_y,dom.max_y), (dom.min_z,dom.max_z), 13000, step)

ps = trace_particles(superthermal_pg!, f, dom, step, N)
kd = KDTree([ustrip.(u"km",p.x) for p in ps], Chebyshev());

function _traced_distribution(ps, kd, x, dx, N)
    indexs = inrange(kd, ustrip.(u"km",x), ustrip(u"km",dx))
    l = length(indexs)
    n = Vector{typeof(1.0u"km^-3")}(undef, l)
    for (i,index) in enumerate(indexs)
        n[i] = N*HybridTools.Distributions.weight(x, ps[index].x, dx)
    end
    Distribution{eltype(ps[1].v), typeof(4m_p), typeof(1e), typeof(1.0u"km^-3")}((
        getproperty.(@view(ps[indexs]),:v),
        fill(4m_p, l),
        fill(1e, l),
        n,
        fill(He_ipui, l)
    ))
end

flyby_ets = utc2et"11:00":30:utc2et"13:00";
ds = _traced_distribution.(Ref(ps), Ref(kd), location.(flyby_ets), s.dx*u"km", macroparticle_N(superthermal_pg!, 400u"km/s"*4.7e10u"km^-3"))
fig, ax = especfigure();
PlottingTools.plot_espec_scatter(fig, ax, flyby_ets, filter.(pepssifov, ds));

nothing
