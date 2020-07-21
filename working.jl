using HybridTools
using Dates
using PyCall
using PyPlot
#using Plots
#plotlyjs()
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

const np = pyimport("numpy")
const plt = pyimport("matplotlib.pyplot")
const st = pyimport("spice_tools")
const s = Simulation(pwd(), n=12)

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

const prefix = joinpath(pwd(),"data")
const para = ParameterSet(prefix)
const dt = Float64(para.dt)*u"s"
#const xinit = SA[8.0, 0.0, 0.0]u"Rp"
const E,B = loadfields(prefix, 27)
#const S0_ds = filter.(fov_polygon.("NH_PEPSSI_S0", ets), ipuidist.(ets))

function split(a)
    return (
        getindex.(a,1),
        getindex.(a,2),
        getindex.(a,3)
    )
end

#=
function plotlines(splitpairs)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    n = 100
    for ((x,y,z),(vx,vy,vz)) in splitpairs
        ax.plot(ustrip(x[1:n:end]), ustrip(y[1:n:end]), ustrip(z[1:n:end]))
    end
    ax.scatter([0],[0],[0], s=100)
    ax.set_xlim((-50,20))
    ax.set_ylim((-35,35))
    ax.set_zlim((-35,35))
end
function plotenergies(splitpairs)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    n = 100
    for ((x,y,z),(vx,vy,vz)) in splitpairs
        xx = ustrip(x[1:n:end])
        yy = ustrip(y[1:n:end])
        zz = ustrip(z[1:n:end])
        vxx = ustrip(vx[1:n:end])
        vyy = ustrip(vy[1:n:end])
        vzz = ustrip(vz[1:n:end])
        ax.scatter(xx, yy, zz, c=ustrip.(vxx.^2 .+ vyy.^2 .+ vzz.^2), s=10)
    end
    ax.scatter([0],[0],[0], s=100)
    ax.set_xlim((-50,20))
    ax.set_ylim((-35,35))
    ax.set_zlim((-35,35))
end

function todaysplot(b, c=["b"], N=100)
    #fig, ax = plt.subplots(subplot_kw=Dict("projection"=>"3d"))
    i = 0
    xx = b[1][1]
    curx = xx[1]
    x,y,z = split(ustrip.(xx[1:N:end]))
    p = plot(x,y,z, color=c[i%length(c)+1])
    for (xx,vv) in b
        if xx[1] != curx
            curx = xx[1]
            i += 1
        end
        x,y,z = split(ustrip.(xx[1:N:end]))
        plot!(p, x,y,z, color=c[i%length(c)+1], legend=false)
    end
    scatter!(p, [0],[0],[0], markersize=1, color="beige")
    plot!(p, split(ustrip.(location.(ets)))..., color="grey", linewidth=3)
    p
end
=#
#=
s0(ets) = filter.(fov_polygon.("NH_PEPSSI_S0", ets), ipuidist.(ets))
expandby(a, n) = reduce(vcat, fill.(a,n))
const ds = s0(ets)
const vvs = getproperty.(ds,:v)
const vs = reduce(vcat, vvs)
const xs = expandby(location.(ets), length.(vvs))
const b = boris.(xs, vs, dt, Ref(E), Ref(B), 10000, 4m_p);
moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]
smooth_bulk_v(ds) = moving_average(filter(!isnan, mean.(getproperty.(ds, :v))), 500)
=#

######################################################################################
mutable struct BufferedGenerator{S,F,T,U}
    v_sampler::S
    fields::F
    bufferdomain::Boris.Domain{T}
    N::Int64
    dt::U
    function BufferedGenerator(v_sampler, f, d::Boris.Domain{T}, N, dt) where {T,U,V}
        new{typeof(v_sampler), typeof(f), T, typeof(dt)}(v_sampler,f,d,N,dt)
    end
end
const generator_lock = SpinLock()
const RNGs = [MersenneTwister() for i in 1:nthreads()]
function (bg::BufferedGenerator)(lst)
    length_before = length(lst)
    @threads for i in 1:bg.N
        rng = RNGs[threadid()]
        p = Boris.Particle(rand(rng, bg.bufferdomain), bg.v_sampler(rng), e/(4m_p))
        p.x,p.v = Boris.step(p, bg.fields, bg.dt)
        if p.x[1] < bg.bufferdomain.min_x
            lock(generator_lock) do
                push!(lst, p)
            end
        end
    end
    numberpushed = length(lst) - length_before
    if numberpushed/bg.N > 0.1
        bg.dt = bg.dt/2
    elseif numberpushed/bg.N < 0.01
        bg.dt = 2bg.dt
    end
end
maxwell() = sqrt(1.5u"eV"/(4m_p)).*@SVector(randn(3)) + SA[-400.0,0.0,0.0]u"km/s"
function sphere_sample(rng)
    u = 2rand(rng) - 1
    v = 2rand(rng) - 1
    s = u^2 + v^2
    while s >= 1
       u = 2rand(rng) - 1
       v = 2rand(rng) - 1
       s = u^2 + v^2
    end
    z = sqrt(1-s)
    SA[2u*z, 2v*z, 1-2s]
end
function superthermal(rng=Random.GLOBAL_RNG)
    mag = rand(rng,Pareto(4, 400))u"km/s"
    mag*sphere_sample(rng) + SA[-400.0, 0.0, 0.0]u"km/s"
end
function shell(rng=Random.GLOBAL_RNG)
    mag = 400.0u"km/s"
    mag*sphere_sample(rng) + SA[-400.0, 0.0, 0.0]u"km/s"
end

const f = Boris.Fields(E,B)
const dom = Boris.Domain(extrema.(parent(E).knots))
const buf_dom = Boris.Domain(dom.max_x, dom.min_y, dom.min_z, dom.max_x+2s.dx*u"km", dom.max_y, dom.max_z)
const superthermal_pg! = BufferedGenerator(superthermal, f, buf_dom, 9000000, dt)
finaltime = 10000*dt
step = 10dt
N = finaltime/step
const ps = trace_particles(superthermal_pg!, f, dom, step, N, saveat=N);
const kd = KDTree([ustrip.(u"km",p.x) for p in ps], Chebyshev())
const ts = inrange.(Ref(kd), ustrip.(u"km",location.(ets)), s.dx)
Distribution(ps::Vector{<:Boris.Particle}) = Distribution(getproperty.(ps,:v),
                                                        fill(4m_p, length(ps)),
                                                        fill(e, length(ps)),
                                                        fill(1u"km^-3", length(ps)),
                                                        fill(He_ipui, length(ps)))

function getparticlehistory(res,i)
    getindex.(res[length.(res) .>= i], i)
end

function plot_particlehistory(ax,h)
    xs = asunitless(u"Rp", getproperty.(h,:x))
    ax.plot(split(xs)...)
end



nothing
