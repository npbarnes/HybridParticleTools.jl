using HybridTools
using Dates
using PyCall
using PyPlot
using Plots
plotlyjs()
using StaticArrays
using Unitful
using LinearAlgebra
using IterTools: firstrest
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
const st = pyimport("spice_tools")
const s = Simulation(pwd(), n=12)

const pepssifov = fov_polygon("NH_PEPSSI_S0", st.nh_in_wake)
const swapfov = fov_polygon("NH_SWAP", st.nh_in_wake)

const et_start = utc2et"11:30"
const et_end   =  utc2et"14:01"
const ets = collect(et_start:(10*60):et_end)

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
s0(ets) = filter.(fov_polygon.("NH_PEPSSI_S0", ets), ipuidist.(ets))
expandby(a, n) = reduce(vcat, fill.(a,n))
const ds = s0(ets)
const vvs = getproperty.(ds,:v)
const vs = reduce(vcat, vvs)
const xs = expandby(location.(ets), length.(vvs))
const b = boris.(xs, vs, dt, Ref(E), Ref(B), 10000, 4m_p);
moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]
smooth_bulk_v(ds) = moving_average(filter(!isnan, mean.(getproperty.(ds, :v))), 500)

nothing
