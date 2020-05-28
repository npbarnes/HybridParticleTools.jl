module PlottingTools
export plt, circlegrid, map_projection, vec_coords, geodetic, mapcoords, ccrs, plotdist, plotshape

using PyCall
using LinearAlgebra
using Unitful
using ..Utility
using ..SphericalShapes

const plt = PyNULL()
const Axes3D = PyNULL()
const ccrs = PyNULL()
const unit_sphere = PyNULL()
const map_projection = PyNULL()
const vec_coords = PyNULL()
const geodetic = PyNULL()

function __init__()
   copy!(plt, pyimport("matplotlib.pyplot"))
   copy!(Axes3D, pyimport("mpl_toolkits.mplot3d").Axes3D)
   copy!(ccrs, pyimport("cartopy.crs"))
   copy!(unit_sphere, ccrs.Globe(semimajor_axis=1., semiminor_axis=1., ellipse=nothing))
   copy!(map_projection, ccrs.Mollweide(globe=unit_sphere))
   copy!(vec_coords, map_projection.as_geocentric())
   copy!(geodetic, map_projection.as_geodetic())
end

function mapcoords(vs; from_crs=vec_coords, to_crs=map_projection)
   mc = to_crs.transform_points(from_crs,
      x=getindex.(vs,1),
      y=getindex.(vs,2),
      z=getindex.(vs,3)
   )
   @views (mc[:,1], mc[:,2])
end

py"""
import numpy as np
def plot_sphere(ax, r, x_offset=0, **kwargs):
   u=np.linspace(0,2*np.pi,100)
   v=np.linspace(0,np.pi,100)
   x = r * np.outer(np.cos(u), np.sin(v)) - x_offset
   y = r * np.outer(np.sin(u), np.sin(v))
   z = r * np.outer(np.ones(np.size(u)), np.cos(v))
   return ax.plot_surface(x,y,z,**kwargs)
"""
plot_sphere(ax, r, x_offset=0; kwargs...) = py"plot_sphere"(ax, r, xoffset, kwargs...)

function plot_3d_dist(ax, d)
    ax.scatter(getindex.(d.v, 1), getindex.(d.v,2), getindex.(d.v,3))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plot_sphere(ax, 100, 400, alpha=0.6)
end

function plotshape(ax, sp; kwargs...)
   es = collect(Iterators.flatten(edges(sp)))
   ax.plot(mapcoords(es, to_crs=geodetic)...; transform=geodetic, kwargs...)
end

function circlepoints(a,N=100)
   [(cos(a), sin(a)*sin(t), -sin(a)*cos(t)) for t in 2π/N*(0:N)]
end

function circlegrid(ax, N=5)
   for a in π/(N+1)*(1:N)
      ax.plot(mapcoords(circlepoints(a), to_crs=geodetic)..., transform=geodetic, color="gray", linewidth=0.5, alpha=0.75)
   end
end

function mapfigure()
    fig, ax = plt.subplots(subplot_kw=Dict("projection"=>map_projection))
    ax.set_global()
    circlegrid(ax)
    ax.background_patch.set_facecolor("lightgray")
    return fig, ax
end

pepssi_fov_polys(et) = [fov_polygon("NH_PEPSSI_S$(i)", et) for i in 0:5]
function plot_pepssi(ax, et)
    polys = pepssi_fov_polys(et)
    plotshape(ax, polys[1], color="blue")
    for p in @view polys[2:end]
        plotshape(ax, p, color="darkgray")
    end
end

function plot_dist(fig, ax, d; marker=".", s=30.0, kwargs...)
   l = -asunitless(d.v)
   E = ustrip(uconvert.(u"keV", energy.(d)))
   mappable = ax.scatter(mapcoords(l)...; marker=marker, s=s, c=E, kwargs...)
   cb = fig.colorbar(mappable)
   cb.set_label("Energy (keV)")
end
plot_sun(ax) = ax.scatter(mapcoords([[1.,0.,0.]])..., marker="*", edgecolors="k", color="gold", s=200)
plot_pluto(ax,pos::AbstractArray) = ax.scatter(mapcoords([-ustrip(pos)])..., marker=raw"$♥$",  edgecolors="k", color="chocolate", s=100)
plot_pluto(ax,et::Number) = plot_pluto(location(et))

end # module
