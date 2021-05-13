module PlottingTools
export  plt, circlegrid, map_projection, vec_coords, geodetic, mapcoords, ccrs, plotdist, plotshape,
        plot_colored_line, mapfigure, plot_dist, plot_dist_hist, plot_pepssi, plot_sun, plot_pluto, especfigure,
        plot_espec_scatter, plot_pepssi_view, datefigure, plotshape_patch, diff_inten_ongrid, platecarree, polygoncolor, LogNorm

using PyCall
using PyPlot
using LinearAlgebra
using Unitful
using IterTools
using StaticArrays
using ..Utility
using ..SphericalShapes
using ..Sensors
using ..Spacecraft
using ..Distributions
using ..Simulations

const Axes3D = PyNULL()
const ccrs = PyNULL()
const unit_sphere = PyNULL()
const map_projection = PyNULL()
const vec_coords = PyNULL()
const geodetic = PyNULL()
const platecarree = PyNULL()
const Rectangle = PyNULL()
const LogNorm = PyNULL()
const mcoll = PyNULL()

function __init__()
    copy!(Axes3D, pyimport("mpl_toolkits.mplot3d").Axes3D)
    copy!(ccrs, pyimport("cartopy.crs"))
    copy!(unit_sphere, ccrs.Globe(semimajor_axis=1., semiminor_axis=1., ellipse=nothing))
    copy!(map_projection, ccrs.Mollweide(globe=unit_sphere))
    copy!(vec_coords, map_projection.as_geocentric())
    copy!(geodetic, map_projection.as_geodetic())
    copy!(platecarree, ccrs.PlateCarree(globe=unit_sphere))
    copy!(Rectangle, pyimport("matplotlib.patches").Rectangle)
    copy!(LogNorm, pyimport("matplotlib.colors").LogNorm)
    copy!(mcoll, pyimport("matplotlib.collections"))

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

end # __init__

plot_sphere(ax, r, x_offset=0; kwargs...) = py"plot_sphere"(ax, r, xoffset, kwargs...)

function mapcoords(vs; from_crs=vec_coords, to_crs=map_projection)
    mc = to_crs.transform_points(from_crs,
        x=getindex.(vs,1),
        y=getindex.(vs,2),
        z=getindex.(vs,3)
    )
    # Use the negative of the longitudes because we're looking from the inside of the sphere
    @views (-mc[:,1], mc[:,2])
end

function plot_3d_dist(ax, d)
    ax.scatter(getindex.(d.v, 1), getindex.(d.v,2), getindex.(d.v,3))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plot_sphere(ax, 100, 400, alpha=0.6)
end

function plotshape(ax, sp::SphericalPolygon; kwargs...)
    es = collect(Iterators.flatten(edges(sp)))
    ax.plot(mapcoords(es, to_crs=geodetic)...; transform=geodetic, kwargs...)
end
function plotshape(ax, sc::SCircle; kwargs...)
    points = circlepoints(sc.center, sc.angle)
    ax.plot(mapcoords(points, to_crs=geodetic)...; transform=geodetic, kwargs...)
end

function plotshape(ax, r::LLAlignedRectangle; kwargs...)
    ax.plot(rad2deg.([-r.lon[1], -r.lon[2], -r.lon[2], -r.lon[1]]), [r.lat[1], r.lat[1], r.lat[2], r.lat[2]]; transform=platecarree, kwargs...)
end
function plotshape_patch(ax, r::LLAlignedRectangle; kwargs...)
    rect = Rectangle((-r.lon[1], r.lat[1]), r.lon[1]-r.lon[2], r.lat[2]-r.lat[1]; transform=platecarree, kwargs...)
    ax.add_patch(rect)
end

function circlepoints(c::AbstractVector,a::Number,N::Int=100)
    # QQ is an othogonal matrix that rotates (1,0,0) to c
    A = [c I]
    Q,R = qr(A)
    QQ = [c Q[:,2:end]]
    # Get the circle points around (1,0,0) and then rotate them using QQ
    return [QQ*v for v in circlepoints(a,N)]
end
function circlepoints(a::Number,N::Int=100)
    [SA[cos(a), sin(a)*sin(t), -sin(a)*cos(t)] for t in 2π/N*(0:N)]
end

function circlegrid(ax, N=5)
    for a in π/(N+1)*(1:N)
        ax.plot(mapcoords(circlepoints(a), to_crs=geodetic)..., transform=geodetic, color="gray", linewidth=0.5, alpha=0.75)
    end
end

function mapfigure(showgrid=true)
    fig, ax = plt.subplots(subplot_kw=Dict("projection"=>map_projection))
    ax.set_global()
    if showgrid
        circlegrid(ax)
    end
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

plot_dist(d;kwargs...) = plot_dist(mapfigure()..., d; kwargs...)
function plot_dist(fig, ax, d; marker=".", s=30.0, kwargs...)
    l = -asunitless(d.v)
    E = ustrip(uconvert.(u"keV", energy.(d)))
    mappable = ax.scatter(mapcoords(l)...; marker=marker, s=s, c=E, kwargs...)
    cb = fig.colorbar(mappable)
    cb.set_label("Energy (keV)")
    fig, ax
end
function plot_dist_hist(fig, ax, d;kwargs...)
    np = pyimport("numpy")
    l = -asunitless(d.v)
    H, xs, ys = np.histogram2d(mapcoords(l)..., bins=100, weights=ustrip.(d.n))
    ax.pcolormesh(xs, ys, H, transform=map_projection)
    #h, xs, ys, mappable = ax.hist2d(mapcoords(l)..., bins=100, weights=ustrip.(d.n); kwargs...)
    ax.set_global()
end

plot_sun(ax) = ax.scatter(mapcoords([[1.,0.,0.]])..., marker="*", edgecolors="k", color="gold", s=200)
plot_pluto(ax,pos::AbstractArray) = ax.scatter(mapcoords([-ustrip(pos)])..., marker=raw"$♥$",  edgecolors="k", color="chocolate", s=100)
plot_pepssi_view(et::Float64, d::Distribution, args...) = plot_pepssi_view(mapfigure()..., et, d, args...)
#=
function plot_pepssi_view(fig, ax, s, et)
    l = location(et)
    #l[3] = zero(l[3])
    d = Distribution(s, l, 2s.dx*u"km")
    d = filter(hastag(He_ipui), d)
    d.v .= Ref(rotmatX(deg2rad(-20))) .* d.v
    plot_dist(fig, ax, d)
    plot_sun(ax)
    plot_pluto(ax, l)
    plot_pepssi(ax, et)
end
=#
function plot_pepssi_view(fig::Figure, ax::PyObject, et::Float64, d::Distribution, x=location(et), R=I)
    d.v .= Ref(R) .* d.v
    plot_dist(fig, ax, d)
    plot_sun(ax)
    plot_pluto(ax, R*x)
    plot_pepssi(ax, et)
    fig, ax
end


"""
Expand x and concat y so that they may be used for a scatter plot. That is,
each y needs to be concatenated and x needs to be expanded by duplicating its
entries.

Examples:
julia> x,y = scatterargs([1,2,3], [[1,2],[1,2,3],[1]])
([1,1,2,2,2,3], [1,2,1,2,3,1])
julia> scatter(x,y)

julia> x,y,s = scatterargs([1,2,3], [[1,2],[1,2,3],[1]], [[20,10],[30,20,10],[10]])
([1,1,2,2,2,3], [1,2,1,2,3,1], [20,10,30,20,10,10])
julia> scatter(x,y, markersize=s)
"""
function scatterargs(x, eachy...)
    # Check that arguments are well formed
    # recall that eachy is a tuple of Vectors of Vectors.
    # 1. Each y in `eachy` should have the same length as `x`
    # (i.e. y contains one list for each element of `x`)
    @assert all(length.(eachy) .== length(x))

    # 2. The sequence of lengths of all the lists in each y in `eachy` must be
    # the same. (i.e. length.(y_1) == length.(y_2) for each y_i in `eachy`)
    # this sequence of lengths, `N`, is used below to expand x.
    fy, ry = firstrest(eachy)
    N = length.(fy)
    for y in ry
        @assert length.(y) == N
    end

    # The arguments should be well formed, so we expand x and concatinate the ys.
    x = expandby(x, N)
    return x, reduce.(vcat, eachy)...
end
expandby(a, n) = reduce(vcat, fill.(a,n))

"Prepare a figure with dates on the x axis"
function datefigure()
    fig, ax = plt.subplots()
    # Setting xlim to work around an error that occurs when you add a Datetime
    # locator and formatter to an axis that contains zero. (An empty axis
    # contains zero by default.)
    ax.set_xlim(2,3)
    # setting xlim turns off autoscale, so turn it back on.
    ax.autoscale()
    mdates = pyimport("matplotlib.dates")
    locator = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    ax.set_xlabel("Time")
    fig, ax
end
"Prepare a spectrogram figure"
function especfigure()
    fig, ax = datefigure()
    ax.set_ylabel("Energy (keV)")

    return fig,ax
end

function _espec_basevalues(xs, ds)
    ys = energies.(ds)
    fs = fluxes.(ds)
    scatterargs(xs, ys, fs)
end

"Returns H, xedges, yedges"
function _espec_histogram(xs, ds; kw...)
    np = pyimport("numpy")
    xs, ys, fs = _espec_basevalues(xs, ds)
    np.histogram2d(xs, ustrip.(u"keV", ys), bins=25, weights=ustrip.(u"cm^-2*s^-1",fs); kw...)
end

function plot_espec_histogram(fig, ax, xs, ds; hist_kw=Dict(), pcolor_kw=Dict())
    H,xedges,yedges = _espec_histogram(xs, ds; hist_kw...)
    mappable = ax.pcolormesh(xedges, yedges, transpose(H); pcolor_kw...)
    cbar = fig.colorbar(mappable)
    cbar.set_label("Ion flux (cm^-2 s^-1)")
    return H, xedges, yedges, mappable
end

function plot_espec_scatter(fig, ax, xs, ds)
    xs, ys, fs = _espec_basevalues(xs, ds)
    st = pyimport("spice_tools")
    ax.scatter(st.et2pydatetime.(xs), ustrip.(u"keV",ys), s=20fs./maximum(fs))
end

function plot_dist_projections(d, time, et=utc2et(time))
    fig, axs = plt.subplots(ncols=3)
    ax1,ax2,ax3 = axs

    ax1.hist2d(ustrip.(u"km/s", getindex.(d.v,1)), ustrip.(u"km/s", getindex.(d.v,2)), bins=50, weights=ustrip.(d.n))
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")

    ax2.hist2d(ustrip.(u"km/s", getindex.(d.v,1)), ustrip.(u"km/s", getindex.(d.v,3)), bins=50, weights=ustrip.(d.n))
    ax2.set_xlabel("x")
    ax2.set_ylabel("z")

    ax3.hist2d(ustrip.(u"km/s", getindex.(d.v,2)), ustrip.(u"km/s", getindex.(d.v,3)), bins=50, weights=ustrip.(d.n))
    ax3.set_xlabel("y")
    ax3.set_ylabel("z")

    fig.suptitle(time)
    for ax in axs
        ax.set_aspect("equal")
    end
    xbounds = flatten(ax.get_xlim() for ax in axs)
    xl, xr = extrema(xbounds)
    ybounds = flatten(ax.get_ylim() for ax in axs)
    yl, yr = extrema(ybounds)

    for ax in axs
        ax.set_xlim([xl, xr])
        ax.set_ylim([yl, yr])
    end
    fig, ax
end

#=
function plot_view(x)
    fig, ax = plot_dist(ipuidist(x))
    plot_sun(ax)
    plot_pluto(ax, x)
end
=#

function diff_inten_ongrid(ax, grid, d, bin, cmap=plt.get_cmap("viridis"))
    _dis = differential_intensity.(grid, Ref(bin), Ref(d))
    @assert all(length.(_dis) .== 1)
    dis = ustrip.(u"keV^-1*cm^-2*sr^-1*s^-1", getindex.(_dis, 1))
    dismissing = [di != 0 ? di : missing for di in dis]
    ldis = log10.(dismissing)
    @show M = 4#maximum(skipmissing(ldis))
    @show m = -2#minimum(skipmissing(ldis))
    normalized_ldis = (ldis .- m) ./ (M - m)
    for (square, value) in zip(grid, normalized_ldis)
        plotshape(ax, square, color="gray", linewidth=0.5, alpha=0.75)
        if value === missing
            continue
        end
        plotshape_patch(ax, square, color=cmap(value))
    end
end

polygoncolor(verts, C; kwargs...) = polygoncolor(plt.gca(), verts, C; kwargs...)
function polygoncolor(ax::PyObject, verts, C::AbstractArray{<:Number, 1}; kwargs...)
    if size(verts,1) != length(C)
        error("Size mismatch. Got $(size(verts,1)) polygons and $(length(C)) colors.")
    end
    collection = mcoll.PolyCollection(verts; kwargs...)
    collection.set_array(C)
    ax.add_collection(collection)
    ax.autoscale_view()
    collection
end

end # module
