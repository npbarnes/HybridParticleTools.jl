using HybridTools
using Dates
using PyCall
import PyPlot
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

const np = pyimport("numpy")
const st = pyimport("spice_tools")
const s = Simulation(pwd(), n=12)

const pepssifov = fov_polygon("NH_PEPSSI_S0", st.nh_in_wake)
const swapfov = fov_polygon("NH_SWAP", st.nh_in_wake)

const et_start = utc2et"11:30"
const et_end   =  utc2et"13:59"
const ets = collect(et_start:30:et_end)

ipuidist(d::Distribution) = filter(x->x.t==He_ipui, d)
ipuidist(pos::AbstractVector) = ipuidist(Distribution(s,pos))
ipuidist(et::Number) = ipuidist(location(et))

plot_distmap(arg) = plot_distmap(mapfigure()..., arg)
function plot_distmap(fig, ax, et::Number)
    plot_distmap(fig, ax, location(et))
    plot_pepssi(ax, et)
    fig, ax
end
function plot_distmap(fig, ax, pos::Vector)
    plot_distmap(fig, ax, ipuidist(pos))
    plot_pluto(ax,pos)
    coord = round.(pos./1189; digits=1)
    ax.set_title("Distribution at coordinates:\n$(repr(coord))")
    fig, ax
end
function plot_distmap(fig, ax, d::Distribution)
    plotdist(fig, ax, d)
    plot_sun(ax)
    fig, ax
end

# Coordinates are always numbers.
_as_xcoords(ets::AbstractVector{<:Number}) = ets
_as_xcoords(locations::AbstractVector{<:Vector}) = getindex.(locations, 1)
_as_xcoords(arg) = 1:length(ds)

# an axis is usually numbers, but may be something else plottable with matplotlib
# e.g. DateTime objects
_as_xaxis(ets::AbstractVector{<:Number}) = st.et2pydatetime.(ets)
_as_xaxis(arg) = _ax_xcoords(arg)

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

function especfigure()
    fig, ax = plt.subplots()
    # Setting xlim to work around an error that occurs when you add a Datetime
    # locator and formatter to an axis that contains zero. (An empty axis
    # contains zero by default.) Turn autoscale back on after.
    ax.set_xlim(2,3)
    ax.autoscale()
    mdates = pyimport("matplotlib.dates")
    locator = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    ax.set_xlabel("Time")
    ax.set_ylabel("Energy (keV)")

    return fig,ax
end

plot_espec_histogram(xs, fovs=fov_polygon.("NH_PEPSSI_S0", ets)) = plot_espec_histogram(especfigure()..., xs, fovs)
function plot_espec_histogram(fig, ax, xs, fovs=fov_polygon.("NH_PEPSSI_S0", ets))
    ds = filter.(fovs, ipuidist.(xs))
    ys = energies.(ds)
    Ns = getproperty.(ds, :n)

    h,xedges,yedges,mappable = _histogram2d(ax, scatterargs(xs, ys, Ns)...)
    cbar = fig.colorbar(mappable)
    cbar.set_label("Apparent ion density")

    fig, ax
end

"""
Most of the time scatterargs() can be used to generate arguments for plt.hist2d()
However DateTime objects don't play nice with plt.hist2d(), so we use this helper.
The helper function _as_xcoords() and _as_xaxis() are used to convert first to a
coordinate that can be histogrammed, and then to an axis for plotting.
"""
function _histogram2d(ax,xs,ys,Ns; hist_kw=Dict(), pcolor_kw=Dict())
    xcoords = _as_xcoords(xs)
    H,xedges_coord,yedges = np.histogram2d(xcoords, ustrip.(u"keV", ys), bins=25, weights=ustrip.(u"km^-3",Ns); hist_kw...)
    xedges_axis = _as_xaxis(xedges_coord)
    artist = ax.pcolormesh(xedges_axis, yedges, transpose(H); pcolor_kw...)
    return H, xedges_axis, yedges, artist
end

function plot_espec_scatter(fig, ax, xs, fovs=fov_polygon.("NH_PEPSSI_S0", ets))
    # get the distributions
    ds = filter.(fovs, ipuidist.(xs))

    # compute the relavent quantities
    ys = energies.(ds)
    ns = getproperty.(ds, :n)
    vs = [norm.(vs) for vs in getproperty.(ds,:v)]

    # expand to uniform arrays
    xs, ys, ns, vs = scatterargs(xs, ys, ns, vs)

    # convert to something appropriate for plotting
    xs = _as_xaxis(xs)
    ys = ustrip.(u"keV", ys)
    fs = ustrip.(u"s^-1*cm^-2", ns .* vs)

    # plot
    ax.scatter(xs, ys, s=20fs/maximum(fs))
    ax.set_xlim(extrema(xs)...)
    ax.set_ylim(extrema(ys)...)

    fig, ax
end

#const prefix = joinpath(pwd(),"data")
#const para = ParameterSet(prefix)
#const dt = Float64(para.dt)*u"s"
#const xinit = SA[8.0, 0.0, 0.0]u"Rp"
#const E,B = loadfields(prefix, 27)
#const d = ipuidist(xinit)

nothing
