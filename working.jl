include("use_stuff.jl")

# approx bins for PEPSSI Helium (I assume that means signly ionized):
const pepssi_bin_names = ["L13", "L11", "L09", "L07", "L05", "L03", "L01"]
const pepssi_bin_edges = [2.2, 4.79, 9.485, 18.05, 33.4, 59.8, 106, 183]u"keV"
const pepssi_S0_area = 0.092677

#folder = "/media/nathan/DATAPART11/2020-Thu-Jan-23/pluto-2"
#folder = "/media/nathan/DATAPART11/2020-Mon-Aug-10/pluto-1"
#folder = "/home/nathan/data/pre-2019/2017-Mon-Nov-13/pluto-7"
#const np = pyimport("numpy")
#const st = pyimport("spice_tools")
#const s = Simulation(folder, n=6)

#const pepssifov = fov_polygon("NH_PEPSSI_S0", st.nh_in_wake)
#const swapfov = fov_polygon("NH_SWAP", st.nh_in_wake)

#const et_start = utc2et"11:30"
#const et_end   =  utc2et"14:01"
#ets = collect(et_start:(10*60):et_end)

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

ch4dist(d::Distribution) = filter(x->x.t==CH4_chex || x.t==CH4_photo || x.t==CH4_stagnant, d)
ch4dist(pos::AbstractVector) = ch4dist(Distribution(s,pos, s.dx*u"km"))
ch4dist(et::Number) = ch4dist(location(et))


#prefix = joinpath(folder,"data")
#const para = ParameterSet(prefix)
#const dt = Float64(para.dt)*u"s"
#E,B = loadfields(prefix)
#f = Boris.Fields(E,B)
#const dom = Boris.Domain(Tuple(t.*u"km" for t in extrema.(E.xgrid.nodes)))

function todaysplot(time::String)
    et = utc2et(time)
    d = _traced_distribution(ps, kd, inv(R)*location(et), 3*700u"km", macro_N)
    d.v .= Ref(R) .* d.v
    fig, ax = mapfigure(false)
    diff_inten_ongrid(ax, gs, d, [9.5, 18.0]u"keV")
    ax.set_title(time)
    plot_sun(ax)
    plot_pluto(ax, location(et))
end

nothing
