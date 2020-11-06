include("working.jl")
finaltime = 25000*dt
step = 20dt
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

function plot_ppc(ax, kd, ets)
    ax.plot(ets, length.(inrange.(Ref(kd), [ustrip.(u"km", l) for l in location.(ets)], 500)))
end
