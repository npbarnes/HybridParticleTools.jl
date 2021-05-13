include("use_stuff.jl")
include("samplers.jl")
include("chex_defs.jl")
include("dipoles.jl")
B_dip = Dipole(SA[0.0,0.0, 4π*ustrip(u"km", 1u"Rp")^3*ustrip(u"T", 30u"nT")/ustrip(u"kg*km*s^-2*A^-2", μ_0)])

dom = Boris.Domain(-100.0u"Rp", -50u"Rp", -50u"Rp", 20u"Rp", 50u"Rp", 50u"Rp")

function test_particle_simulation(path, fields,  dom, sampler, σ, nn)
    L = dom.max_x - dom.min_x
    transittime = L/400u"km/s"
    #gyroperiod = 2pi*4m_p/(e*0.1u"nT")
    finaltime = 2transittime
    timestep = 0.333u"s"
    N = ceil(Int, finaltime/timestep)
    dx = ([sampler() for i in 1:1000000] |> a->getindex.(a,1) |> a->quantile(a,0.99)) * timestep
    pg! = UnBufferedGenerator(sampler, dom.max_x, dx, (dom.min_y,dom.max_y), (dom.min_z,dom.max_z), 3500, timestep)

    ps = trace_particles(pg!, fields, dom, timestep, N, σ, nn)
    serialize(path, ps)
    ps
end
base_path = "/media/nathan/DATAPART11/testparticles"

high_imf = "/home/nathan/data/pre-2019/2017-Mon-Nov-13/pluto-7/data"
medium_imf = "/home/nathan/data/pre-2019/2018-Fri-Jan-26/pluto-1/data"
low_imf = "/home/nathan/data/pre-2019/2017-Mon-Nov-13/pluto-7/data"

B_hyb = loadB(high_imf)
B = B_hyb + Covariant_fromfunction(B_hyb.maingrid, B_dip)
E = loadE(high_imf; B)
f = Boris.Fields(E,B)
for sampler in [:shell, :superthermal]
    (σ, nn) = (:σ_standard, :N_power_law)
    name = join([sampler, nn, σ], "_")
    filename = "$name.jls"
    @eval test_particle_simulation(joinpath(base_path, "dipole_added", "high_IMF", $filename), f, dom, $sampler, $σ, $nn)
end

B_hyb = loadB(medium_imf)
B = B_hyb + Covariant_fromfunction(B_hyb.maingrid, B_dip)
E = loadE(medium_imf; B)
f = Boris.Fields(E,B)
for sampler in [:shell, :superthermal]
    (σ, nn) = (:σ_standard, :N_power_law)
    name = join([sampler, nn, σ], "_")
    filename = "$name.jls"
    @eval test_particle_simulation(joinpath(base_path, "dipole_added", "medium_IMF", $filename), f, dom, $sampler, $σ, $nn)
end

B_hyb = loadB(low_imf)
B = B_hyb + Covariant_fromfunction(B_hyb.maingrid, B_dip)
E = loadE(low_imf; B)
f = Boris.Fields(E,B)
for sampler in [:shell, :superthermal]
    (σ, nn) = (:σ_standard, :N_power_law)
    name = join([sampler, nn, σ], "_")
    filename = "$name.jls"
    @eval test_particle_simulation(joinpath(base_path, "dipole_added", "low_IMF", $filename), f, dom, $sampler, $σ, $nn)
end

#=
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
#fig, ax = especfigure();
#PlottingTools.plot_espec_scatter(fig, ax, flyby_ets, filter.(pepssifov, ds));

function plot_ppc(ax, kd, ets)
    ax.plot(ets, length.(inrange.(Ref(kd), [ustrip.(u"km", l) for l in location.(ets)], 500)))
end


function diff_inten_by_angle(d, θs=range(0,pi, length=2000); bin_edges=pepssi_bin_edges, cone_angle=0.0)
    fovs = [SCircle([cos(θ); rotmat(-cone_angle)*[0.0, sin(θ)]], acos(1-pepssi_S0_area/(2pi))) for θ in θs]
    dis = differential_intensity.(fovs, Ref(bin_edges), Ref(d))
    θs, dis
end

function sample_dist(p=5, n=4.7e10u"km^-3", l=Int(1e6))
    Distribution(
        [superthermal(p) for i in 1:l],
        fill(4m_p, l),
        fill(1e, l),
        fill(n/l, l),
        fill(He_ipui, l)
    )
end

function plot_diff_inten_by_angle(θs, dis; bins=1:length(pepssi_bin_edges)-1)
    fig, ax = plt.subplots()
    for i in bins
        ax.plot(rad2deg.(θs), ustrip.(u"keV^-1*cm^-2*s^-1*sr^-1", getindex.(dis,i)), label="$(pepssi_bin_edges[i]) - $(pepssi_bin_edges[i+1]) ($(pepssi_bin_names[i]))")
    end

    ax.set_yscale("log")
    ax.set_xlim([0,180])
    ticker = pyimport("matplotlib.ticker")
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))

    ax.plot([30, 40, 50], [40, 13, 1], color="black", label="uncertainty (eyeballed from Kollmann et al. 2019 Figure 1)")
    ax.plot([40, 50, 60], [500, 230, 20], color="black")

    ax.plot([30, 40], [500,500], color="black", linestyle="--")
    ax.plot([50, 60], [1, 1], color="black", linestyle="--")

    ax.set_ylabel("Differential Intensity\n(1/(keV cm\$^2\$ s sr)")
    ax.set_xlabel("Cone angle")
    fig.legend()
    fig, ax
end
=#
