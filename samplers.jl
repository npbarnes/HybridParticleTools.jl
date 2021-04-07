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

"""UnBufferedGenerators are for sampling N particles at a time from the flux of a distribution sampled by v_sampler."""
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

"""Helper function to figure out the number of micro particles represented by each macro particle.
Check the answer, though. I think there's something wrong with it.
"""
function macroparticle_N(ubg::UnBufferedGenerator, trueflux)
    A = (ubg.y_extent[2] - ubg.y_extent[1])*(ubg.z_extent[2] - ubg.z_extent[1])
    macroflux = ubg.N/(A*ubg.dt)
    uconvert(Unitful.NoUnits, trueflux/macroflux)
end
