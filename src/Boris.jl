module Boris

export boris, borisstep

using LinearAlgebra
using ..Distributions
import PhysicalConstants.CODATA2018: m_p, e

function vstep!(v,dt,E,B)
    v .= vstep(v,dt,E,B)
end
function vstep(v,dt,E,B)
    t = B*dt/2
    s = 2t/(1+t⋅t)

    v⁻ = v + E*dt/2
    v′ = v⁻ + v⁻ × t
    v⁺ = v⁻ + v′ × s
    return v⁺ + E*dt/2
end
function xstep!(x,v)
    x .= xstep(x,v)
end
function xstep(x,dt,v)
    x + v*dt
end

function borisstep(x, v, dt, E, B)
    v = next_v(v,dt,E,B)
    x = next_x(x,v)
    return x,v
end

function boris(xinit, vinit, dt, E, B, N;m=m_p,q=e)
    ret_x = Vector{typeof(xinit)}(undef, N+1)
    ret_v = Vector{typeof(vinit)}(undef, N+1)
    ret_x[1] = xinit
    ret_v[1] = vinit
    qm = q/m
    @inbounds for i in 1:N
        EE = qm*E(ret_x[i][1], ret_x[i][2], ret_x[i][3])
        BB = qm*B(ret_x[i][1], ret_x[i][2], ret_x[i][3])
        ret_v[i+1] = vstep(ret_v[i], dt, EE, BB)
        ret_x[i+1] = xstep(ret_x[i], dt, ret_v[i+1])
    end
    return ret_x, ret_v
end

function boris(d::Distribution, xinit, dt, E, B, N)
    broadcast((x,v,dt,E,B,N,m,q) -> boris(x,v,dt,E,B,N,m=m,q=q),
        Ref(xinit),
        d.v,
        dt,
        (x,y,z)->E(x,y,z),
        (x,y,z)->B(x,y,z),
        N, d.m, d.q)
end

end # module
