module Boris

export boris, borisstep

using LinearAlgebra
import PhysicalConstants.CODATA2018: m_p, e

function vstep(v,dt,E,B)
    t = B*dt/2
    s = 2t/(1+t⋅t)

    v⁻ = v + E*dt/2
    v′ = v⁻ + v⁻ × t
    v⁺ = v⁻ + v′ × s
    return v⁺ + E*dt/2
end

function xstep(x,dt,v)
    x + v*dt
end

function boris(xinit, vinit, dt, E, B, N, m=m_p, q=e)
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

end # module
