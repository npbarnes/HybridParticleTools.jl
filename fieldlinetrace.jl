include("./working.jl")

using DifferentialEquations

struct MyODEParameters{T,U}
    B::T
    dom::U
end
p = MyODEParameters(B, dom)

function BB(v,x,p,t)
    v .= ustrip.(u"nT", p.B((x*u"km")...))
end

domaincheck = DiscreteCallback(
    (x,t,integrator) -> x*u"km" ∉ integrator.p.dom,
    integrator -> terminate!(integrator, :Success),
    save_positions=(true,false)
)

struct TwoWayODEProblem
    p1::ODEProblem
    p2::ODEProblem
    function TwoWayODEProblem(f, u_i, tspan, p=NullParameters(), t_i=mean(tspan); callback=nothing)
        @assert tspan[1] <= t_i <= tspan[2]
        p1 = ODEProblem(f, u_i, (t_i, tspan[1]), p, callback=callback)
        p2 = ODEProblem(f, u_i, (t_i, tspan[2]), p, callback=callback)
        new(p1, p2)
    end
end

function twoway_solve(prob::TwoWayODEProblem, args...; kwargs...)
    s1 = solve(prob.p1, args...; kwargs...)
    s2 = solve(prob.p2, args...; kwargs...)
    (s1,s2)
end

prob = TwoWayODEProblem(BB, [13.0*1187, 0.0, 0.0], (-10^12, 10^12), p, 0.0, callback=domaincheck)
(s1,s2) = twoway_solve(prob, Tsit5())

sol(t) = t<0 ? s1(t) : s2(t)
ts = collect(range(minimum(s1.t), maximum(s2.t), length=1000))
plt.plot(ts, sol.(ts))
