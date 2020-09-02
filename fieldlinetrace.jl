include("./working.jl")

using DifferentialEquations

function BB(v,x,B,t)
    v .= ustrip.(u"T", B((x*u"km")...))
end

prob = ODEProblem(BB, [13.0*1187, 0.0, 0.0], (0.0, 1_000_000_000_000_000))
sol = solve(prob)
