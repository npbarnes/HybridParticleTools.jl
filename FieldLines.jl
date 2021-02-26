module FieldLines
using DifferentialEquations
using LinearAlgebra
using Interpolations
using Unitful
using Distributions
using QuadGK
using PhysicalConstants.CODATA2018: e, m_p, m_e
const e_C = ustrip(u"C", e)
const mp_kg = ustrip(u"kg", m_p)
const me_kg = ustrip(u"kg", m_e)

export find_fieldline, make_electron_ensemble, solve_electron_ensemble, make_electron_problem, interpolate_E_parallel, eval_on_fieldline, loopintegral

function in_domain(x,vf)
    vf.zgrid.nodes[1][begin] ≤ x[1] ≤ vf.zgrid.nodes[1][end] &&
    vf.zgrid.nodes[2][begin] ≤ x[2] ≤ vf.zgrid.nodes[2][end] &&
    vf.zgrid.nodes[3][begin] ≤ x[3] ≤ vf.zgrid.nodes[3][end]
end

const domaincheck = DiscreteCallback(
    (x,t,integrator) -> !in_domain(integrator.u, integrator.p.B),
    terminate!,
    save_positions=(false,false)
)


function norm_B(x,p,t)
    normalize(p.B(x[1],x[2],x[3]))
end

function find_fieldline(B, x0)
    prob = ODEProblem(norm_B, x0, (0.0, 1e12), (;B), callback=domaincheck)
    solve(prob, RK4(), reltol=1e-8, abstol=1e-8, save_end=false)
end

function eval_on_fieldline(fieldline, x_1d, F)
    [F(x3d[1], x3d[2], x3d[3]) for x3d in fieldline.(x_1d)]
end

function interpolate_E_parallel(fieldline, x_1d, E, B)
    E_vec_1d_arr = eval_on_fieldline(fieldline, x_1d, E)
    B_vec_1d_arr = eval_on_fieldline(fieldline, x_1d, B)
    E_parallel_1d_arr = dot.(E_vec_1d_arr, normalize.(B_vec_1d_arr))
    extrapolate(interpolate((x_1d,), E_parallel_1d_arr, Gridded(Linear())), Flat())
end

function interpolate_1d(line, field, x_1d=range(line.t[begin], line.t[end], length=10length(line.t)))
    f1d = eval_on_fieldline(line, x_1d, field)
    extrapolate(interpolate((x_1d,), f1d, Gridded(Linear())), Flat())
end


function electron_acceleration(v,x,p,t)
    # This should come out in km/s^2 if E is in V/km
    - e_C/me_kg * (p.E(x)/1000)/1000
end

function make_electron_problem(fieldline, E, B)
    a = fieldline.t[begin]
    b = fieldline.t[end]
    x_1d = range(a, b, length=10*length(fieldline.t))
    E_parallel = interpolate_E_parallel(fieldline, x_1d, E, B)

    domaincheck_1d = DiscreteCallback(
        (x,t,integrator) -> !(a ≤ x[2] ≤ b),
        terminate!,
        save_positions = (false, false)
    )
    SecondOrderODEProblem(electron_acceleration, 0.0, -1.0, (0.0, 1.0), (;E=E_parallel), callback=domaincheck_1d)
end

function gen_prob_func(kT)
    function prob_func(prob, i, repeat)
        a,b = extrema(prob.p.E.itp.knots[1])
        remake(prob; u0=ArrayPartition(
            ustrip(u"km/s", sqrt(kT*u"eV"/m_e))*randn(),
            rand(Uniform(a,b))
            )
        )
    end
    prob_func
end
output_func = (sol,i) -> (sol, false)

function make_electron_ensemble(E, B, x0, kT)
    fl = find_fieldline(B, x0)
    init_prob = make_electron_problem(fl, E, B)
    prob_func = gen_prob_func(kT)
    EnsembleProblem(init_prob; prob_func, output_func)
end

function solve_electron_ensemble(ensemble_prob)
    solve(ensemble_prob, DPRKN6(), reltol=1e-8, abstol=1e-8, trajectories=10000)
end

function loopintegral(E, x, yl, yr, zbot, ztop)
    Iab, err_ab = quadgk(y->E(x, y, zbot)[2], yl, yr)
    Ibc, err_bc = quadgk(z->E(x, yr, z)[3], zbot, ztop)
    Icd, err_cd = quadgk(y->E(x, y, ztop)[2], yr, yl)
    Ida, err_da = quadgk(z->E(x, yl, z)[3], ztop, zbot)
    Iab + Ibc + Icd + Ida
end



end # module
