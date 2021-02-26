using DifferentialEquations
using StaticArrays
using LinearAlgebra
using Interpolations: interpolate, extrapolate, Gridded, Linear, Flat
using Distributions: Uniform
using Random
using PhysicalConstants.CODATA2018: e, m_p, m_e
using Unitful
using Statistics
using IterTools: firstrest

includet("FieldLines.jl")
using .FieldLines

using PyCall
using PyPlot
const plt = pyimport("matplotlib.pyplot")

includet("src/KH3D_utils.jl")
using .KH3d_utils
KH3d_utils.__init__()

const pre_prefix = "/media/nathan/DATAPART11/kronos-KH/"
const step = 16
#prefix_seeded   = joinpath(pre_prefix, "seeded")
const prefix_unseeded = joinpath(pre_prefix, "unseeded")

#coord_seeded   = readcoord(prefix_seeded)
const coord_unseeded = readcoord(prefix_unseeded)

#nodes_seeded   = (coord_seeded.qx,   coord_seeded.qy,   coord_seeded.qz)
const nodes_unseeded = (coord_unseeded.qx, coord_unseeded.qy, coord_unseeded.qz)

#dom_seeded   = Domain(extrema.(nodes_seeded))
const dom_unseeded = Domain(extrema.(nodes_unseeded))

#B_seeded   = extrapolate(interpolate(nodes_seeded,   readgridvector(prefix_seeded,   "b1"), Gridded(Linear())), Flat())
const B_unseeded = extrapolate(interpolate(nodes_unseeded, readgridvector(prefix_unseeded, "b1", step), Gridded(Linear())), Flat())
#E_seeded   = extrapolate(interpolate(nodes_seeded,   readgridvector(prefix_seeded,   "E"), Gridded(Linear())), Flat())
const E_unseeded = extrapolate(interpolate(nodes_unseeded, readgridvector(prefix_unseeded, "E", step), Gridded(Linear())), Flat())

#p_seeded   = (E=E_seeded,   B=B_seeded,   dom=dom_seeded)
const _p_unseeded = (E=E_unseeded, B=B_unseeded, dom=dom_unseeded)

function find_fieldline(p; x0=SA[10000.0, 1000.0, 10000.0])
    prob = ODEProblem(norm_BB, x0, (0.0, 1e12), p, callback=domaincheck)
    solve(prob, RK4(), reltol=1e-12, abstol=1e-12)
end

#p_seeded   = (p_seeded...,   fieldline=find_fieldline(p_seeded))
const p_unseeded = (_p_unseeded..., fieldline=find_fieldline(_p_unseeded))



function plot_dist1d(dist)
    fig, ax = plt.subplots()
    xs = getindex.(dist, 2)
    vs = getindex.(dist, 1)
    ax.scatter(xs, vs, s=1, c=[(0.0, 0.0, 0.0, 0.5)])
    ax.set_xlabel("x (km)")
    ax.set_ylabel("v //")
    fig, ax
end

function EE_1d(v, x, p, t)
    x3d = p.fieldline(x)
    ustrip.(u"V/m", m_p/e * p.E(x3d...)u"km/s^2") ⋅ normalize(p.B(x3d...))
end

function electron_acceleration(v,x,p,t)
    E = EE_1d(v,x,p,t)
    ustrip(u"km/s^2", e/m_e * E*u"V/m")
end

function e_acc(a, v, x, p, t)
    x3d = p.fieldline(x)
    a .= (m_p/m_e) .* p.E(x3d...)
end

function plot_dist_with_E(dist, p)
    fig, ax = plot_dist1d(dist)
    ax2 = ax.twinx()
    ax2.plot(p.fieldline.t, EE_1d.(0, p.fieldline.t, Ref(p), 0))
    ax2.set_ylabel("E // (V/m)")
    fig, (ax, ax2)
end

function plot_dist_histogram(dist)
    xs = getindex.(dist, 2)
    vs = getindex.(dist, 1)

    a,b = extrema(xs)
    l = b-a
    a′ = 0.02 * l + a
    b′ = 0.98 * l + a
    dist′ = filter(u -> a′ ≤ u[2] ≤ b′, dist)
    vs′ = getindex.(dist′, 1)

    kT = ustrip(u"eV", m_e * var(vs*u"km/s"))
    kT′ = ustrip(u"eV", m_e * var(vs′*u"km/s"))

    fig, ax = plt.subplots()
    ax.hist([vs′, vs], histtype="stepfilled", bins="auto", label=["Excluding boundary: $(round(Int,kT′)) eV", "Including boundary: $(round(Int,kT)) eV"])

    ax.set_xlabel("v // (km/s)")
    ax.set_ylabel("f(v//) (count)")
    ax.legend()

    fig, ax
end

function spectrogram(dist)
    xs = getindex.(dist, 2)
    vs = getindex.(dist, 1)
    Es = ustrip.(u"eV", 0.5 .* m_e .* (vs.*u"km/s").^2)
    fig, ax = plt.subplots()
    ax.hist2d(xs, log.(10, Es), bins=(100,40))
    fig, ax
end

function gen_prob_func(kT, endtime)
    function prob_func(prob, i, repeat)
        p = prob.p
        a,b = extrema(p.fieldline.t)
        SecondOrderODEProblem(
            electron_acceleration,
            ustrip(u"km/s", sqrt(kT*u"eV"/m_e))*randn(),
            rand(Uniform(a,b)),
            (0.0, endtime),
            p, callback=domaincheck_1d
        )
    end
    prob_func
end
output_func(sol, i) = (sol, false)

const init_prob = SecondOrderODEProblem(electron_acceleration, 0.0, 0.5*sum(extrema(p_unseeded.fieldline.t)), (0.0,5.0), p_unseeded, callback=domaincheck_1d)

function solve_ensemble(kT, endtime, args...; kwargs...)
    eprob = EnsembleProblem(init_prob; prob_func=gen_prob_func(kT, endtime), output_func=output_func)
    solve(eprob, args...; kwargs...)
end

#e_dist_unseeded_50eV = solve(e_dist_prob_unseeded, DPRKN6(), reltol=1e-8, abstol=1e-8, dt=0.01, trajectories=10000)
function build_results(kTs, endtimes)
    pairs = Iterators.product(kTs, endtimes)
    (kT1, endtime1), pairs = firstrest(pairs)
    println("Solving for kT = ", kT1, " and endtime = ", endtime1)
    r = solve_ensemble(kT1, endtime1, DPRKN6(), reltol=1e-8, abstol=1e-8, dt=0.01, trajectories=10000)
    results = Dict((kT1, endtime1) => r)

    for (kT, endtime) in pairs
        println("Solving for kT = ", kT, " and endtime = ", endtime)
        results[(kT, endtime)] = solve_ensemble(kT, endtime, DPRKN6(), reltol=1e-8, abstol=1e-8, dt=0.01, trajectories=10000)
    end
    results
end

const results = build_results([10, 20, 50], [20.0])

function ensemblelines(N, step, xbounds, zbounds)
    p = p_at(step)
    Ux = Uniform(xbounds...)
    Uz = Uniform(zbounds...)
    pps = [(p..., fieldline=find_fieldline(p, SA[rand(Ux), 1000.0, rand(Uz)])) for i in 1:N]
end

#lines = ensemblelines(100, 16, (10000,20000), (10000,12000))
for i in 1:100
    p = (_p_unseeded..., fieldline=find_fieldline(_p_unseeded, SA[rand(Ux), 1000.0, rand(Uz)]))
    init_prob.p = p
    r = build_results([10,20,50], [20.0])
