module KH3d_utils

using PyCall
using StaticArrays
using LinearAlgebra
using Interpolations: interpolate, extrapolate, Gridded, Flat, Linear
using DifferentialEquations

export KH_Para_Dat, KH_Coord, readpara, readcoord, readgridscalar, readgridvector, Domain, domaincheck, domaincheck_1d
export BB, norm_BB, plot_background, plot_pointspread, plot_separation, plotline
export plot_B_stuff, line_projections, p_at
#export TwoWayODEProblem, twoway_solve, twoway_sol_func

const np = PyNULL()
const ff = PyNULL()
const hr = PyNULL()

function __init__()
    copy!(np, pyimport("numpy"))
    copy!(ff, pyimport("FortranFile"))
    copy!(hr, pyimport("HybridReader2").HybridReader2)
end

## Reading files stuff
struct KH_Para_Dat
    nx::Int32
    ny::Int32
    nz::Int32
    dx::Float32
    dy::Float32
    delz::Float32
end
struct KH_Coord
    nx::Int32
    ny::Int32
    nz::Int32
    qx::Vector{Float32}
    qy::Vector{Float32}
    qz::Vector{Float32}
    function KH_Coord(nx,ny,nz,qx,qy,qz)
        @assert length(qx) == nx
        @assert length(qy) == ny
        @assert length(qz) == nz
        new(nx,ny,nz,qx,qy,qz)
    end
end

readpara() = readpara(pwd())
function readpara(prefix)
    filename = joinpath(prefix, "para.dat")
    p = ff.FortranFile(filename)
    record = p.readOther([
        ("nx", np.int32),
        ("ny", np.int32),
        ("nz", np.int32),
        ("dx", np.float32),
        ("dy", np.float32),
        ("delz", np.float32)
    ])
    KH_Para_Dat(only(record)...)
end

readcoord() = readcoord(pwd())
function readcoord(prefix)
    filename = joinpath(prefix, "c.coord.dat")
    f = ff.FortranFile(filename)
    nx = only(f.readOther(np.int32))
    ny = only(f.readOther(np.int32))
    nz = only(f.readOther(np.int32))
    qx = f.readReals("f")
    qy = f.readReals("f")
    qz = f.readReals("f")
    KH_Coord(nx,ny,nz,qx,qy,qz)
end

function seeklastrecord!(file)
    index = file.index()
    file.seek(index[end-1])
end
function seektimestep!(file, step)
    index = file.index()
    # hybrid code output timesteps are a pair of records: internal step number then data.
    # Multiply `step` (which is numbered by output number) by 2 to get the record number we want.
    i = 2*step
    file.seek(index[i])
end

function _read_hybrid_datafile(file, step=nothing; index=file.index())
    i = step === nothing ? lastindex(index)-1 : 2*step
    file.seek(index[i])
    file.readReals()
end

readgridscalar(name) = readgridscalar(pwd(), name)
function readgridscalar(prefix, name, step=nothing)
    file = ff.FortranFile(joinpath(prefix, "c.$(name).dat"))
    raw = _read_hybrid_datafile(file, step)

    para = readpara(prefix)
    reshape(raw, (para.nx, para.ny, para.nz))
end

readgridvector(name) = readgridvector(pwd(), name)
function readgridvector(prefix, name, step=nothing)
    file = ff.FortranFile(joinpath(prefix, "c.$(name).dat"))
    raw = _read_hybrid_datafile(file, step)

    para = readpara(prefix)
    nx,ny,nz = para.nx, para.ny, para.nz
    arr = reshape(raw, (nx, ny, nz, 3))
    ret = Array{SVector{3,eltype(arr)},3}(undef, nx, ny, nz)
    for k in 1:nz, j in 1:ny, i in 1:nx
        ret[i,j,k] = SVector{3}(@view arr[i,j,k,:])
    end
    return ret
end

## Field line stuff

struct Domain{T}
    min_x::T
    min_y::T
    min_z::T
    max_x::T
    max_y::T
    max_z::T
end
Domain(a,b,c,d,e,f) = Domain(promote(a,b,c,d,e,f)...)
Domain(boundaries::Tuple{T,T,T}) where T<:Tuple{U,U} where U = Domain(
    boundaries[1][1], boundaries[2][1], boundaries[3][1],
    boundaries[1][2], boundaries[2][2], boundaries[3][2]
)
Base.in(x, d::Domain) = d.min_x <= x[1] <= d.max_x && d.min_y <= x[2] <= d.max_y && d.min_z <= x[3] <= d.max_z

function BB(x,p,t)
    p.B(x...)
end

function norm_BB(x, p, t)
    normalize(BB(x,p,t))
end

const domaincheck = DiscreteCallback(
    (x,t,integrator) -> x ∉ integrator.p.dom,
    #integrator -> terminate!(integrator, :Success),
    terminate!,
    save_positions=(true,false)
)

const domaincheck_1d = DiscreteCallback(
    (x,t,integrator) -> x.x[2]<0 || x.x[2]>integrator.p.fieldline.t[end],
    terminate!,
    save_positions=(true, false)
)


function plotline(ax, u_i, args...; kwargs...)
    prob = ODEProblem(norm_BB, u_i, (0.0, 1e12), p, callback=domaincheck)
    sol = solve(prob, args...; kwargs...)

    plotline(ax, sol)
end
function plotline(ax, line::ODESolution, vars=(1,3); kwargs...)
    ts = range(extrema(line.t)...; length=1000)
    ax.plot(getindex.(line.(ts), vars[1]), getindex.(line.(ts), vars[2]); kwargs...)
end

plot_background(prefix) = plot_background(plt.subplots()..., prefix)
function plot_background(fig, ax, prefix)
    coord = readcoord(prefix)
    mixed = readgridscalar(prefix, "mixed")
    mesh = ax.pcolormesh(coord.qx, coord.qz, transpose(@view mixed[:,size(mixed,2)÷2,:]))
    return fig, ax, mesh
end

function plot_pointspread(ax, p; starting_point=SA[14195.786, 718.774 + 450, 10961.3], init_spread=0.5)
    B_prob = ODEProblem(norm_BB, starting_point, (0.0, 1e12), p, callback=domaincheck)
    #fieldline = solve(prob, Vern8(), reltol=1e-12, abstol=1e-12)

    prob_func(prob, i, repeat) = ODEProblem(norm_BB, init_spread*randn(3)+starting_point, (0.0, 1e12), p, callback=domaincheck)
    output_func(sol, i) = ((sol[begin], sol[end]), false)

    nearbylines_prob = EnsembleProblem(B_prob; output_func, prob_func)
    nearbylines_sol = solve(nearbylines_prob, RK4(), EnsembleSerial(), reltol=1e-12, abstol=1e-12, trajectories=5000)
    utop = filter(u->u[2][2]>p.B.itp.knots[2][1]+50, nearbylines_sol.u)
    starts = getindex.(utop, 1)
    ends = getindex.(utop,2)

    ax.plot(getindex.(starts,1), getindex.(starts, 3), color="white", marker=".", markersize=1.0, linestyle="None")
    ax.plot(getindex.(ends,1), getindex.(ends, 3), color="black", marker=".", markersize=1.0, linestyle="None")
    nothing
end

function trace_retrace(alg, p; trace_init=starting_point, solution_kwargs...)
    trace_prob = ODEProblem(norm_BB, trace_init, (0.0, 1e12), p, callback=domaincheck)
    trace = solve(trace_prob, alg; solution_kwargs...)
    retrace_init = trace[end-1]
    retrace_prob = ODEProblem(norm_BB, retrace_init, (1e12, 0.0), p, callback=domaincheck)
    retrace = solve(retrace_prob, alg; solution_kwargs...)
    endpoints = (trace[1], trace[end-1], retrace[1], retrace[end-1])
    (getindex.(endpoints, 1), getindex.(endpoints, 3))
end
function plot_separation(ax, p, trace_init)
    x,y = trace_retrace(RK4(), p; trace_init, reltol=1e-12, abstol=1e-12)
    ax.plot(x[1], y[1], marker="o", color="red", markersize=4)
    ax.plot(x[2], y[2], marker="x", color="black", markersize=4)
    @assert x[2] == x[3]
    @assert y[2] == y[3]
    ax.plot(x[4], y[4], marker="o", color="blue", markersize=4)
end

function p_at(step)
    (
        B   = extrapolate(interpolate(nodes_unseeded, readgridvector(prefix_unseeded, "b1", step), Gridded(Linear())), Flat()),
        E   = extrapolate(interpolate(nodes_unseeded, readgridvector(prefix_unseeded, "E", step), Gridded(Linear())), Flat()),
        dom = dom_unseeded
    )
end

function plot_B_stuff(B)
    cb_args = (fraction=0.08, pad=0.0, aspect=40)

    fig, axs = plt.subplots(ncols=3, nrows=3, sharex="row", sharey="row")#, subplot_kw=Dict(:aspect=>"equal"))
    ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = eachrow(axs)
    m1 = ax1.pcolormesh(qx, qy, transpose( getindex.(B[:,:, size(B,3)÷2], 1) ), vmin=-0.4, vmax=0.4)
    m2 = ax2.pcolormesh(qx, qy, transpose( getindex.(B[:,:, size(B,3)÷2], 2) ), vmin=0.0, vmax=0.8)
    m3 = ax3.pcolormesh(qx, qy, transpose( getindex.(B[:,:, size(B,3)÷2], 3) ), vmin=-0.4, vmax=0.4)
    fig.colorbar(m1, ax=ax1; label="B_x", cb_args...)
    fig.colorbar(m2, ax=ax2; label="B_y", cb_args...)
    fig.colorbar(m3, ax=ax3; label="B_z", cb_args...)
    ax1.set_ylabel("Y")
    ax1.set_xlabel("X")
    ax2.set_xlabel("X")
    ax3.set_xlabel("X")

    m4 = ax4.pcolormesh(qz, qy, getindex.(B[size(B,1)÷2, :, :], 1), vmin=-0.4, vmax=0.4)
    m5 = ax5.pcolormesh(qz, qy, getindex.(B[size(B,1)÷2, :, :], 2), vmin=0.0, vmax=0.8)
    m6 = ax6.pcolormesh(qz, qy, getindex.(B[size(B,1)÷2, :, :], 3), vmin=-0.4, vmax=0.4)
    fig.colorbar(m4, ax=ax4; label="B_x", cb_args...)
    fig.colorbar(m5, ax=ax5; label="B_y", cb_args...)
    fig.colorbar(m6, ax=ax6; label="B_z", cb_args...)
    ax4.set_ylabel("Y")
    ax4.set_xlabel("Z")
    ax5.set_xlabel("Z")
    ax6.set_xlabel("Z")

    m7 = ax7.pcolormesh(qx, qz, transpose(getindex.(B[:, size(B,2)÷2, :], 1)), vmin=-0.4, vmax=0.4)
    m8 = ax8.pcolormesh(qx, qz, transpose(getindex.(B[:, size(B,2)÷2, :], 2)), vmin=0.0, vmax=0.8)
    m9 = ax9.pcolormesh(qx, qz, transpose(getindex.(B[:, size(B,2)÷2, :], 3)), vmin=-0.4, vmax=0.4)
    fig.colorbar(m7, ax=ax7; label="B_x", cb_args...)
    fig.colorbar(m8, ax=ax8; label="B_y", cb_args...)
    fig.colorbar(m9, ax=ax9; label="B_z", cb_args...)
    ax7.set_ylabel("Z")
    ax7.set_xlabel("X")
    ax8.set_xlabel("X")
    ax9.set_xlabel("X")

    fig.subplots_adjust(top=.98, hspace=0.4)

    fig, axs
end

function line_projections(step; N_lines=100, scale=1.0)
    fig, axs = plot_B_stuff(readgridvector(prefix_unseeded, "b1", step))
    fls = find_fieldline.(Ref(p_at(step)), [SA[(qx[1]+qx[end])/2, 1000.0, (qz[1]+qz[end])/2] + scale*randn(3) for i in 1:N_lines])
    for fl in fls
        for (row, vars) in zip(eachrow(axs), [(1,2),(3,2),(1,3)])
            for ax in row
                plotline(ax, fl, vars; color="black")
            end
        end
    end
    fig, axs
end

function find_fieldline(p, x0=SA[10000.0, 1000.0, 10000.0])
    prob = ODEProblem(norm_BB, x0, (0.0, 1e12), p, callback=domaincheck)
    solve(prob, RK4(), reltol=1e-8, abstol=1e-8)
end

end # module
