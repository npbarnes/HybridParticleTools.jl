module WorkingModule
export plt, currentsim, issimnull, dist, ipuidist, He_ipuidist, H_ipuidist, pepssi_bin_edges, pepssi_bin_names, pepssi_S0_area

using Rexport

@reexport using HybridTools
@reexport using Dates
@reexport using PyCall
@reexport using PyPlot
@reexport using StaticArrays
@reexport using Unitful
@reexport using LinearAlgebra
@reexport using Statistics
@reexport using NearestNeighbors
@reexport using Random
@reexport using Base.Threads
using IterTools: firstrest
using Distributions: Pareto
@reexport import PhysicalConstants.CODATA2018: m_p, e
const plt = PyNULL()

function __init__()
    copy!(plt, pyimport("matplotlib.pyplot"))
end

# Save a current working simulation that will be used with some helper methods defined below.
mutable struct CurrentSimulation
    s::Union{Simulation, Nothing}
    E::Union{VectorField, Nothing}
    B::Union{VectorField, Nothing}
end
const CURRENT_SIMULATION = CurrentSimulation(nothing, nothing, nothing)

function issimnull()
    props = propertynames(CURRENT_SIMULATION)
    nullities = [getproperty(CURRENT_SIMULATION, p) === nothing for p in props]
    if all(nullities)
        return true
    elseif any(nullities)
        return false
    else
        nulls = [p for (i,p) in enumerate(props) if nullities[i] === nothing]
        error("Either all the fields should be nothing or none of them. Only property/ies $(nulls) were found to be nothing.")
    end
end

function currentsim()
    if issimnull()
        error("No current simulation selected")
    end
    CURRENT_SIMULATION.s
end
function currentfields()
    if issimnull()
        error("No current simulation selected")
    end
    Boris.Fields(CURRENT_SIMULATION.E, CURRENT_SIMULATION.B)
end

function currentsim(path, n=6, step=-1)
    s = Simulation(path; n, step)
    E, B = loadfields(joinpath(path, "data"); step)
    currentsim(s, E, B)
end
function currentsim(s, E, B)
    CURRENT_SIMULATION.s = s
    CURRENT_SIMULATION.E = E
    CURRENT_SIMULATION.B = B
    CURRENT_SIMULATION
end

dist(args..., kwargs...) = dist(currentsim(), args...; kwargs...)
dist(s::Simulation, pos::AbstractVector, dx=s.dx*u"km") = Distribution(s, pos, dx)
dist(s::Simulation, et::Number) = dist(s, location(et))
dist(s::Simulation, time::String) = dist(s, utc2et(time))

ipuidist(args...; kwargs...) = filter(hasanytag(H_ipui, He_ipui), dist(args...; kwargs...))
H_ipuidist(args...; kwargs...) = filter(hastag(H_ipui), dist(args...; kwargs...))
He_ipuidist(args...; kwargs...) = filter(hastag(He_ipui), dist(args...; kwargs...))

# approx bins for PEPSSI Helium (I assume that means signly ionized):
const pepssi_bin_names = ["L13", "L11", "L09", "L07", "L05", "L03", "L01"]
const pepssi_bin_edges = [2.2, 4.79, 9.485, 18.05, 33.4, 59.8, 106, 183]u"keV"
const pepssi_S0_area = 0.092677

plot_view(x) = plot_view(currentsim(), x)
function plot_view(s::Simulation, x)
    fig, ax = plot_dist(ipuidist(s, x))
    plot_sun(ax)
    plot_pluto(ax, x)
end
end
