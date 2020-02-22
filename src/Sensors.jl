module Sensors
export density, bulkvelocity, pressuretensor, pressure, thermalenergy, energyspectrogram, flythrough

using StaticArrays
using LinearAlgebra
using StatsBase
using Unitful

using ..Simulations
using ..Spacecraft
using ..SphericalPolygons

# Get moments from a distribution
density(f::Distribution) = sum(f.ns)
bulkvelocity(f::Distribution) = sum(f.vs .* f.ns)/density(f)
function pressuretensor(f::Distribution)
    uf = bulkvelocity(f)
    sum(zip(f.vs, f.ns, f.ms)) do (v,n,m)
        w = v-uf
        m * w*w' * n
    end
end
pressure(f::Distribution) = tr(pressuretensor(f))/3
"Also called kinetic temperature"
thermalenergy(f::Distribution) = pressure(f)/density(f)

"Compute E/q spectrum for the given distribution"
function energyspectrogram(f::Distribution) where T
    Eoverq = (1/2 .* f.ms .* norm.(f.vs).^2) ./ f.qs
    Eoverq_ = ustrip.(u"J*C^-1", Eoverq)
    G = 0.15u"cm^2*sr" * 10u"s" / (4π*u"sr") # hypothetical 4π sr particle detector
    counts = ustrip.(Unitful.NoUnits, G .* norm.(f.vs) .* f.ns)
    fit(Histogram, Eoverq_, weights(log10.(counts)), exp10.(range(1, stop=4, length=64)))
end

flythrough(s::Simulation, traj::Trajectory, sensor, sensor_args...) = flythrough(s, traj.pos, sensor, sensor_args...)
function flythrough(s::Simulation, traj_arr::AbstractVector, sensor, sensor_args...)
    ds = Distribution.(s, traj_arr)
    T = typeof(sensor(ds[1],sensor_args...))
    Vector{Union{T,Nothing}}[length(d)==0 ? nothing : sensor(d,sensor_args...) for d in ds]
end

end # module
