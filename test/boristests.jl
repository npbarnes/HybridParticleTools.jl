using Test

using LinearAlgebra
using StaticArrays
using Plots

using HybridParticleTools.Boris
using HybridParticleTools.Utility

"x should trace out a circle in xy coordinates, z should be zero"
function _estimate_gyroradius(x)
    x = viewasarray(x)
    ex = extrema(x, dims=2)
    radii = [(e[2] - e[1])/2 for e in ex]
    @test radii[3] ≈ zero(radii[3])
    return (radii[1] + radii[2])/2
end

function constant_velocity_args()
    xinit = @SVector zeros(3)
    vinit = @SVector [1.0, 0.0, 0.0]
    dt = 0.001
    B = @SVector [0.0, 1.0, 0.0]
    E = -vinit × B # E = [0, 0, -1]
    N = 10000
    args = (xinit,vinit,dt,(x,y,z)->E,(x,y,z)->B,N)
    kwargs = Dict(:m=>1,:q=>1)
    return args,kwargs
end

function E_zero_args()
    xinit = @SVector zeros(3)
    vinit = @SVector [1.0, 0.0, 0.0]
    dt = 0.001
    B = @SVector [0.0, 0.0, 1.0]
    E = @SVector zeros(3)
    N = 10000
    args = (xinit,vinit,dt,(x,y,z)->E,(x,y,z)->B,N)
    kwargs = Dict(:m=>1,:q=>1)
    return args,kwargs
end

function add_units(gen_func)
    (xinit,vinit,dt,E,B,N), kw = gen_func()
    xinit = xinit*u"m"
    vinit = vinit*u"m/s"
    dt = dt*u"s"
    Ẽ = (x,y,z)->E(x,y,z)*u"V/m"
    B̃ = (x,y,z)->B(x,y,z)*u"T"
    kw = Dict(:m=>kw[:m]u"kg", :q=>kw[:q]u"C")
    return (xinit,vinit,dt,Ẽ,B̃,N), kw
end
function real_electron_args()
    xinit = @SVector(zeros(3))u"m"
    vinit = SA[0.0, 1e5, 0.0]u"m/s"
    dt = 3e-11u"s"
    E = @SVector(zeros(3))u"V/m"
    B = @SVector(zeros(3))u"T"
    N = 1000
    args = (xinit,vinit,dt,(x,y,z)->E,(x,y,z)->B,N)
    kw = Dict(:m=>m_e, :q=>-e)
    return args, kw
end

function test_const_velocity(x,v)
    Δx = diff(x)
    @test all(v .≈ Ref(v[1]))
    @test all(Δx .≈ Ref(Δx[1]))
end
function test_pure_gyromotion(x,v)
    vmag = norm.(v)
    gr = _estimate_gyroradius(x)
    @test gr ≈ oneunit(gr) atol=0.0001*oneunit(gr)
    @test all(vmag .≈ Ref(vmag[1]))
end

@testset "All Tests" begin
    @testset "Constant Velocity solution" begin
        args, kwargs = constant_velocity_args()
        x,v = boris(args...; kwargs...)
        test_const_velocity(x,v)
        @testset "with units" begin
            args, kwargs = add_units(constant_velocity_args)
            x,v = boris(args...; kwargs...)
            test_const_velocity(x,v)
        end
    end
    @testset "Pure Gyromotion, E=0" begin
        args, kwargs = E_zero_args()
        x,v = boris(args...; kwargs...)
        test_pure_gyromotion(x,v)
        @testset "with units" begin
            args, kwargs = add_units(E_zero_args)
            x,v = boris(args...; kwargs...)
            test_pure_gyromotion(x,v)
        end
    end
end
