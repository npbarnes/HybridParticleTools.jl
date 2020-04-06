module ParameterSets
export ParameterSet

using PyCall

struct ParameterSet
    nx::Int32
    ny::Int32
    nz::Int32
    dx::Float32
    dy::Float32
    delz::Float32
    nt::Int32
    dtsub_init::Float32
    ntsub::Int32
    dt::Float32
    nout::Int32
    out_dir::String
    vtop::Float32
    vbottom::Float32
    Ni_max::Int32
    mproton::Float32
    m_pu::Float32
    m_heavy::Float32
    np_top::Float32
    np_bottom::Float32
    b0_top::Float32
    b0_bottom::Float32
    vth_top::Float32
    vth_bottom::Float32
    alpha::Float64
    beta::Float32
    RPluto::Float32
    RIo::Float32
    b0_init::Float32
    ion_amu::Int32
    mpu::Float32
    nf_init::Float32
    dt_frac::Float32
    vsw::Float32
    vth::Float32
    Ni_tot_frac::Float32
    dx_frac::Float32
    nu_init_frac::Float32
    mrestart::Int32
    pluto_offset::Int64
    ri0::Int32
    part_nout::Int32
    qx::Array{Float32,1}
    qy::Array{Float32,1}
    qz::Array{Float32,1}
    dz_grid::Array{Float32,1}
    dz_cell::Array{Float32,1}
    num_proc::Int64
    zrange::Int64
    saved_steps::Float64
    simulation_height::Float64
    qzrange::Array{Float64,1}
    grid_points::Tuple{Array{Float32,1},Array{Float32,1},Array{Float64,1}}
end
ParameterSet(o::PyObject) = ParameterSet(convert(Dict{Symbol,Any}, o))
function ParameterSet(d::Dict{Symbol,Any})
    ParameterSet(
        d[:nx],
        d[:ny],
        d[:nz],
        d[:dx],
        d[:dy],
        d[:delz],
        d[:nt],
        d[:dtsub_init],
        d[:ntsub],
        d[:dt],
        d[:nout],
        d[:out_dir],
        d[:vtop],
        d[:vbottom],
        d[:Ni_max],
        d[:mproton],
        d[:m_pu],
        d[:m_heavy],
        d[:np_top],
        d[:np_bottom],
        d[:b0_top],
        d[:b0_bottom],
        d[:vth_top],
        d[:vth_bottom],
        d[:alpha],
        d[:beta],
        d[:RPluto],
        d[:RIo],
        d[:b0_init],
        d[:ion_amu],
        d[:mpu],
        d[:nf_init],
        d[:dt_frac],
        d[:vsw],
        d[:vth],
        d[:Ni_tot_frac],
        d[:dx_frac],
        d[:nu_init_frac],
        d[:mrestart],
        d[:pluto_offset],
        d[:ri0],
        d[:part_nout],
        d[:qx],
        d[:qy],
        d[:qz],
        d[:dz_grid],
        d[:dz_cell],
        d[:num_proc],
        d[:zrange],
        d[:saved_steps],
        d[:simulation_height],
        d[:qzrange],
        d[:grid_points]
    )
end
Base.show(io::IO, hp::ParameterSet) = print(io, "ParameterSet(...)")
function Base.show(io::IO, ::MIME"text/plain", hp::ParameterSet)
    print(io, "HybridParams:")
    for f in fieldnames(typeof(hp))
        print(io, "\n\t$f: $(getfield(hp,f))")
    end
end

function ParameterSet(path)
    pymod = pyimport("HybridParams")
    pyhp = pymod.HybridParams(path)
    return ParameterSet(pyhp."para")
end

end #module
