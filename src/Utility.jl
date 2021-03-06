module Utility
export  listofvectors, listofmatrices, geomspace, asunitless, viewasarray,
        @utc2et_str, @utc2datetime_str, getcolumns, r_transform, utc2et, rotmat, rotmatX
using StaticArrays
using StructArrays
using Unitful
using PyCall

const st = PyNULL()

function __init__()
    copy!(st, pyimport("spice_tools"))
end

"""
    listofvectors(array)
Convert Array{T,2} to Vector{SVector{3,T}}.
If size(arr,1) != 3, there will be an error.
A copy is made.
"""
function listofvectors(arr::AbstractArray{T,2}) where T
    SVector{3,T}.(eachrow(arr))
end
"Convert Array{T,3} with size (3,3,:) to Vector{SMatrix{3,3,T}}"
function listofmatrices(arr::AbstractArray{T,3}) where T
    SMatrix{3,3,T,9}.(eachslice(arr,dims=1))
end

function geomspace(start, stop, N)
    exp10.(range(log10(start), stop=log10(stop), length=N))
end

"No-copy ustrip for Arrays of StaticArrays of Quantities"
function asunitless(u::Unitful.Units, A::AbstractArray{SA}) where {Size, Q<:Quantity, SA<:StaticArray{Size,Q}}
    ret = similar(A, similar_type(SA, Unitful.numtype(Q)))
    for (i,a) in enumerate(A)
        ret[i] = ustrip.(u, a)
    end
    return ret
end

asunitless(A::AbstractArray{SA}) where {Size, Q<:Quantity, SA<:StaticArray{Size, Q}} = reinterpret(similar_type(SA,Unitful.numtype(Q)), A)


viewasarray(x::AbstractArray{SA,1}) where {S,T,SA<:StaticArray{S,T}} = reshape(reinterpret(T,x), (size(SA)...,:))
getcolumns(x::AbstractArray{SA<:StaticArray,1}) = Tuple(getindex.(x, i) for i in 1:length(SA))

utc2et(tstr) = st.sp.str2et("2015-7-14T$tstr")
macro utc2et_str(t)
    :(utc2et($t))
end
macro utc2datetime_str(t)
    :(st.et2pydatetime(@utc2et_str $t))
end

# This specific method is needed for filter! to work on StructArrays
function deleteat!(sa::StructArray, r::UnitRange{<:Integer})
    for field in fieldarrays(sa)
        deleteat!(field, r)
    end
end

"""2d rotation matrix"""
rotmat(t) = SA[
    cos(t) -sin(t);
    sin(t) cos(t)
]

"""3d rotation matrix round the X axis"""
rotmatX(t) = SA[
    1    0     0;
    0 cos(t) -sin(t);
    0 sin(t) cos(t);
]

"""r_transform is the IMF rotation transformation for vectors"""
r_transform(v, θ=180.0) = rotmatX(θ) * v


end # module
