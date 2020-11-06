using NearestNeighbors
using Unitful

trystrip(unit, value::Quantity) = ustrip(unit, value)
trystrip(unit, values::Vector{<:Quantity}) = ustrip.(unit, values)
trystrip(unit, value) = value

touching(tree::NNTree, x, dx) = inrange(tree, trystrip.(u"km", x), trystrip(u"km", dx))

# Fallback is to assume it acts like a StructArray
buildtree(ps) = KDTree([trystrip.(u"km", p.x) for p in ps], Chebyshev())
buildtree(ps::Vector) = KDTree(trystrip.(ps.x), Chebyshev())
