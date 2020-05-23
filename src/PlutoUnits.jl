module PlutoUnits
using Unitful
@unit Rp "Rp" PlutoRadius 1189.0u"km" false

# localunits must be defined after all new units are made
const localunits = Unitful.basefactors
## Only if @dimension is used
#const localpromotion = Unitful.promotion
function __init__()
    merge!(Unitful.basefactors, localunits)
    
    ## Only if @dimension is used
    # merge!(Unitful.promotion, localpromotion)

    Unitful.register(PlutoUnits)
end
end # module
