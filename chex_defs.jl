import PhysicalConstants.CODATA2018: G, k_B, m_p, e
using SpecialFunctions

P(α,x) = gamma_inc(α,x)[1]
Q(α,x) = gamma_inc(α,x)[2]
γ(α,x) = gamma_inc(α,x)[1]*gamma(α)

const M = 1.309e22u"kg"
const m = 16m_p
const r_c = 2900u"km"
const T_c = 66u"K"
const N_c = 3e6u"cm^-3"
λ(r) = ustrip(Unitful.NoUnits, G*M*m/(k_B*T_c*r))
const λ_c = λ(r_c)
ψ₁(λ) = λ^2/(λ + λ_c)

ζ_bal_sat(λ) = P(3/2, λ)
ζ_esc(λ) = 1/2 * (Q(3/2,λ) - sqrt(λ_c^2 - λ^2)/λ_c * exp(-ψ₁(λ))*Q(3/2,λ-ψ₁(λ)))

ζ_bal(λ) = (
    P(3/2,λ) - sqrt(λ_c^2 - λ^2)/λ_c * exp(-ψ₁(λ)) * P(3/2, λ - ψ₁(λ))
)

ζ_bal_sat_esc(λ) = ζ_bal_sat(λ) + ζ_esc(λ)
ζ_bal_esc(λ) = ζ_bal(λ) + ζ_esc(λ)

N_b(r) = N_c*exp(-(λ_c - λ(r)))
N_bal_sat_esc(r) = N_b(r)*ζ_bal_sat_esc(λ(r))
N_bal_esc(r) = N_b(r)*ζ_bal_esc(λ(r))

N_power_law(r) = 1e15*(1e15*(1188u"km"/r)^25 + 5e9*(1188u"km"/r)^8)u"km^-3"

σ_standard(v) = 1e-15u"cm^2"
σ_high(v) = 1e-14u"cm^2"
