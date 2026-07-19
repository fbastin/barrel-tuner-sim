# =============================================================================
# PCP CONTRE ARME À FEU : L'EXCITATION EST-ELLE VRAIMENT NÉGLIGEABLE ?
#
# La documentation du site affirmait qu'un PCP est « quasi sans recul » et que
# le modèle y « devient quasi muet ». Objection reçue le 2026-07-19 : la
# conservation de la quantité de mouvement s'applique aussi à un PCP, qui
# recule donc bel et bien. Ce script mesure l'écart réel au lieu de le supposer.
#
# CE QUI EST EN JEU. Deux rapports très différents circulent :
#   • poussée CRÊTE à la culasse  : ~26 : 1   (200 MPa × 24,6 mm² contre
#                                              ~12 MPa × 15,9 mm²)
#   • quantité de MOUVEMENT       : ~3 à 6 : 1
# Ils diffèrent parce que la détente d'un PCP est plate et longue là où la
# combustion est brève et pointue. Lequel gouverne la vibration dépend du
# RÉGIME : le mode fondamental est à ~40 Hz (période 25 ms) pour un chargement
# de ~3 ms, donc en régime impulsionnel la réponse suit l'IMPULSION, pas le pic.
# Si c'est le cas, « quasi muet » est faux d'environ un facteur 5 à 8.
#
# Usage :   julia pcp_vs_firearm.jl
# =============================================================================

include(joinpath(@__DIR__, "simulation.jl"))

using Printf

const H = H_OFFSET_EFF

# Intégration trapézoïdale sur grille fine — les autres scripts du dépôt sont
# sans dépendance externe, on ne tire pas QuadGK pour une intégrale lisse 1-D.
# Pas de 1 µs sur 30 ms : l'erreur est très inférieure aux incertitudes physiques.
function integrate(f, a, b; n = 30_000)
    h = (b - a) / n
    s = 0.5 * (f(a) + f(b))
    for i in 1:n-1
        s += f(a + i * h)
    end
    return s * h
end

# -----------------------------------------------------------------------------
# 1. BUDGETS DE QUANTITÉ DE MOUVEMENT (indépendants du modèle vibratoire)
# -----------------------------------------------------------------------------
# Arme à feu .22 LR : balle + gaz de poudre (~0,1 g éjectés à ~1,5·v₀).
const P_FIREARM = 2.6e-3 * 318.0 + 0.10e-3 * 1.5 * 318.0

# PCP .177 haute puissance : plomb + air chassé. L'air compte pour la moitié du
# bilan, ce qui n'est pas intuitif — le plomb seul sous-estime le recul d'un
# facteur 2. Densité de l'air prise à ~36 kg/m³ (résiduel ~30 bar dans le canon).
const M_PELLET, V_PELLET = 0.55e-3, 260.0
const A_BORE_177 = π/4 * (4.5e-3)^2
const M_AIR      = 0.60 * A_BORE_177 * 36.0
const P_PCP_FULL = M_PELLET * V_PELLET + M_AIR * 1.5 * V_PELLET
# Part délivrée AVANT la sortie du plomb : l'air continue de pousser après, et
# cette fraction-là n'excite plus rien d'utile pour la compensation.
const P_PCP_EXIT = M_PELLET * V_PELLET

# -----------------------------------------------------------------------------
# 2. PROFILS DE MOMENT À LA CULASSE
# -----------------------------------------------------------------------------
# Forme « arme à feu » : l'ANCIEN profil autonome de simulation.jl (gamma, pic
# 200 MPa à 0,5 ms). Il est conservé ICI, et nulle part ailleurs, parce que ce
# script est celui qui documente son défaut : son impulsion valait 5,19× le
# recul physique. simulation.jl ne le porte plus depuis le 2026-07-19, sa
# pression étant désormais dérivée de la balistique intérieure couplée.
const LEGACY_PMAX, LEGACY_TPEAK, LEGACY_ALPHA = 200e6, 0.5e-3, 2.0
legacy_pressure(t) = (t <= 0 || t >= 5 * LEGACY_TPEAK) ? 0.0 :
    (u = t / LEGACY_TPEAK; LEGACY_PMAX * u^LEGACY_ALPHA * exp(LEGACY_ALPHA * (1 - u)))
shape_firearm(t) = legacy_pressure(t) / LEGACY_PMAX

# Forme « PCP » : montée rapide à l'ouverture de soupape puis décroissance
# lente — plate et longue, à l'opposé du pic de combustion.
const τ_RISE, τ_DECAY = 0.15e-3, 2.0e-3
shape_pcp(t) = t <= 0 ? 0.0 : (1 - exp(-t / τ_RISE)) * exp(-t / τ_DECAY)

# Met une forme à l'échelle pour que ∫M dt = h · p_cible (impulsion angulaire).
function scaled_moment(shape, p_target; h = H)
    I = integrate(shape, 0.0, 30e-3)
    k = h * p_target / I
    return t -> k * shape(t)
end

# -----------------------------------------------------------------------------
# 3. MESURE
# -----------------------------------------------------------------------------
function run_case(label, moment, v_p, m_proj; m_tuner = 0.0, d = 0.0)
    r = simulate_shot(m_tuner; d_overhang = d, h_offset = H,
                      moment_of_t = moment, v_p = v_p, m_proj = m_proj,
                      verbose = false)
    pic = maximum(abs.(r.θdot_L)) |> moa_per_ms
    return (; label, θ = r.θ_at_tb, θdot = r.θdot_MOAms, pic, t_b = r.t_b)
end

println("="^78)
println(" PCP contre arme à feu — l'excitation est-elle négligeable ?")
println("="^78)
@printf("\nBudgets de quantité de mouvement (conservation, hors modèle) :\n")
@printf("  arme à feu .22 LR      p = %.3f kg·m/s\n", P_FIREARM)
@printf("  PCP .177 (total)       p = %.3f kg·m/s   → rapport %.1f : 1\n",
        P_PCP_FULL, P_FIREARM / P_PCP_FULL)
@printf("  PCP .177 (avant sortie) p = %.3f kg·m/s   → rapport %.1f : 1\n",
        P_PCP_EXIT, P_FIREARM / P_PCP_EXIT)

# Contrôle : le profil de pression du modèle est-il cohérent avec le budget ?
I_model = integrate(t -> legacy_pressure(t) * A_bore, 0.0, 30e-3)
@printf("\nAncien profil autonome (retiré du modèle) : ∫p·A dt = %.3f kg·m/s", I_model)
@printf("  (budget balistique %.3f)\n", P_FIREARM)
@printf("  → il portait %.2f× le recul physique. Le modèle actuel dérive p(t) de la\n", I_model / P_FIREARM)
@printf("    balistique intérieure : la conservation y est exacte par construction.\n")

m_fa  = scaled_moment(shape_firearm, P_FIREARM)
m_pcp = scaled_moment(shape_pcp,     P_PCP_FULL)
m_pcp_exit = scaled_moment(shape_pcp, P_PCP_EXIT)
# Le contrefactuel : forme d'arme à feu réduite au rapport des PICS (26:1),
# soit l'hypothèse implicite du « quasi muet ».
m_peak_scaled = let f = scaled_moment(shape_firearm, P_FIREARM)
    t -> f(t) / 25.8
end

cases = [
    run_case("arme à feu .22 LR",              m_fa,          318.0, 2.6e-3),
    run_case("PCP (impulsion totale)",         m_pcp,         V_PELLET, M_PELLET),
    run_case("PCP (impulsion avant sortie)",   m_pcp_exit,    V_PELLET, M_PELLET),
    run_case("contrefactuel : pic ÷ 25,8",     m_peak_scaled, 318.0, 2.6e-3),
]

println("\n" * "-"^78)
@printf("%-32s | %9s | %12s | %10s\n", "cas (canon nu)", "t_b (ms)", "θ̇(t_b)", "pic |θ̇|")
println("-"^78)
for c in cases
    @printf("%-32s | %9.2f | %+9.3f    | %9.2f\n",
            c.label, c.t_b * 1e3, c.θdot, c.pic)
end
println("-"^78)

fa, pcp, pcp_x, cf = cases
@printf("\nRapport des PICS de |θ̇| (ce que le modèle répond réellement) :\n")
@printf("  arme à feu / PCP total        = %5.1f : 1\n", fa.pic / pcp.pic)
@printf("  arme à feu / PCP avant sortie = %5.1f : 1\n", fa.pic / pcp_x.pic)
@printf("  arme à feu / contrefactuel    = %5.1f : 1   (par construction 25,8)\n",
        fa.pic / cf.pic)

println("\n" * "="^78)
println("LECTURE")
println("="^78)
r_imp = fa.pic / pcp.pic
if r_imp < 10
    println("Le modèle suit l'IMPULSION, non la poussée crête : l'écart mesuré est de")
    @printf("l'ordre de %.0f : 1, non de 26 : 1. L'affirmation « quasi muet » surestime\n", r_imp)
    println("donc l'écart d'environ un facteur ", round(25.8 / r_imp, digits = 1), ".")
    println()
    println("RAISON. Le chargement (~3 ms) est court devant la période du mode")
    println("fondamental (~25 ms à 40 Hz) : on est en régime impulsionnel, où la")
    println("réponse modale est fixée par ∫M dt et non par max M. La forme du profil")
    println("n'intervient qu'au second ordre — c'est bien ce que montre la ligne")
    println("« contrefactuel », qui garde la forme d'arme à feu et ne diffère que par")
    println("l'échelle.")
else
    println("Le modèle suit la poussée crête davantage que l'impulsion : l'écart")
    @printf("mesuré (%.0f : 1) reste du même ordre que le rapport des pics.\n", r_imp)
end
println()
println("CE QUE CELA NE DIT PAS. Qu'un PCP mérite un modèle distinct reste vrai,")
println("mais pour une autre raison : il possède une excitation que l'arme à feu")
println("n'a pas — la frappe marteau/soupape — dont le calage temporel diffère")
println("(cf. pcp_tuner.jl). L'argument central de ce dernier, à savoir que la")
println("POSITION des sweet spots ne dépend pas de l'amplitude d'excitation, est")
println("indépendant de tout ce qui précède.")
println("="^78)
