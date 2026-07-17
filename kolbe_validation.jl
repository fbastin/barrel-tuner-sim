# =============================================================================
# VALIDATION contre les mesures de G. Kolbe (Border Barrels)
#   « Using barrel vibrations to tune a barrel »
#   http://www.geoffrey-kolbe.com/articles/rimfire_accuracy/tuning_a_barrel.htm
#
# POURQUOI KOLBE PLUTÔT QUE HARRAL. Harral ne publie pas les données de son
# arme (un seul chiffre : 10.5 lb) et — Kolbe le dit — « his work has lacked
# the experimental confirmation needed to verify his computer modelling ».
# Nous avons donc passé le début de ce chantier à courir après un modèle que
# personne n'a validé. Kolbe, LUI, A MESURÉ, et il publie une chaîne de
# nombres dont les premiers maillons sont de la cinématique pure —
# REPRODUCTIBLES SANS AUCUNE INCONNUE. C'est la validation qui nous manquait.
#
# La chaîne de Kolbe pour la compensation positive COMPLÈTE à 50 m :
#   (a) chute naturelle              : 0.016 MOA par ft/s de vitesse
#   (b) sensibilité du temps de sortie : 375 ft/s par ms  (canon 26")
#   (c) taux requis = (a) × (b)      : 6.0 MOA/ms
# et ses mesures sur banc :
#   canon nu  : −9.4 MOA/ms (bouche DESCENDANTE) → dispersion verticale
#   + 200 g   : +6.0 MOA/ms (bouche montante)   → groupements ronds
#
# ⚠️ Le MÉCANISME, chez Kolbe : « The transverse vibrations are due to the
# recoil forces in the rifle imparting a moment on the back of the barrel AS
# THE RIFLE ROTATES ABOUT ITS CENTRE OF GRAVITY. » L'amplitude dépend donc de
# deux grandeurs MESURABLES sur l'arme : son poids et la distance CG↔âme.
# Corollaire pratique : « The model assumes that the rifle is essentially in
# free space and unconstrained by any rests or clamps as it recoils. Shooting
# a rifle off bags is a fair approximation to this. However, if small calibre
# rifle is gripped tightly or pulled hard into the shoulder then the recoil
# dynamics could be affected. »
#
# Usage :  julia kolbe_validation.jl
# =============================================================================

using Printf

const IN    = 0.0254
const FPS   = 0.3048              # 1 ft/s → m/s
const G     = 9.80665             # m/s²
const D50M  = 50.0                # portée de référence de Kolbe : 50 MÈTRES
const MOA   = (1/60) * π/180      # 1 MOA en radians

moa_at(d) = d * MOA               # taille d'1 MOA à la distance d (m)

# -----------------------------------------------------------------------------
# (a) CHUTE NATURELLE : 0.016 MOA par ft/s
#   Argument de Kolbe : « there will naturally be a drop at the target of 0.016
#   MOA for every ft/sec. drop in muzzle velocity, JUST BECAUSE THE BULLET
#   TAKES A LITTLE LONGER TO GET THERE. » C'est de la cinématique pure, sans
#   traînée : t = d/v, chute = ½gt² ⇒ d(chute)/dv = −g·d²/v³.
# -----------------------------------------------------------------------------
drop_sensitivity_moa_per_fps(v, d) = (G * d^2 / v^3) * FPS / moa_at(d)

# -----------------------------------------------------------------------------
# (b) SENSIBILITÉ DU TEMPS DE SORTIE : dv/dt_b
#   Notre profil de balistique intérieure (celui de harral_a22lr.jl) :
#     x(t) = v_m·(t − τ(1 − e^{−t/τ}))   ⇒  t_b tel que x = L
# -----------------------------------------------------------------------------
function exit_time(v_m, L; τ = 0.40e-3)
    x(t) = v_m * (t - τ*(1 - exp(-t/τ)))
    lo, hi = 1e-6, 20e-3
    for _ in 1:200
        mid = 0.5*(lo+hi)
        x(mid) < L ? (lo = mid) : (hi = mid)
    end
    return 0.5*(lo+hi)
end

function dv_dtb_fps_per_ms(v_m, L; dv = 0.5)
    t1 = exit_time(v_m - dv, L); t2 = exit_time(v_m + dv, L)
    return (2dv / FPS) / ((t1 - t2) * 1e3)      # ft/s par ms (positif)
end

# =============================================================================
function main()
println("="^78)
println(" Validation contre les mesures publiées de G. Kolbe (50 m, .22 LR)")
println("="^78)

# --- (a) ---------------------------------------------------------------------
println("\n(a) CHUTE NATURELLE — cinématique pure, AUCUNE inconnue")
println("    Kolbe : 0.016 MOA par ft/s à 50 m")
@printf("    %-14s | %-12s | %s\n", "v (ft/s)", "v (m/s)", "MOA par ft/s")
println("    " * "-"^48)
for vfps in (1000.0, 1035.0, 1050.0, 1085.0, 1100.0)
    v = vfps * FPS
    @printf("    %-14.0f | %-12.1f | %.4f\n", vfps, v, drop_sensitivity_moa_per_fps(v, D50M))
end
# Vitesse qui redonne EXACTEMENT le 0.016 de Kolbe
let lo = 200.0, hi = 500.0
    for _ in 1:200
        mid = 0.5*(lo+hi)
        drop_sensitivity_moa_per_fps(mid, D50M) > 0.016 ? (lo = mid) : (hi = mid)
    end
    v = 0.5*(lo+hi)
    @printf("    → 0.0160 MOA/ft/s correspond à v = %.1f m/s (%.0f ft/s) : plausible pour de l'Eley.\n",
            v, v/FPS)
end
println("    ✅ REPRODUIT : le chiffre de Kolbe tombe dans la plage .22 LR réelle.")

# --- (b) ---------------------------------------------------------------------
println("\n(b) SENSIBILITÉ DU TEMPS DE SORTIE — dépend de NOTRE balistique intérieure")
println("    Kolbe : 375 ft/s par ms (canon 26\")")
L26 = 26.0 * IN
@printf("    %-14s | %-12s | %s\n", "v (ft/s)", "t_b (ms)", "dv/dt_b (ft/s par ms)")
println("    " * "-"^54)
for vfps in (1035.0, 1050.0, 1085.0)
    v = vfps * FPS
    @printf("    %-14.0f | %-12.3f | %.0f\n", vfps, exit_time(v, L26)*1e3,
            dv_dtb_fps_per_ms(v, L26))
end
v50 = 1050.0 * FPS
ours = dv_dtb_fps_per_ms(v50, L26)
@printf("    ⚠️ ÉCART : nous %.0f contre 375 chez Kolbe (facteur %.2f).\n", ours, ours/375)
println("    Notre profil donne asymptotiquement dv/dt_b = v²/L, indépendant de τ.")
@printf("    Le 375 de Kolbe implique t_b = v/375 ≈ %.2f ms ; nous trouvons %.2f ms.\n",
        (v50/FPS)/375, exit_time(v50, L26)*1e3)
println("    ⇒ Notre balistique intérieure sort la balle TROP TÔT. Ce n'est pas un")
println("      détail : (c) est le PRODUIT de (a) et (b), donc l'écart s'y propage.")

# --- (c) ---------------------------------------------------------------------
println("\n(c) TAUX REQUIS POUR COMPENSATION COMPLÈTE À 50 m = (a) × (b)")
println("    Kolbe : 0.016 × 375 = 6.0 MOA/ms")
a = drop_sensitivity_moa_per_fps(v50, D50M)
@printf("    Nous  : %.4f × %.0f = %.1f MOA/ms\n", a, ours, a*ours)
@printf("    Avec le (b) de Kolbe : %.4f × 375 = %.1f MOA/ms  ← l'écart vient bien de (b)\n",
        a, a*375)

println("\n" * "-"^78)
println("BILAN")
println("  (a) chute naturelle           : ✅ reproduit exactement (cinématique pure)")
println("  (b) sensibilité temps de sortie : ❌ écart ~35 % — NOTRE balistique intérieure")
println("  (c) taux requis                : hérite de l'écart de (b)")
println()
println("  Le maillon faible est identifié et il est INTERNE : le profil v(t) de")
println("  `harral_a22lr.jl` (τ = 0.40 ms, calé sur rien de publié) fait sortir la")
println("  balle trop tôt. C'est réparable — et contrairement aux inconnues de")
println("  Harral, ça ne dépend que de nous.")
println()
println("  Repères de Kolbe pour la suite (canon 26\", 50 m) :")
println("    canon nu  : −9.4 MOA/ms (bouche descendante) → dispersion verticale")
println("    + 200 g   : +6.0 MOA/ms → groupements ronds, compensation complète")
println("    ⚠️ Ces deux-là dépendent de la souplesse de son banc, NON publiée :")
println("       hors de portée, comme le tableau de Harral. Ne pas s'y caler.")
println("="^78)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
