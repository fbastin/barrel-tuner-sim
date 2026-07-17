# =============================================================================
# AMPLITUDE DU TAUX ANGULAIRE DE BOUCHE — le fusil libre atteint-il l'ordre de Kolbe ?
#
# Contexte. Le modèle à culasse ENCASTRÉE (harral_a22lr.jl, simulation.jl) donne
# le bon SIGNE de θ̇(t_b) mais une amplitude ~3× trop faible pour les 6 MOA/ms
# mesurés par Kolbe. La roadmap l'attribuait à « la rotation de corps rigide
# manquante ». Ce script teste cette hypothèse quantitativement, sur le modèle
# FUSIL LIBRE de harral_rifle_sweep.jl (canon + boîte + crosse + avant-bras sur
# deux sacs), avec :
#   - excitation par le moment de recul PHYSIQUE (∫A_bore·p dt = (1+β)·m_p·v) ;
#   - amortissement d'acier léger calé sur les modes de WHIP (pas les sacs).
# Les deux corrections de 2026-07-17 qui ont assaini `harral_rifle_sweep.jl`.
#
# Ce qu'on lit : θ̇(L, t_b) et le PIC de |θ̇|, canon nu vs tuner, en MOA/ms.
# Référence Kolbe (canon 26", 50 m) : −9,4 MOA/ms nu, +6,0 tuned.
#
# RÉSULTAT (voir BILAN en fin) : l'amplitude (le PIC) est de l'ordre de Kolbe,
# et h_bore la commande. Mais la rotation de corps rigide seule n'y contribue
# que ~0,2 MOA/ms : c'est le WHIP élastique qu'elle EXCITE qui porte l'amplitude.
# En revanche θ̇ À t_b (ce qui fixe le lancement) reste petit et de structure
# opposée à Kolbe : la sortie tombe près d'un rebroussement du whip, et caler
# cet instant exige la géométrie/le t_b réels, non publiés.
#
# Usage :  julia kolbe_amplitude.jl
# =============================================================================

include(joinpath(@__DIR__, "harral_rifle_sweep.jl"))
using Printf

const V0 = 320.0                 # m/s (.22 LR ~1050 fps)
const MOAms = (180*60/π) * 1e-3  # rad/s → MOA/ms

# Config d'arme représentative (milieu des plages plausibles du balayage)
const CFG = (x_breech = 0.5, L_fore = 0.25, x_rear = 0.038,
             EI_stock = 1e4, k_rest = 1e5)

# -----------------------------------------------------------------------------
# (0) Estimation analytique : taux angulaire de CORPS RIGIDE dû au recul
#   ω = J·h / I,  J = (1+β)·m_p·v (impulsion),  I = inertie /CG.
# -----------------------------------------------------------------------------
function rigid_body_rate(; h_bore = 0.0254, m_rifle = M_RIFLE_TOTAL, L_rifle = 1.1)
    J = (1 + β_GAZ) * m_p * V0
    I = m_rifle * L_rifle^2 / 12
    return J * h_bore / I * MOAms          # MOA/ms
end

# -----------------------------------------------------------------------------
# (1) Fusil libre : θ̇(t_b) et pic de |θ̇|, via shoot_rifle (excitation physique
#     + amortissement de whip déjà intégrés).
# -----------------------------------------------------------------------------
shoot(m_tuner; h_bore = 0.0254) =
    shoot_rifle(V0; m_tuner = m_tuner, h_bore = h_bore, CFG...)

function main()
    println("="^76)
    println(" Amplitude du taux angulaire de bouche — fusil libre vs Kolbe")
    println("="^76)

    @printf("\n(0) Corps rigide SEUL (analytique), h_bore = 1\" : ω = %.2f MOA/ms\n",
            rigid_body_rate())
    println("    ⇒ négligeable devant les 6-9 de Kolbe : la rotation de corps rigide")
    println("      n'est pas une source directe d'amplitude, elle EXCITE le whip.\n")

    println("(1) Fusil libre sur deux sacs, recul physique, amortissement d'acier :")
    @printf("    %-8s | %-10s | %-15s | %s\n",
            "tuner", "f_whip Hz", "θ̇(t_b) MOA/ms", "pic |θ̇| MOA/ms")
    println("    " * "-"^58)
    for (lbl, m) in (("nu", 0.0), ("100 g", 0.1), ("200 g", 0.2), ("300 g", 0.3))
        r = shoot(m)
        if r === nothing
            @printf("    %-8s | (config instable — l'arme bascule)\n", lbl)
        else
            @printf("    %-8s | %-10.0f | %+-15.2f | %.1f\n",
                    lbl, r.f1, r.θdot_tb, r.θdot_peak)
        end
    end
    println("\n    Kolbe (mesuré) : nu −9,4 / 200 g +6,0 MOA/ms.")
    println("    Clamped (harral_a22lr) : nu −1,1 / 200 g +1,8 MOA/ms (pic ~10, calibré).")

    println("\n(2) Sensibilité du PIC de |θ̇| à la hauteur d'âme h_bore (canon nu) :")
    for hb in (0.0254, 0.038, 0.051)
        r = shoot(0.0; h_bore = hb)
        r === nothing || @printf("    h_bore = %.0f mm : pic |θ̇| = %.1f MOA/ms\n",
                                 hb*1e3, r.θdot_peak)
    end

    println("\n" * "-"^76)
    println("BILAN")
    println("  • Le PIC de |θ̇| atteint l'ORDRE de Kolbe (~5 MOA/ms à h_bore=1\",")
    println("    ~10 à 2\"), à partir du seul recul physique : l'amplitude n'est plus")
    println("    un mystère ni un paramètre de calage. ✅ hypothèse d'amplitude validée.")
    println("  • Mécanisme précisé : la rotation de corps rigide (~0,2 MOA/ms) ne")
    println("    fournit pas l'amplitude ; elle EXCITE le whip élastique, qui l'amplifie.")
    println("  • MAIS θ̇ À t_b reste petit et de structure opposée à Kolbe (sortie près")
    println("    d'un rebroussement du whip). Caler cet instant exige la fréquence de")
    println("    whip ET le t_b réels — donc la géométrie de l'arme, non publiée.")
    println("    ⇒ l'amplitude est récupérable ; le calage EXACT sur Kolbe ne l'est pas.")
    println("="^76)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
