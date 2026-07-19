# =============================================================================
# CADRE v2 — les trois pièces validées, réunies pour la première fois
#
# L'enquête de juillet 2026 a invalidé successivement chaque pièce du cadre
# initial, sans jamais rapprocher des mesures de Kolbe. Le diagnostic final
# n'était pas « telle constante est fausse » mais « le cadre a atteint sa
# limite ». Ce script le rebâtit — non pas de zéro, mais en COMPOSANT les trois
# pièces que l'enquête a validées séparément et qui n'avaient jamais tourné
# ensemble :
#
#   1. STRUCTURE  — fusil entier sur deux sacs à contact unilatéral (crosse,
#      canon, avant-bras). Seule configuration atteignant l'ordre de grandeur
#      de Kolbe : pic |θ̇| ≈ 7 contre 9,4 mesurés (free_boundary.jl a montré
#      qu'un canon libre mais sans crosse plafonne à 0,5).
#   2. EXCITATION — balistique intérieure couplée : la cinématique est intégrée
#      depuis la pression, donc ∫p·A dt = m_eff·v exact par construction. Le
#      profil autonome antérieur délivrait 5,19× le recul physique.
#   3. TUNER      — tube ÉLASTIQUE maillé en poutre, masse coulissante à
#      l'intérieur, en aluminium. La masse ponctuelle rigide gonflait l'effet
#      d'un ordre de grandeur, et l'acier ne laissait aucune masse au curseur.
#
# CE QUI N'EST TOUJOURS PAS DANS LE CADRE, et qu'il faut avoir en tête :
#   • le mouvement est PLANAIRE (vertical seul) ; la bouche décrit en réalité
#     une orbite 2D, dont la composante horizontale n'est compensée par rien ;
#   • l'amortissement à 2-3 kHz n'est pas mesuré, et c'est le paramètre dont
#     dépend le plus la conclusion (cf. le balayage de ζ) ;
#   • l'excitation reste le seul moment de recul : ni gravure, ni frottement
#     du projectile, ni frappe de percuteur.
#
# Usage :   julia frame_v2.jl
# =============================================================================

include(joinpath(@__DIR__, "rifle_coupled.jl"))
using Printf

# Arme de référence (mêmes inconnues que kolbe_amplitude.jl et rifle_coupled.jl)
const CFG_V2 = (x_breech = 0.35, L_fore = 0.15, x_rear = 0.038,
                EI_stock = 1e4, k_rest = 1e5)

# Budget de masse d'un tuner à tube ALUMINIUM : le tube consomme ρA·L, le
# curseur reçoit le reste. C'est ce budget qui rend le réglage en position
# physiquement possible — en acier, il ne resterait rien à faire coulisser.
function tube_budget(m_total, L_tube)
    m_tube = TUBE_RHOA_R * L_tube
    return (m_tube = m_tube, m_slider = max(m_total - m_tube, 0.0))
end

function shoot_v2(bal; m_total = 0.0, L_tube = 0.20, d_slider = 0.10, ζ = 0.002)
    if m_total <= 0
        return shoot_rifle(V0_C; m_tuner = 0.0, h_bore = 0.0254, CFG_V2...,
                           ζ1 = ζ, ζ2 = ζ, p_of_t = bal.p, x_of_t = bal.x,
                           t_b_override = bal.t_b)
    end
    b = tube_budget(m_total, L_tube)
    return shoot_rifle(V0_C; m_tuner = 0.0, h_bore = 0.0254, CFG_V2...,
                       ζ1 = ζ, ζ2 = ζ, p_of_t = bal.p, x_of_t = bal.x,
                       t_b_override = bal.t_b,
                       L_tube = L_tube, m_slider = b.m_slider, d_slider = d_slider)
end

function main()
    bal = coupled_ballistics(V0_C)
    println("="^78)
    println(" CADRE v2 — structure fusil × balistique couplée × tube élastique")
    println("="^78)
    @printf("\nBalistique : t_b = %.3f ms, v = %.0f m/s, pic = %.1f MPa\n",
            bal.t_b*1e3, bal.v_b, bal.p_peak/1e6)
    b = tube_budget(0.200, 0.20)
    @printf("Tuner 200 g, tube alu 200 mm : tube %.0f g + curseur %.0f g\n\n",
            b.m_tube*1e3, b.m_slider*1e3)

    nu = shoot_v2(bal)
    @printf("Canon nu : θ̇(t_b) = %+.2f MOA/ms, pic %.2f   (Kolbe : −9,4)\n\n",
            nu.θdot_tb, nu.θdot_peak)

    println("Réglage en POSITION du curseur (tube 200 mm, 200 g au total) :")
    println("  position |  θ̇(t_b)  |  pic  | écart au canon nu")
    println("  " * "-"^52)
    vals = Float64[]
    for d in 0.02:0.03:0.20
        r = shoot_v2(bal; m_total = 0.200, d_slider = d)
        r === nothing && continue
        push!(vals, r.θdot_tb)
        @printf("  %6.0f mm | %+8.2f | %5.2f | %+17.2f\n",
                d*1e3, r.θdot_tb, r.θdot_peak, r.θdot_tb - nu.θdot_tb)
    end

    println()
    println("-"^78)
    if !isempty(vals)
        @printf("Débattement obtenu par la position : %.2f MOA/ms   (Kolbe : 15,4)\n",
                maximum(vals) - minimum(vals))
        @printf("Plage atteinte : %+.2f à %+.2f            (Kolbe : −9,4 à +6,0)\n",
                minimum(vals), maximum(vals))
    end
    println("-"^78)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
