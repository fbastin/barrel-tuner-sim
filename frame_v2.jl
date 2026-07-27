# =============================================================================
# CADRE v2 — MODÈLE UNIFIÉ DE RÉFÉRENCE (physique la plus complète du dossier)
#
# Rôle dans le dossier (fixé le 21 juillet 2026). C'est le modèle de RÉFÉRENCE au
# sens de la physique la plus riche : il compose les trois pièces validées de
# l'enquête et sert de VALIDATEUR qui borne l'erreur du modèle encastré
# `simulation.jl`. Il ne le REMPLACE pas : `simulation.jl` reste la source des
# grandeurs PUBLIÉES (wiki, PDF, simulateur web) parce que, encastré et paramétré,
# ses chiffres sont propres, reproductibles et indépendants de toute géométrie
# d'arme non mesurée. Ce cadre-ci, lui, dépend d'une géométrie d'arme inconnue :
# ses nombres ne sont robustes qu'en BALAYAGE (`sensitivity_coupled.jl`), pas à
# géométrie figée (voir le garde-fou en pied de `main`). Référence de physique,
# donc, non référence de chiffres — les deux rôles sont distincts et assumés.
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
#   3. TUNER      — modélisé des DEUX façons, à masse totale égale : masse
#      ponctuelle RIGIDE déportée (idéalisation, réglée par le porte-à-faux) et
#      tube ÉLASTIQUE maillé en poutre avec curseur coulissant (réaliste, alu,
#      réglé par la position du curseur). Les deux tournent ici côte à côte pour
#      montrer que le cadre les porte l'un comme l'autre.
#
# CE QUE CE CADRE FAIT — ET NE FAIT PAS. Il DÉMONTRE que les trois pièces
# tournent ensemble, à une géométrie d'arme représentative (CFG_V2). Il ne
# QUANTIFIE PAS la de-autorité du tube : à géométrie FIGÉE, l'ordre rigide/tube
# n'est pas robuste et peut s'inverser (c'est le piège récurrent du dossier). La
# de-autorité chiffrée vient du BALAYAGE de géométrie de `sensitivity_coupled.jl`
# (au meilleur cas θ̇max 8,7 rigide → 6,1 tube), pas d'un tir unique.
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

# Un tir dans le cadre unifié, avec le tuner selon l'un des DEUX modèles :
#   :tube   — tube élastique maillé en poutre + curseur coulissant (réaliste) ;
#             le réglage en position est celui du curseur (d_slider).
#   :rigide — masse ponctuelle déportée affectée de son inertie (idéalisation) ;
#             le réglage en position est le porte-à-faux (d_overhang).
# À masse TOTALE égale : c'est la seule comparaison honnête, l'écart mesurant
# alors ce que l'idéalisation rigide gonfle. `d` est la position (curseur ou
# porte-à-faux selon le mode), en mètres.
function shoot_v2(bal; m_total = 0.0, mode = :tube, L_tube = 0.20, d = 0.10, ζ = 0.002)
    if m_total <= 0
        return shoot_rifle(V0_C; m_tuner = 0.0, h_bore = 0.0254, CFG_V2...,
                           ζ1 = ζ, ζ2 = ζ, p_of_t = bal.p, x_of_t = bal.x,
                           t_b_override = bal.t_b)
    end
    if mode === :rigide
        return shoot_rifle(V0_C; m_tuner = m_total, h_bore = 0.0254, CFG_V2...,
                           ζ1 = ζ, ζ2 = ζ, p_of_t = bal.p, x_of_t = bal.x,
                           t_b_override = bal.t_b, d_overhang = d)
    end
    b = tube_budget(m_total, L_tube)
    return shoot_rifle(V0_C; m_tuner = 0.0, h_bore = 0.0254, CFG_V2...,
                       ζ1 = ζ, ζ2 = ζ, p_of_t = bal.p, x_of_t = bal.x,
                       t_b_override = bal.t_b,
                       L_tube = L_tube, m_slider = b.m_slider, d_slider = d)
end

# Balayage en position d'un mode : renvoie le vecteur des θ̇(t_b) obtenus et le
# débattement (max − min), pour la comparaison rigide vs tube.
function position_sweep(bal, mode; m_total = 0.200, ζ = 0.002)
    vals = Tuple{Float64,Float64,Float64}[]   # (position, θ̇, pic)
    for d in 0.02:0.03:0.20
        r = shoot_v2(bal; m_total = m_total, mode = mode, d = d, ζ = ζ)
        r === nothing && continue
        push!(vals, (d, r.θdot_tb, r.θdot_peak))
    end
    return vals
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

    # Les DEUX idéalisations du tuner, dans le même cadre et à masse totale égale.
    sw = Dict(m => position_sweep(bal, m) for m in (:rigide, :tube))
    labels = ((:rigide, "masse ponctuelle RIGIDE (porte-à-faux)"),
              (:tube,   "tube ÉLASTIQUE + curseur (position curseur)"))

    for (mode, title) in labels
        println("Réglage en POSITION — tuner : ", title)
        println("  position |  θ̇(t_b)  |  pic  | écart au canon nu")
        println("  " * "-"^52)
        for (d, θ, pic) in sw[mode]
            @printf("  %6.0f mm | %+8.2f | %5.2f | %+17.2f\n",
                    d*1e3, θ, pic, θ - nu.θdot_tb)
        end
        println()
    end

    println("-"^78)
    for (mode, title) in labels
        vals = [θ for (_, θ, _) in sw[mode]]
        isempty(vals) && continue
        @printf("%-8s | débattement %.2f MOA/ms | plage %+.2f à %+.2f\n",
                mode === :rigide ? "RIGIDE" : "TUBE",
                maximum(vals) - minimum(vals), minimum(vals), maximum(vals))
    end
    println("Référence mesurée par Kolbe : 15,4 MOA/ms de débattement (−9,4 à +6,0).")
    println()
    println("⚠ NE PAS lire ces deux lignes comme une mesure de de-autorité du tube.")
    println("  Ce cadre fige UNE géométrie d'arme (CFG_V2), et à géométrie figée l'ordre")
    println("  rigide/tube n'est pas robuste — il peut même s'inverser, comme ici. C'est")
    println("  précisément le piège que le dossier documente. La de-autorité du tube est")
    println("  établie sur le BALAYAGE de géométrie (sensitivity_coupled.jl) : au meilleur")
    println("  cas, θ̇max passe de 8,7 (rigide) à 6,1 (tube). Ce cadre-ci DÉMONTRE que les")
    println("  trois pièces tournent ensemble ; il ne QUANTIFIE pas, une seule config ne le")
    println("  permet pas.")
    println("-"^78)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
