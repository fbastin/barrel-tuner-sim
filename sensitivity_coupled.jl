# =============================================================================
# LE MUR EST-IL UN ARTEFACT DE L'ENCASTREMENT ?  —  carte de sensibilité portée
# sur le modèle COUPLÉ (crosse + appuis unilatéraux + balistique intérieure).
#
# `sensitivity_map.jl` a établi, sur `simulation.jl` (canon ENCASTRÉ, tuner en
# masse ponctuelle), qu'atteindre la compensation exigerait de 3,4 à 13,8 fois le
# bras de levier physique — un « mur » structurel. Objection légitime : ce modèle
# OMET la rotation de corps rigide de l'arme sous le recul, c'est-à-dire le
# mécanisme même que Kolbe désigne comme la source de la vibration. Le mur
# pourrait n'être qu'un artefact de l'amputation.
#
# `rifle_coupled.jl` contient ce mécanisme (crosse, deux sacs à contact
# unilatéral, balistique couplée) mais ne balayait que la MASSE du tuner, jamais
# le porte-à-faux — le vrai levier d'accord. Ce script ajoute ce balayage et pose
# la MÊME question que la carte encastrée : à h_bore PHYSIQUE, quel multiple du
# bras de levier faudrait-il pour que θ̇(t_b) atteigne la cible dérivée ?
#
# MÉTHODE. θ̇(t_b) est linéaire en h_bore (vérifié : contact fermé sous la
# précharge de gravité, écart < 0,1 %). Pour chaque amortissement ζ, on optimise
# θ̇(t_b) sur (masse × porte-à-faux) à h_bore = 25,4 mm — hauteur d'âme au-dessus
# de la ligne des sacs, grandeur GÉOMÉTRIQUE réelle, sans le facteur d'échelle
# qu'exigeait l'encastrement — et le facteur requis vaut cible / θ̇max.
#
# DEUX RAFFINEMENTS, parce que la réponse change à chacun.
#   • Géométrie d'arme : la raideur de crosse, la position des appuis et de la
#     culasse sont INCONNUES (Kolbe ne les publie pas). À géométrie figée le facteur
#     reste O(10) ; dès qu'on la balaie il tombe à O(1). On rapporte donc, par ζ, la
#     fourchette entre la MEILLEURE et la PIRE géométrie de la grille.
#   • Modèle du tuner : masse ponctuelle RIGIDE (idéalisation) contre TUBE ÉLASTIQUE
#     (réaliste — cf. flexible_tube.jl : dans la bande 2-3 kHz le tube fléchit et ne
#     transmet plus le moment d'une masse rigide déportée). On calcule les DEUX cartes
#     à budget de masse égal, et on mesure de combien le tube resserre la frontière.
#
# -----------------------------------------------------------------------------
# RÉSULTAT (exécution du 2026-07-21 ; valeurs non recopiées en dur pour éviter
# qu'elles ne divergent du code — voir la sortie). Cible dérivée locale 7,24 MOA/ms,
# h_bore physique 25,4 mm. Facteur requis = cible / θ̇max, MEILLEURE géométrie :
#
#     ζ        tuner RIGIDE      tuner TUBE ÉLASTIQUE
#     0,003     0,8× (atteinte)   1,2×
#     0,005     0,9× (atteinte)   1,4×
#     0,010     1,2×              1,8×
#     0,020     1,7×              2,2×
#     0,030     2,1×              2,5×
#
# LECTURE — en trois temps, parce que la réponse change à chaque raffinement.
#
# (1) Le « mur » 3,4-13,8× de la carte ENCASTRÉE était en grande partie un ARTEFACT
# de l'amputation. Restaurer la rotation de corps rigide fait tomber le facteur
# requis à O(1) : le modèle encastré exigeait un bras EFFECTIF de 62,5 mm (6× le
# physique, sans réalité géométrique) et tombait ENCORE court ; le couplé opère
# avec un bras RÉEL de 25,4 mm. Ce n'était pas une paroi structurelle.
#
# (2) Avec un tuner RIGIDE, la meilleure géométrie ATTEINT même la cible à bras
# physique pour ζ ≲ 0,007. Mais c'est l'idéalisation qui l'y aide : une masse
# ponctuelle déportée transmet à la bouche un moment qu'un tube fléchissant ne
# transmet pas.
#
# (3) Avec le tuner en TUBE ÉLASTIQUE (réaliste : alu maillé en poutre, curseur
# coulissant, MÊME masse totale), la fenêtre se REFERME : la cible n'est atteinte
# à AUCUN ζ. La de-autorité est modeste sur le meilleur cas (×1,2-1,5, la liberté
# de géométrie récupérant l'essentiel du ×3 vu à géométrie figée), mais suffit à
# faire basculer : le manque passe de 1,2× (ζ=0,003, au mieux) à 2,5× (ζ=0,03).
#
# SYNTHÈSE. Ni « mur structurel », ni « compensation atteinte » : le modèle est
# ROBUSTEMENT, MODESTEMENT court. La conclusion publiée « atteint à aucun réglage
# physique » (pas ainsi) SURVIT au modèle le plus riche — mais sa QUANTIFICATION
# change : ce n'est pas un facteur 3,4-13,8× (spécifique à l'encastrement) mais un
# écart de 1,2 à 2,5× selon l'amortissement, assez petit pour tenir dans
# l'incertitude de modèle. La priorité « mesurer ζ » en ressort NUANCÉE : même au ζ
# le plus favorable le tube reste court de 1,2× — ζ ne suffit donc pas seul à
# ouvrir la compensation ; il faudrait que d'autres facteurs (modes supérieurs mal
# maillés, appui de crosse) jouent AUSSI dans le bon sens.
#
# Réserve : le pire cas de géométrie donne toujours θ̇(t_b) < 0 (mauvais signe,
# cible inatteignable par principe, pas par échelle) — la géométrie d'arme reste
# donc déterminante. Et rifle_coupled reproduit l'ORDRE de Kolbe, jamais son couple
# exact −9,4/+6,0 : incertitude de modèle réelle, du même ordre que l'écart de 1,2×.
#
# Usage :   julia sensitivity_coupled.jl    (~70-80 min : 2 modes × géométrie)
# =============================================================================

include(joinpath(@__DIR__, "rifle_coupled.jl"))
using Printf

const G_C  = 9.81
const D_C  = 50.0
const R2M  = 180/π*60
const H_PHYS_C = 0.0254            # hauteur d'âme physique au-dessus des sacs (m)

target_MOAms(τv_SI) = G_C * D_C / (V0_C^3 * τv_SI) * R2M / 1000

# Cible dérivée de la cinématique LOCALE (évite l'erreur de cible empruntée) :
# τ_v = |d t_b / d v| sur la balistique couplée de ce canon.
function derived_target()
    dv = 5.0
    tb_p = coupled_ballistics(V0_C + dv).t_b
    tb_m = coupled_ballistics(V0_C - dv).t_b
    τv   = abs(tb_p - tb_m) / (2dv)
    return target_MOAms(τv), τv
end

function assert_linearity(bal)
    kw = (p_of_t = bal.p, x_of_t = bal.x, t_b_override = bal.t_b)
    ref = nothing; worst = 0.0
    for hb in (0.0254, 0.0508, 0.1016)
        r = shoot_rifle(V0_C; m_tuner=0.2, h_bore=hb, d_overhang=0.065, CFG_C..., kw...)
        r === nothing && continue
        ref === nothing && (ref = r.θdot_tb / hb)
        worst = max(worst, abs(r.θdot_tb/hb/ref - 1))
    end
    return worst
end

# Grille de géométries d'arme plausibles (Kolbe ne publie pas la sienne) — même
# étendue que le balayage de rifle_coupled.main.
const XB = (0.35, 0.45, 0.55, 0.65)
const LF = (0.15, 0.25, 0.35)
const EI = (3e3, 1e4, 4e4)
const MASS = (0.1, 0.2, 0.3)
const OVER = 0.0:0.025:0.20

# Tuner en TUBE ÉLASTIQUE. build_rifle maille une poutre alu (Ø40×1,25 mm,
# TUBE_*_R) en avant de la bouche, avec un curseur coulissant. RÈGLE DU BUDGET DE
# MASSE (flexible_tube.jl) : à masse totale fixée, le tube consomme ρA·L_tube et
# le curseur reçoit le reste — sinon on comparerait deux objets de masses
# différentes. Tube alu de 200 mm ≈ 82 g, comme le Starik/Centra de référence.
const L_TUBE = 0.20
const M_TUBE = TUBE_RHOA_R * L_TUBE            # masse propre du tube (alu)

# θ̇(t_b) maximal sur (masse × réglage) pour une géométrie, un ζ et un MODE de
# tuner donnés. :rigid = masse ponctuelle déportée (porte-à-faux) ; :tube = tube
# élastique + curseur (position du curseur), à budget de masse égal.
function peak_over_tuner(cfg, ζ, kw; mode = :rigid)
    best = -Inf
    for m in MASS, d in OVER
        r = if mode === :rigid
            shoot_rifle(V0_C; m_tuner=m, h_bore=H_PHYS_C, d_overhang=d,
                        ζ1=ζ, ζ2=ζ, cfg..., kw...)
        else
            shoot_rifle(V0_C; m_tuner=0.0, h_bore=H_PHYS_C, L_tube=L_TUBE,
                        m_slider=max(m - M_TUBE, 0.0), d_slider=d,
                        ζ1=ζ, ζ2=ζ, cfg..., kw...)
        end
        r === nothing && continue
        r.θdot_tb > best && (best = r.θdot_tb)
    end
    return best
end

function main()
    bal = coupled_ballistics(V0_C)
    tgt, τv = derived_target()
    lin = assert_linearity(bal)
    kw = (p_of_t = bal.p, x_of_t = bal.x, t_b_override = bal.t_b)

    println("="^80)
    println(" Carte de sensibilité — modèle COUPLÉ (crosse + appuis + balistique)")
    println("="^80)
    @printf("\nBalistique : t_b = %.3f ms, v = %.1f m/s, pic = %.1f MPa\n",
            bal.t_b*1e3, bal.v_b, bal.p_peak/1e6)
    @printf("Cible dérivée (τ_v = %.2f µs/(m/s), local) : %.2f MOA/ms\n", τv*1e6, tgt)
    @printf("h_bore physique : %.1f mm  ;  linéarité θ̇/h : écart max %.3f %%\n\n",
            H_PHYS_C*1e3, 100*lin)

    @printf("Balayage géométrie d'arme (%d configs) × masse × réglage, par ζ et par mode.\n",
            length(XB)*length(LF)*length(EI))
    @printf("Tube alu : L = %.0f mm, masse propre %.0f g (curseur = masse totale − %.0f g).\n\n",
            L_TUBE*1e3, M_TUBE*1e3, M_TUBE*1e3)

    ζs = (0.003, 0.005, 0.010, 0.015, 0.020, 0.030)
    frontier = Dict{Symbol,Any}()          # meilleur cas θ̇max par ζ, pour la synthèse

    for mode in (:rigid, :tube)
        title = mode === :rigid ? "TUNER RIGIDE (masse ponctuelle déportée)" :
                                   "TUNER TUBE ÉLASTIQUE (alu, budget de masse égal)"
        println("── ", title, " ", "─"^(60 - length(title)))
        @printf("  %-7s | %-13s | %-11s | %-13s | %-11s | %s\n",
                "ζ", "θ̇max meilleur", "fact. min", "θ̇max pire", "fact. max", "cible atteinte ?")
        println("  " * "-"^80)
        bestθ_by_ζ = Float64[]
        for ζ in ζs
            best_θ = -Inf; worst_θ = Inf
            for xb in XB, lf in LF, ei in EI
                cfg = (x_breech=xb, L_fore=lf, x_rear=0.038, EI_stock=ei, k_rest=1e5)
                θ = peak_over_tuner(cfg, ζ, kw; mode = mode)
                θ == -Inf && continue
                best_θ  = max(best_θ, θ)
                worst_θ = min(worst_θ, θ)
            end
            push!(bestθ_by_ζ, best_θ)
            fmin = tgt / best_θ
            # Pire géométrie : θ̇max déjà négatif ⇒ cible inatteignable par SIGNE,
            # pas par échelle (le rapport tgt/θ sortirait négatif, sans sens).
            fmax_str = worst_θ > 0 ? @sprintf("%6.1f×", tgt / worst_θ) : "  signe ✗"
            reached  = best_θ >= tgt ? "OUI (bras physique)" : "non"
            @printf("  %-7.3f | %+-13.2f | %6.1f×     | %+-13.2f | %-9s | %s\n",
                    ζ, best_θ, fmin, worst_θ, fmax_str, reached)
        end
        frontier[mode] = bestθ_by_ζ
        println()
    end

    println("─ RESSERREMENT DE LA FRONTIÈRE (meilleur cas, facteur = cible/θ̇max) ", "─"^14)
    @printf("  %-7s | %-16s | %-16s | %s\n", "ζ", "rigide", "tube élastique", "de-autorité")
    println("  " * "-"^70)
    for (i, ζ) in enumerate(ζs)
        θr = frontier[:rigid][i]; θt = frontier[:tube][i]
        fr = θr > 0 ? @sprintf("%.1f×", tgt/θr) : "signe✗"
        ft = θt > 0 ? @sprintf("%.1f×", tgt/θt) : "signe✗"
        rr = θr >= tgt ? " (atteinte)" : ""
        tt = θt >= tgt ? " (atteinte)" : ""
        @printf("  %-7.3f | %+6.2f → %-7s%s | %+6.2f → %-7s%s | ×%.1f\n",
                ζ, θr, fr, rr, θt, ft, tt, θr/θt)
    end
    println()
    println("-"^80)
    println("Rappel carte ENCASTRÉE (sensitivity_map.jl) : 3,4 à 13,8× — mur uniforme.")
    println("-"^80)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
