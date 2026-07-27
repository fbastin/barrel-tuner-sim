# =============================================================================
# CARTE DE SENSIBILITÉ DU VERDICT DE COMPENSATION  —  (L, ζ, h_offset)
#
# POURQUOI CE SCRIPT EXISTE. Entre le 19 et le 20 juillet 2026, la conclusion
# publiée sur la compensation positive a changé trois fois — « atteinte »,
# « atteinte à aucun réglage », « atteinte si le canon fait 28" ». À chaque fois
# le modèle n'avait pas changé : c'est un paramètre non mesuré qui avait bougé.
# Trois d'entre eux basculent le verdict INDÉPENDAMMENT :
#     L         longueur de canon    — choix jamais justifié (26" vs 28")
#     ζ         amortissement modal  — non mesuré ; les 4 estimations publiées
#                                      encadrent le seuil de bascule
#     h_offset  bras de levier       — CALÉ, ~6× la seule cote mesurée (Vaughn)
# Ce script cartographie le verdict sur les trois, pour répondre à la seule
# question qui vaille : ce verdict est-il un résultat du modèle, ou l'ombre d'un
# choix de paramètre ?
#
# ASTUCE QUI REND LA CARTE BON MARCHÉ. θ̇ est LINÉAIRE en h_offset (l'excitation
# lui est proportionnelle et la structure est linéaire). h_offset n'a donc pas
# besoin d'être un axe : on balaie (L, ζ) à un h_offset de référence, puis on
# INVERSE — pour chaque cellule, quel bras de levier faudrait-il pour atteindre
# la cible ? C'est une grandeur continue et interprétable, là où un verdict
# oui/non perd toute l'information.
#
# -----------------------------------------------------------------------------
# RÉSULTAT (2026-07-20). Bras de levier REQUIS, en multiples de la cote physique
# de Vaughn (~10,5 mm), seule valeur adossée à une mesure :
#
#            ζ=0,5%   ζ=1%    ζ=1,5%   ζ=2%    ζ=3%
#   610 mm    7,5×    8,1×     9,3×   10,6×   13,4×
#   660 mm    4,5×    6,4×     8,3×    9,9×   12,6×
#   711 mm    3,6×    5,5×     7,7×    9,9×   13,4×
#   762 mm    3,4×    4,9×     7,1×    9,4×   13,8×
#
# LE RÉSULTAT EST UN PLANCHER, ET IL EST ROBUSTE. Sur TOUT le domaine plausible,
# le modèle exige entre 3,4× et 13,8× le bras de levier physique. Il n'existe
# AUCUNE combinaison (L, ζ) où la compensation soit atteinte avec une excitation
# physiquement défendable. Ce n'est donc pas une question de réglage : c'est un
# MÉCANISME D'EXCITATION QUI MANQUE — très probablement la rotation de corps
# rigide de l'arme sous le recul, que l'ossature encastrée ne peut pas produire.
#
# CE QUI EST INDÉCIDABLE, ET CE QUI NE L'EST PAS.
#   * Indécidable : « le modèle atteint-il la cible ? » à h_offset = 62,5 mm.
#     Le verdict bascule sur L ET sur ζ dans leurs plages plausibles —
#     OUI à (762, 0,5 %), non à (610, 3 %). Mais cette question est CIRCULAIRE :
#     62,5 mm est précisément le nombre qui avait été calé pour que ça marche.
#   * Décidé, et négativement : « le modèle reproduit-il la compensation avec une
#     excitation physique ? » Non, partout, d'un facteur 3,4 au moins.
# La seconde formulation est la bonne, et c'est la seule à publier.
#
# EFFET DE ζ CONTRE EFFET DE L. ζ domine : à L fixé, passer de 0,5 % à 3 %
# multiplie l'exigence par 3,5 à 4. L aide, mais deux fois moins (facteur ~2,2 à
# ζ faible) — et son bénéfice S'ÉVANOUIT à fort amortissement (140 mm à 610 contre
# 145 à 762 pour ζ = 3 % : l'ordre s'inverse). Conséquence pratique : mesurer ζ
# prime sur trancher L, et trancher L à ζ inconnu ne règle rien.
#
# ATTENTION AUX COTES PUBLIÉES. Le porte-à-faux optimal se déplace énormément
# avec L : 170 mm à 610, 110 à 660, 70 à 711, 40 à 762 (tuner de 100 g). Toute
# republication de cotes suppose L tranché — et ζ, qui les déplace aussi de
# 10 à 30 mm.
#
# Usage :  julia sensitivity_map.jl
# =============================================================================

include(joinpath(@__DIR__, "simulation.jl"))

const H_REF   = 0.0625   # bras de levier calé (valeur publiée), référence linéaire
const H_PHYS  = 0.0105   # cote physique de Vaughn (9-12 mm), milieu de fourchette
const M_PROBE = 0.100    # tuner de 100 g
const ZETAS   = [0.005, 0.01, 0.015, 0.02, 0.03]
const D_GRID  = 0.0:0.01:0.20

# NOTE : L est une `const` de simulation.jl. Ce script cartographie donc ζ à la
# longueur COMPILÉE, et documente les autres longueurs depuis l'en-tête ci-dessus.
# Pour refaire la carte complète, relancer en éditant `const L` (les valeurs de
# l'en-tête ont été obtenues ainsi, par variantes du fichier).

function peak_thetadot(ζ)
    θmax = -Inf; dbest = 0.0
    for d in D_GRID
        r = simulate_shot(M_PROBE; d_overhang = d, ζ1 = ζ, ζ2 = ζ,
                          h_offset = H_REF, verbose = false)
        if r.θdot_MOAms > θmax
            θmax = r.θdot_MOAms; dbest = d
        end
    end
    return θmax, dbest
end

function main()
    cible = θdot_optimum_MOAms
    @printf("Longueur compilée L = %.0f mm — cible dérivée %.2f MOA/ms (τ_v = %.2f µs/(m/s))\n\n",
            L * 1e3, cible, projectile_kinematics(v_muzzle, L).τ_v * 1e6)
    @printf("%6s | %10s | %8s | %12s | %10s\n",
            "ζ", "θ̇ max", "d optim", "h_offset req", "× physique")
    println("-"^60)
    for ζ in ZETAS
        θmax, dbest = peak_thetadot(ζ)
        h_req = H_REF * cible / θmax          # linéarité en h_offset
        @printf("%5.1f%% | %10.2f | %6.0f mm | %9.0f mm | %9.1f×\n",
                ζ * 100, θmax, dbest * 1e3, h_req * 1e3, h_req / H_PHYS)
    end
    println("-"^60)
    println("Lecture : « h_offset req » est le bras de levier qu'il FAUDRAIT pour")
    println("atteindre la cible. La cote physique mesurée vaut ~10,5 mm (Vaughn).")
    println("Un facteur > 1 mesure ce que le modèle doit inventer pour y arriver.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
