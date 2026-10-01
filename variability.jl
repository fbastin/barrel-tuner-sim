# =============================================================================
# variability.jl — Dispersion prédite quand l'excitation NE se répète PAS
#
# POURQUOI CE SCRIPT
# simulation.jl applique coup après coup exactement le même moment de recul :
# une même cartouche y produit une même vibration, donc un même angle de bouche.
# Vaughn (Rifle Accuracy Facts, 1998, ch. 4) mesure qu'il n'en est rien — sur
# une .270 Win instrumentée à l'anneau de culasse, le moment crête varie de
# 300 à 600 in-lb AVEC LE MÊME LOT, soit ~±30 %, et c'est cette variabilité,
# non la vibration moyenne, qui produit chez lui ~0,8" de dispersion à 100 yd.
#
# CE QUE ÇA CHANGE, ET C'EST LE POINT
# La compensation positive fonctionne parce que l'angle de bouche à la sortie
# est CORRÉLÉ à la vitesse : une balle lente sort plus tard, la bouche a monté
# davantage, le tir plus haut compense la chute supplémentaire. La variabilité
# d'amplitude, elle, est DÉCORRÉLÉE de la vitesse : aucun réglage de tuner ne
# peut l'annuler. Elle constitue donc un PLANCHER, au même titre que la
# dispersion propre de la munition.
#
# Ce script chiffre ce plancher, et le compare au gain que le tuner retire
# effectivement. Il produit une dispersion en cible (mm), directement
# comparable à un groupement mesuré — ce que le modèle déterministe ne savait
# pas faire, ne produisant qu'un angle.
#
# MODÈLE
# Pour chaque coup on tire deux aléas indépendants :
#   δv ~ N(0, σ_v)   dispersion de vitesse du lot (chronographe)
#   k  ~ N(1, σ_k)   facteur d'échelle de l'excitation (Vaughn)
# et l'on forme la hauteur d'impact, relative au coup nominal :
#   y = D·θ(t_b(v)) + g·D²·δv/v₀³
#       └─ angle de lancement ─┘   └─ chute différentielle ─┘
# Le premier terme dépend de l'instant de sortie t_b(v), le second de la vitesse
# seule ; leur annulation EST la compensation positive.
#
# APPROXIMATION ASSUMÉE. L'historique θ(t) est calculé UNE fois (deux passes
# Newmark, cf. plus bas) puis échantillonné à des t_b différents, au lieu d'être
# recalculé pour chaque v. C'est licite au premier ordre : l'excitation est la
# pression p(t) issue de la balistique intérieure, fonction du TEMPS et non de la
# vitesse. Seule la charge mobile du projectile dépend de v — terme minuscule
# (2,6 g) isolé ci-dessous.
#
# Usage :   julia variability.jl
# =============================================================================

include("simulation.jl")

using Random
using Statistics

# -----------------------------------------------------------------------------
# 1. PARAMÈTRES DE VARIABILITÉ
# -----------------------------------------------------------------------------
# σ_k : Vaughn rapporte une plage de 300 à 600 in-lb autour de 450, établie sur
# « several hundred records ». Une plage quasi complète se lit ~±3σ, d'où
# σ_k ≈ 10 %. C'est l'hypothèse RETENUE, et la plus prudente des trois :
# lire ±30 % comme ±2σ (σ = 15 %) ou comme ±1σ gonflerait le plancher. Le
# balayage final montre la sensibilité à ce choix.
const σ_K_DEFAULT = 0.10

# σ_v : l'écart-type de vitesse du lot. Pour une série de ~20 coups, ES ≈ 3,7σ ;
# l'ES de 10 m/s utilisé comme exemple de référence sur le site correspond donc
# à σ_v ≈ 2,7 m/s.
const σ_V_DEFAULT = 2.7

const N_SHOTS = 20_000
const SEED    = 20260718

# -----------------------------------------------------------------------------
# 2. DÉCOMPOSITION DE LA RÉPONSE
#
# Le système est linéaire, et l'excitation se compose de deux termes dont UN
# SEUL porte h_offset :
#     θ(t)  =  k · θ_rec(t)  +  θ_proj(t)
# On les sépare par deux passes : une avec h_offset nominal (total), une avec
# h_offset = 0 (charge mobile seule). C'est ce qui rend le Monte-Carlo gratuit —
# 2 résolutions de Newmark au lieu de N_SHOTS.
# -----------------------------------------------------------------------------
function response_split(m_tuner; d_overhang = 0.0, h_offset = H_OFFSET_EFF,
                        J_tuner = tuner_inertia(m_tuner))
    tot  = simulate_shot(m_tuner; d_overhang, h_offset, J_tuner, verbose = false)
    proj = simulate_shot(m_tuner; d_overhang, h_offset = 0.0, J_tuner, verbose = false)
    return (ts = tot.ts, θ_rec = tot.θ_L .- proj.θ_L, θ_proj = proj.θ_L)
end

# Interpolation linéaire de θ à un instant quelconque (pas de temps constant).
function interp(ts, ys, t)
    t <= ts[1]   && return ys[1]
    t >= ts[end] && return ys[end]
    Δ = ts[2] - ts[1]
    i = clamp(floor(Int, (t - ts[1]) / Δ) + 1, 1, length(ts) - 1)
    w = (t - ts[i]) / Δ
    return (1 - w) * ys[i] + w * ys[i + 1]
end

# -----------------------------------------------------------------------------
# 3. MONTE-CARLO
# -----------------------------------------------------------------------------
# Renvoie l'écart-type et l'extrême spread de la hauteur d'impact, en mm.
# La graine est FIXE et commune à tous les appels : les configurations comparées
# subissent ainsi exactement les mêmes tirages (variables aléatoires communes),
# de sorte que les écarts entre lignes du tableau viennent du réglage et non du
# bruit d'échantillonnage.
const KIN_REF = projectile_kinematics(v_muzzle, L)

function dispersion(resp; σ_v = σ_V_DEFAULT, σ_k = σ_K_DEFAULT,
                    D = D_target, v0 = v_muzzle, n = N_SHOTS, seed = SEED)
    rng = MersenneTwister(seed)
    ys  = Vector{Float64}(undef, n)
    for i in 1:n
        δv = σ_v * randn(rng)
        k  = 1.0 + σ_k * randn(rng)
        v   = v0 + δv
        # t_b(v) : linéarisation autour du nominal via τ_v, DÉRIVÉE du modèle
        # couplé (section 6 de simulation.jl). Remplace l'ancienne formule
        # burnout (1+φ)L/v, supprimée avec la refonte du 2026-07-19 — c'est
        # cette référence morte qui cassait ce script. Recalculer la balistique
        # intérieure complète pour chacun des 20 000 tirs serait prohibitif, et
        # inutile : τ_v est par définition ∂t_b/∂v, donc exact au premier ordre
        # sur des écarts de vitesse de quelques m/s.
        t_b = KIN_REF.t_b - KIN_REF.τ_v * δv
        θ   = k * interp(resp.ts, resp.θ_rec, t_b) + interp(resp.ts, resp.θ_proj, t_b)
        # angle de lancement + chute différentielle (balle rapide = impact haut)
        ys[i] = (D * θ + g_accel * D^2 * δv / v0^3) * 1e3
    end
    return (sd = std(ys), es = maximum(ys) - minimum(ys))
end

# -----------------------------------------------------------------------------
# 4. EXÉCUTION
# -----------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__

println("="^78)
println(" Dispersion prédite avec une excitation NON reproductible")
println(" Variabilité mesurée par Vaughn (1998, ch. 4) : ±30 % sur le moment crête")
println("="^78)
println()
@printf("Distance %.0f m | v₀ %.0f m/s | σ_v = %.1f m/s (ES ≈ %.0f) | σ_k = %.0f %%\n",
        D_target, v_muzzle, σ_V_DEFAULT, 3.7 * σ_V_DEFAULT, 100 * σ_K_DEFAULT)
@printf("%d tirs simulés par configuration, graine %d.\n", N_SHOTS, SEED)
println()

# Deux configurations : canon nu, et l'accord retenu dans la documentation.
# Les libellés servent de clés plus bas : on les nomme une seule fois.
# Réglage de référence : l'optimum COURANT du modèle, non une cote figée.
# Il a bougé plusieurs fois (100 → 110 → 135 mm) au fil des révisions du
# 2026-07-19 ; le coder en dur laissait ce script décrire un réglage périmé.
const M_DOC = 0.100
# Le critère « au plus proche de θ̇ = 6,0 » appliqué au balayage. Rendu FONCTION
# de la masse : la cote du tuner 200 g était restée un littéral en dur (0.085)
# pendant que celle du 100 g se recalculait, si bien qu'une révision du modèle
# n'en corrigeait qu'une des deux.
function d_criterion(m; ds = 0.0:0.005:0.20)
    rs = [simulate_shot(m; d_overhang = d, h_offset = H_OFFSET_EFF,
                        verbose = false).θdot_MOAms for d in ds]
    ds[argmin([abs(x - θdot_optimum_MOAms) for x in rs])]
end
const D_DOC = d_criterion(M_DOC)
const LBL_BARE  = "canon nu (sans tuner)"
const LBL_TUNED = @sprintf("accordé (%d g à %d mm)", round(Int, M_DOC*1e3), round(Int, D_DOC*1e3))
configs = [
    (LBL_BARE,  0.0,   0.0),
    (LBL_TUNED, M_DOC, D_DOC),
]

println("-"^78)
@printf("%-26s | %-22s | %8s | %8s\n", "configuration", "source d'aléa", "SD (mm)", "ES (mm)")
println("-"^78)

results = Dict{String,Any}()
for (label, m, d) in configs
    resp = response_split(m; d_overhang = d)
    only_v = dispersion(resp; σ_k = 0.0)                    # vitesse seule
    only_k = dispersion(resp; σ_v = 0.0)                    # amplitude seule
    both   = dispersion(resp)                               # les deux
    results[label] = (; only_v, only_k, both)
    @printf("%-26s | %-22s | %8.2f | %8.2f\n", label, "vitesse seule",   only_v.sd, only_v.es)
    @printf("%-26s | %-22s | %8.2f | %8.2f\n", "",    "amplitude seule", only_k.sd, only_k.es)
    @printf("%-26s | %-22s | %8.2f | %8.2f\n", "",    "les deux",        both.sd,   both.es)
    println("-"^78)
end

println()
println("LECTURE")
nu, ac = results[LBL_BARE], results[LBL_TUNED]
@printf("  • Le tuner écrase la composante de VITESSE : %.2f → %.2f mm d'écart-type\n",
        nu.only_v.sd, ac.only_v.sd)
@printf("    (c'est la compensation positive, seul mécanisme que ce modèle décrit).\n")
@printf("  • La composante d'AMPLITUDE, elle, n'est pas compensée : %.2f → %.2f mm.\n",
        nu.only_k.sd, ac.only_k.sd)
@printf("    Elle baisse tout de même, mais par un tout autre chemin — non parce que\n")
@printf("    le tuner l'annule (décorrélée de la vitesse, elle échappe par construction\n")
@printf("    à la compensation), mais parce que ce réglage se trouve à un angle absolu\n")
@printf("    plus faible. C'est une propriété de la POSITION choisie, pas du principe.\n")
@printf("  • Dispersion résiduelle après accord : %.2f mm d'écart-type, dont\n", ac.both.sd)
@printf("    %.0f %% imputables à la seule variabilité de l'excitation.\n",
        100 * ac.only_k.sd^2 / ac.both.sd^2)
println()

# ---------------------------------------------------------------------------
# PRÉDICTION NOUVELLE, et c'est l'apport de ce script.
#
# La compensation positive ne contraint que la DÉRIVÉE θ̇(t_b). La dispersion
# due à la variabilité d'excitation, elle, est proportionnelle à l'angle
# ABSOLU θ(t_b) — puisque y_ampl = D·θ·σ_k. Ces deux critères sont
# indépendants : parmi les réglages qui satisfont θ̇ ≈ 6 MOA/ms, ceux dont
# |θ| est petit sont strictement meilleurs. Le modèle déterministe ne pouvait
# pas voir cette hiérarchie, n'ayant qu'un seul critère.
# ---------------------------------------------------------------------------
scan(m, ds) = map(ds) do d
    res = simulate_shot(m; d_overhang = d, h_offset = H_OFFSET_EFF, verbose = false)
    (d = d, θdot = res.θdot_MOAms, θ = res.θ_at_tb,
     sd = dispersion(response_split(m; d_overhang = d)).sd)
end

# Les deux masses documentées, et leur réglage publié au critère θ̇ seul.
# Réglages et plages révisés le 2026-07-19 : avec J dérivé de la masse, le
# réglage au critère θ̇ passe de 90 à 100 mm à 100 g (inchangé à 65 mm à 200 g,
# la section du tube étant calée pour redonner k = 5,02 cm à cette masse).
# Les plages sont étendues en conséquence pour encadrer le nœud, qui recule lui
# aussi (110 → 120 mm à 100 g).
scans = [(0.100, 0.080:0.005:0.180, D_DOC),
         (0.200, 0.040:0.005:0.130, d_criterion(0.200; ds = 0.040:0.005:0.130))]
bests = Dict{Float64,Any}()

for (m, ds, d_pub) in scans
    @printf("OÙ ACCORDER, UNE FOIS L'ALÉA PRIS EN COMPTE  (tuner %d g)\n", round(Int, m*1e3))
    println("-"^78)
    @printf("%14s | %12s | %12s | %10s | %10s\n",
            "porte-à-faux", "θ̇ (MOA/ms)", "θ (µrad)", "SD (mm)", "critère θ̇")
    println("-"^78)
    rows = scan(m, ds)
    for r in rows
        @printf("%11.0f mm | %+12.3f | %+12.1f | %10.2f | %8s%s\n",
                r.d * 1e3, r.θdot, r.θ * 1e6, r.sd,
                abs(r.θdot - θdot_optimum_MOAms) < 1.0 ? "✓" : "",
                abs(r.d - d_pub) < 1e-9 ? "  ← publié" : "")
    end
    println("-"^78)
    b = rows[argmin([r.sd for r in rows])]
    pub = rows[argmin([abs(r.d - d_pub) for r in rows])]
    bests[m] = (; best = b, pub, rows)
    @printf("Minimum de dispersion à %.0f mm (%.2f mm) contre %.2f mm au réglage publié\n",
            b.d * 1e3, b.sd, pub.sd)
    @printf("de %.0f mm — facteur %.1f. L'angle absolu y passe de %.0f à %.0f µrad.\n\n",
            pub.d * 1e3, pub.sd / b.sd, abs(pub.θ) * 1e6, abs(b.θ) * 1e6)
end

best = bests[0.100].best
best_d, best_sd = best.d, best.sd
println("LE CRITÈRE, ET POURQUOI LES DEUX NE COÏNCIDENT PAS")
println("-"^78)
println("Une version antérieure affirmait ici que les deux critères — viser θ̇ maximal,")
println("viser l'angle neutre θ = 0 — COÏNCIDENT, « θ et θ̇ étant en quadrature ». La")
println("quadrature est réelle, mais elle porte sur le TEMPS : à porte-à-faux fixé,")
println("θ̇(t) est bien maximal quand θ(t) passe par zéro. Elle ne dit RIEN du balayage")
println("en POSITION, qui change à la fois l'amplitude et la phase du système. Les")
println("tables ci-dessus le montrent : θ̇ culmine bien avant que θ ne s'annule.")
println("C'était une conflation entre deux quadratures, l'une temporelle et vraie,")
println("l'autre positionnelle et fausse.")
println()
println("Les chiffres ci-dessous sont RECALCULÉS à chaque exécution : les figer en dur")
println("les aurait laissés mentir à la première révision du modèle (ce qui est arrivé).")
@printf("La cible est désormais DÉRIVÉE de la cinématique du modèle : %.2f MOA/ms\n",
        θdot_optimum_MOAms)
@printf("(contre %.1f repris de Kolbe, soit %+.0f %%). Elle n'est atteignable à AUCUN\n",
        θdot_KOLBE_MOAms, 100*(θdot_optimum_MOAms/θdot_KOLBE_MOAms-1))
println("réglage — θ̇ plafonne en dessous — de sorte que « au plus proche » retombe sur")
println("le MAXIMUM de θ̇ aux deux masses. Ce n'est pas un artefact : c'est le modèle")
println("qui dit que la compensation n'est pas atteinte à ce h_offset.")
let b2 = bests[0.200], b1 = bests[0.100]
    # La formulation elle-même est DÉRIVÉE du balayage. Elle affirmait « θ̇
    # plafonne et n'atteint jamais 6,0 » : vrai de l'ancien modèle, faux du
    # nouveau, où θ̇ culmine au-delà de 6,0 et redescend. Seuls les nombres
    # étaient recalculés, si bien qu'une phrase fausse encadrait des chiffres
    # justes — et 5,98, qui n'est que la valeur AU POINT PUBLIÉ, se lisait
    # comme un maximum.
    θdot_max2 = maximum(r.θdot for r in b2.rows)
    if θdot_max2 <= θdot_optimum_MOAms
        @printf("  • à 200 g, θ̇ plafonne à %.2f et n'atteint jamais la cible de %.2f ;\n",
                θdot_max2, θdot_optimum_MOAms)
        @printf("    « au plus proche » retombe donc sur le MAXIMUM de θ̇, à %.0f mm — et non\n",
                b2.pub.d * 1e3)
        println("    sur le nœud, qui est ailleurs.")
    else
        @printf("  • à 200 g, θ̇ culmine à %.2f vers %.0f mm, DÉPASSE donc la cible, puis\n",
                θdot_max2, b2.rows[argmax([r.θdot for r in b2.rows])].d * 1e3)
        @printf("    redescend ; « au plus proche de 6,0 » retient %.0f mm sur la branche\n",
                b2.pub.d * 1e3)
        @printf("    DESCENDANTE, où θ̇ = %.2f et l'angle n'est plus qu'à %.0f µrad du neutre.\n",
                b2.pub.θdot, abs(b2.pub.θ) * 1e6)
    end
    let r2 = b2.pub.sd / b2.best.sd
        if r2 < 1.1
            @printf("    Le %.0f mm publié est déjà optimal — facteur %.1f, rien à corriger.\n",
                    b2.pub.d * 1e3, r2)
        else
            @printf("    Le minimum de dispersion est ailleurs, à %.0f mm : facteur %.1f.\n",
                    b2.best.d * 1e3, r2)
        end
    end
    @printf("  • à 100 g, le minimum de dispersion est à %.0f mm (θ̇ = %.2f) contre %.0f mm\n",
            b1.best.d*1e3, b1.best.θdot, b1.pub.d*1e3)
    @printf("    au critère θ̇ : %.0f µrad du neutre au lieu de %.0f, soit un facteur %.1f.\n",
            abs(b1.pub.θ) * 1e6, abs(b1.best.θ) * 1e6, b1.pub.sd / b1.best.sd)
    println()
    println("  ⚠️ HISTOIRE DE CE FACTEUR, QUI A FAILLI ÊTRE ENTERRÉ À TORT. Il valait 2,8,")
    println("  puis a semblé s'effondrer à 1,1 après la refonte en balistique intérieure")
    println("  couplée — au point qu'on a publié que son « enjeu chiffré avait disparu ».")
    println("  C'ÉTAIT UN ARTEFACT. La cible restait figée à 6,0, empruntée à Kolbe, quand")
    println("  la cinématique du modèle en implique ~6,9 : le critère tombait par hasard")
    @printf("  près du nœud. Cible rendue auto-cohérente, le facteur remonte à %.1f à 100 g\n",
            b1.pub.sd / b1.best.sd)
    @printf("  et %.1f à 200 g, du même ordre que le 2,8 d'origine. La recommandation\n",
            b2.pub.sd / b2.best.sd)
    println("  « viser le nœud plutôt que le chiffre » retrouve donc sa force, et pour la")
    println("  bonne raison cette fois : les deux critères ne coïncident PAS (ci-dessus).")
end
println()
println("PRÉDICTION TESTABLE : viser le passage de θ par zéro, et non le maximum de θ̇")
println("ni une valeur nominale de θ̇ — les deux critères ne coïncident pas (ci-dessus).")
println()

# Sensibilité à l'hypothèse sur σ_k — le paramètre le moins assuré.
println("SENSIBILITÉ À L'INTERPRÉTATION DES ±30 % DE VAUGHN")
println("-"^78)
@printf("%-34s | %10s | %14s\n", "lecture de la plage 300-600 in-lb", "σ_k", "SD accordé (mm)")
println("-"^78)
resp_ac = response_split(0.100; d_overhang = 0.080)
for (lbl, σk) in (("±30 % ≈ ±3σ  (retenu, prudent)", 0.10),
                  ("±30 % ≈ ±2σ", 0.15),
                  ("±30 % ≈ ±1σ  (majorant)", 0.30))
    @printf("%-34s | %9.0f %% | %14.2f\n", lbl, 100σk, dispersion(resp_ac; σ_k = σk).sd)
end
println("-"^78)
println()
println("RÉSERVE. Les ±30 % sont mesurés sur une carabine CENTERFIRE de chasse ;")
println("leur magnitude ne se transpose pas telle quelle à la .22 LR de match.")
println("C'est le mécanisme — une excitation qui ne se répète pas — qui se")
println("transpose, et le plancher qu'il impose à tout accord de tuner.")
println("="^78)

end
