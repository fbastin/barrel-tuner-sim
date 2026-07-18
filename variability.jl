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
# Le premier terme dépend de l'instant de sortie t_b = (1+φ)·L/v, le second de
# la vitesse seule ; leur annulation EST la compensation positive.
#
# APPROXIMATION ASSUMÉE. L'historique θ(t) est calculé UNE fois (deux passes
# Newmark, cf. plus bas) puis échantillonné à des t_b différents, au lieu d'être
# recalculé pour chaque v. C'est licite au premier ordre : l'excitation est
# chamber_pressure(t), fonction du TEMPS et non de la vitesse. Seule la charge
# mobile du projectile dépend de v — terme minuscule (2,6 g) isolé ci-dessous.
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
function response_split(m_tuner; d_overhang = 0.0, h_offset = H_OFFSET_EFF)
    tot  = simulate_shot(m_tuner; d_overhang, h_offset, verbose = false)
    proj = simulate_shot(m_tuner; d_overhang, h_offset = 0.0, verbose = false)
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
function dispersion(resp; σ_v = σ_V_DEFAULT, σ_k = σ_K_DEFAULT,
                    D = D_target, v0 = v_muzzle, n = N_SHOTS, seed = SEED)
    rng = MersenneTwister(seed)
    ys  = Vector{Float64}(undef, n)
    for i in 1:n
        δv = σ_v * randn(rng)
        k  = 1.0 + σ_k * randn(rng)
        v   = v0 + δv
        t_b = (1 + PHI_BURN) * L / v
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
configs = [
    ("canon nu (sans tuner)",        0.0,   0.0),
    ("accordé (100 g à 80 mm)",      0.100, 0.080),
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
nu, ac = results["canon nu (sans tuner)"], results["accordé (100 g à 80 mm)"]
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
scans = [(0.100, 0.055:0.005:0.130, 0.080), (0.200, 0.035:0.005:0.110, 0.065)]
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
    bests[m] = (; best = b, pub)
    @printf("Minimum de dispersion à %.0f mm (%.2f mm) contre %.2f mm au réglage publié\n",
            b.d * 1e3, b.sd, pub.sd)
    @printf("de %.0f mm — facteur %.1f. L'angle absolu y passe de %.0f à %.0f µrad.\n\n",
            pub.d * 1e3, pub.sd / b.sd, abs(pub.θ) * 1e6, abs(b.θ) * 1e6)
end

best = bests[0.100].best
best_d, best_sd = best.d, best.sd
println("LE CRITÈRE, ET POURQUOI LES DEUX MASSES NE RÉAGISSENT PAS PAREIL")
println("-"^78)
println("θ et θ̇ sont en QUADRATURE : θ̇ est maximal là où θ traverse zéro. Les deux")
println("critères ne s'opposent donc pas, ils coïncident — le réglage qui compense le")
println("mieux est aussi le moins sensible à l'aléa d'excitation. Le critère unifié")
println("s'énonce : faire sortir la balle quand la bouche est à son angle NEUTRE.")
println()
println("Reste que viser θ̇ = 6,0 EXACTEMENT n'y conduit pas toujours :")
println("  • à 200 g, le modèle plafonne à 5,91 et n'atteint jamais 6,0 ; « au plus")
println("    proche » retombe donc sur le maximum, tout près du passage par zéro.")
println("    Le 65 mm publié est déjà optimal — facteur 1,0, rien à corriger.")
println("  • à 100 g, θ̇ franchit 6,0 dès 80 mm puis continue de monter jusqu'à 6,18.")
println("    Viser la valeur nominale arrête donc AVANT le maximum, et laisse")
println("    l'angle absolu à 131 µrad au lieu de 16. D'où le facteur 3.")
println()
println("PRÉDICTION TESTABLE : viser le passage de θ par zéro (ou, ce qui revient au")
println("même, le MAXIMUM de θ̇) plutôt qu'une valeur nominale de θ̇. La différence")
println("n'apparaît que lorsque la courbe dépasse la cible — cas du tuner léger.")
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
