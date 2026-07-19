# =============================================================================
# Recréation de l'analyse FEA de A. Harral (varmintal.com/a22lr.htm)
# « 22 Long Rifle Barrel Tuner Analysis -- FEA Dynamic Analysis »
#
# On reproduit la MÉTHODE de Harral dans notre propre cadre MEF :
#   1. Poutre d'Euler-Bernoulli à SECTION VARIABLE = le canon reverse-taper
#      d'« Esten » (contour donné par la page varmintal).
#   2. Réponse transitoire (Newmark-β) au moment de recul + poids du projectile.
#   3. À l'instant de sortie t_b (propre à chaque vitesse), on lit :
#        - la projection de bouche  proj = θ(L,t_b)·D_cible   [pouces]
#        - la vitesse verticale de bouche  vvel = ẏ(L,t_b)     [in/s]
#   4. Point d'impact à 50 yd :  POI = proj + vvel·TOF − chute
#      dispersion verticale = |POI(1035 fps) − POI(1075 fps)|
#   5. On balaie les masses de tuner de Harral (0, 4.9, 8.6, 16.0 oz).
#
# GRAVITÉ : le poids propre du canon et celui du tuner sont appliqués (charges
# cohérentes), et le transitoire DÉMARRE à l'équilibre statique. La colonne
# `proj` de Harral est en effet une flèche STATIQUE, pas une vibration — il
# l'écrit : « Gravity is applied and the stock deforms as well as the barrel
# sag. That is why the muzzle starts out pointing approximately 1.35 inches
# below zero. » Cette colonne est désormais reproduite SANS aucun calage ; le
# reliquat (≈ 0.2-0.4 in) est la flèche de la CROSSE, que Harral inclut et que
# nous ne modélisons pas encore.
#
# RÉSERVE — LA COLONNE `spread` N'EST TOUJOURS PAS REPRODUITE. Avec la culasse
# encastrée, la carabine ne peut pas tourner en corps rigide au recul : il ne
# reste que la flexion élastique, et la compensation positive n'apparaît pas.
# Toutes les configs convergent vers ≈ 0.17 in — la dispersion NON compensée
# (cf. ligne « bouche rigide ») — au lieu de la structure de Harral (0.091 nu,
# croissant à 0.183). Ce n'est pas un défaut d'échelle : le système est
# linéaire, aucun réglage de l'amplitude ne peut créer la compensation. Il
# manque la ROTATION DE CORPS RIGIDE du fusil sur ses appuis. Harral :
# « The stock contacts the simulated sandbag rests with zero friction in the
# calculation. The rests are fixed in space. » — soit un APPUI UNILATÉRAL et
# NON un encastrement : réaction verticale seule, ni effort horizontal, ni
# moment ; l'arme est libre de reculer et de tourner, mais reste posée sur des
# appuis fixes dans l'espace (c'est cette réaction qui permet la flèche).
# Chantier suivant : remplacer l'encastrement par deux appuis sans frottement
# (avant-bras + talon) + d.d.l. de corps rigide, puis flexibilité de la crosse.
#
# PIÈGE : chez Harral, tableau et graphiques n'ont pas la même référence —
# « The sag due to gravity was in the calculations, but subtracted out for the
# chart so the curves can more easily be compared. » La flèche est retirée des
# COURBES, pas de la colonne `proj` du tableau 1 (qui part bien de −1.32 in).
#
# Table 1 de Harral (référence, pouces) :
#   Config              proj(1035)  proj(1075)  vvel(1035) vvel(1075)  spread
#   Reverse Taper nu     -1.3183     -1.3869     1.3303     1.2790     0.0910
#   + 4.9 oz             -1.7178     -1.7434     0.6895     0.5327     0.1218
#   + 8.6 oz             -1.9640     -1.9682     0.5265     0.4377     0.1554
#   + 16.0 oz            -2.4183     -2.4022     0.4378     0.4076     0.1832
#   Bouche rigide         0.0         0.0        0.0        0.0        0.1737
#
# Usage :  julia harral_a22lr.jl
# =============================================================================

using LinearAlgebra
using Printf

const PLOTS_AVAILABLE = try
    using Plots; true
catch
    println("→ Plots.jl absent : les tracés seront ignorés (tableau seul).")
    false
end

# -----------------------------------------------------------------------------
# 0. CONVERSIONS
# -----------------------------------------------------------------------------
const IN   = 0.0254            # 1 pouce → m
const OZ   = 0.0283495         # 1 once  → kg
const M2IN = 1 / IN            # m → pouces

# -----------------------------------------------------------------------------
# 1. GÉOMÉTRIE : canon reverse-taper d'Esten (contour varmintal, en pouces)
#   - cylindre de culasse : 1.1" dia sur 1"
#   - raccordement descendant vers 0.750" (rayon, ici linéaire sur ~1")
#   - taper MONTANT 0.750" → 0.915"
#   - cylindre de bouche : 0.915" sur les 2 derniers pouces
#   Alésage .224" (.22 LR). Longueur totale 24.75".
# -----------------------------------------------------------------------------
const L_IN    = 24.75
const L       = L_IN * IN
const D_bore  = 0.224 * IN

# Diamètre extérieur (m) en fonction de la position axiale x (m)
function D_out_of_x(x)
    xi = x * M2IN               # position en pouces
    if xi <= 1.0
        return 1.100 * IN                                   # cyl. culasse
    elseif xi <= 2.0
        return (1.100 + (0.750 - 1.100) * (xi - 1.0)) * IN  # raccord ↓
    elseif xi <= 22.75
        return (0.750 + (0.915 - 0.750) * (xi - 2.0) / (22.75 - 2.0)) * IN  # taper ↑
    else
        return 0.915 * IN                                   # cyl. bouche
    end
end

# -----------------------------------------------------------------------------
# 2. MATÉRIAU ET MAILLAGE
# -----------------------------------------------------------------------------
const E        = 200e9         # module d'Young acier (Pa)
const ρ_steel  = 7850.0        # masse volumique acier (kg/m³)

const N_elements = 48
const N_nodes    = N_elements + 1
const ndof       = 2 * N_nodes
const L_e        = L / N_elements

const A_bore = π/4 * D_bore^2

# Propriétés par élément (évaluées au milieu de l'élément)
function element_section(e)
    x_mid = (e - 0.5) * L_e
    D = D_out_of_x(x_mid)
    A = π/4  * (D^2 - D_bore^2)
    I = π/64 * (D^4 - D_bore^4)
    return A, I
end

# -----------------------------------------------------------------------------
# 3. MATRICES ÉLÉMENTAIRES (Euler-Bernoulli, Hermite cubique)
# -----------------------------------------------------------------------------
function element_matrices(L_e, EI, ρA)
    Ke = (EI / L_e^3) * [
        12.0     6*L_e       -12.0    6*L_e   ;
        6*L_e    4*L_e^2     -6*L_e   2*L_e^2 ;
       -12.0    -6*L_e        12.0   -6*L_e   ;
        6*L_e    2*L_e^2     -6*L_e   4*L_e^2 ]
    Me = (ρA * L_e / 420.0) * [
        156.0    22*L_e       54.0   -13*L_e   ;
        22*L_e   4*L_e^2      13*L_e  -3*L_e^2 ;
        54.0     13*L_e      156.0   -22*L_e   ;
       -13*L_e  -3*L_e^2     -22*L_e   4*L_e^2 ]
    return Ke, Me
end

# -----------------------------------------------------------------------------
# 4. INERTIE DU TUNER (cylindre annulaire acier serré sur la bouche)
#   OD ≈ 1.4", ID ≈ 0.915" (dia bouche) ; la longueur découle de la masse.
#   J transverse ≈ m·(3(Ro²+Ri²)+ℓ²)/12 (CoM supposé au nœud de bouche).
#
# ARCHITECTURE : BAGUE COMPACTE, et c'est un choix différent de celui de
# simulation.jl — divergence assumée, à ne pas « corriger » par alignement.
# Ici J est DÉRIVÉ de la géométrie, donc le rayon de giration équivalent
# k = √(J/m) croît avec la masse : 1,1 cm à 50 g, 1,7 cm à 200 g, 2,8 cm à
# 400 g. simulation.jl modélise au contraire un ENSEMBLE TUBE (Starik/Centra)
# avec k = 5,0 cm fixe, seule architecture autorisant son porte-à-faux de
# ~10 cm. À 200 g : J = 5,6e-5 ici contre 5,0e-4 là-bas, facteur ~9.
#
# CONSÉQUENCE. Les θ̇ et les optima des deux familles ne sont PAS comparables
# entre eux — ils décrivent deux produits réels distincts, pas deux estimations
# du même. Seules les grandeurs cinématiques (τ_v, taux de compensation requis,
# cf. kolbe_validation.jl) se recoupent, précisément parce qu'elles ne dépendent
# d'aucune inertie de tuner.
# -----------------------------------------------------------------------------
function tuner_inertia(m)
    m <= 0 && return 0.0
    Ro = 0.5 * 1.40 * IN
    Ri = 0.5 * 0.915 * IN
    A  = π * (Ro^2 - Ri^2)
    ℓ  = m / (ρ_steel * A)
    return m * (3*(Ro^2 + Ri^2) + ℓ^2) / 12
end

# -----------------------------------------------------------------------------
# 5. ASSEMBLAGE + ENCASTREMENT CULASSE + TUNER À LA BOUCHE
# -----------------------------------------------------------------------------
function build_system(m_tuner)
    K = zeros(ndof, ndof)
    M = zeros(ndof, ndof)
    for e in 1:N_elements
        A, I = element_section(e)
        Ke, Me = element_matrices(L_e, E * I, ρ_steel * A)
        idx = (2*e - 1):(2*e + 2)
        @views K[idx, idx] .+= Ke
        @views M[idx, idx] .+= Me
    end
    active = 3:ndof                       # encastrement : y₁ = θ₁ = 0
    Ka = K[active, active]
    Ma = copy(M[active, active])
    Ma[end-1, end-1] += m_tuner           # masse ponctuelle à la bouche
    Ma[end,   end]   += tuner_inertia(m_tuner)
    return Ka, Ma
end

# -----------------------------------------------------------------------------
# 6. MODES PROPRES
# -----------------------------------------------------------------------------
# `want_modes` : renvoie en plus les déformées, NORMALISÉES EN MASSE
# (Φᵀ M Φ = I), sans quoi la projection modale d'une réponse n'a pas de sens.
# Sert à identifier QUEL mode porte θ̇ à l'instant de sortie de la balle.
function modal_analysis(Ka, Ma; n_modes = 5, want_modes = false)
    decomp = eigen(Ka, Ma)
    λ = real.(decomp.values); Φ = real.(decomp.vectors)
    keep = λ .> 1e-3
    λ = λ[keep]; Φ = Φ[:, keep]
    order = sortperm(λ); λ = λ[order]; Φ = Φ[:, order]
    ω = sqrt.(λ)
    nm = min(n_modes, length(ω))
    if want_modes
        Φn = Φ[:, 1:nm]
        for k in 1:nm
            Φn[:, k] ./= sqrt(Φn[:, k]' * Ma * Φn[:, k])
        end
        return ω[1:nm] ./ (2π), ω[1:nm], Φn
    end
    return ω[1:nm] ./ (2π), ω[1:nm]
end

# -----------------------------------------------------------------------------
# 7. BALISTIQUE INTERNE : profil « burnout » — accélération constante sur une
#   fraction φ du canon, puis coast à v_muzzle.
#
#   POURQUOI CE PROFIL. Le rôle de la balistique interne ici est de fixer le
#   temps de sortie t_b ET sa sensibilité à la vitesse τ_v = −∂t_b/∂v₀, qui
#   place la sortie dans le cycle vibratoire. Kolbe mesure τ_v = 8,8 µs/(m/s).
#   L'ancien profil (lag exponentiel v = v_m(1−e^{−t/τ})) est « front-loaded » :
#   toute l'accélération à la culasse, puis coast → τ_v plafonne à L/v² ≈ 6,5 µs,
#   quel que soit τ. Il ne pouvait PAS atteindre 8,8.
#
#   Le profil burnout donne, analytiquement :
#     t_b  = (1+φ)·L / v_muzzle
#     τ_v  = (1+φ)·L / v_muzzle²            (entre L/v² à φ=0 et 2L/v² à φ=1)
#   φ = 0,35 reproduit τ_v = 8,8 µs/(m/s) sur le canon 26" de Kolbe (v=318).
#   Physiquement : la balle accélère sur le premier tiers du canon, puis coaste
#   sous frottement — caricature raisonnable d'une .22 LR.
# -----------------------------------------------------------------------------
const PHI_BURN = 0.35        # fraction de canon en phase d'accélération

exit_time(v_muzzle; φ = PHI_BURN) = (1 + φ) * L / v_muzzle

function projectile_pos(v_muzzle; φ = PHI_BURN)
    x_bo  = φ * L                       # position de fin d'accélération
    a     = v_muzzle^2 / (2 * x_bo)     # accélération constante
    t_acc = v_muzzle / a                # = 2φL/v_muzzle
    return function (t)
        t <= 0     && return 0.0
        t < t_acc  && return 0.5 * a * t^2
        return x_bo + v_muzzle * (t - t_acc)
    end
end

# -----------------------------------------------------------------------------
# 8. EXCITATION : moment de recul (culasse) + poids du projectile (charge mobile)
# -----------------------------------------------------------------------------
const p_max   = 200e6          # pic de pression .22 LR (Pa)
const t_peak  = 0.30e-3        # pic ~0.3 ms (impulsion brève du rimfire)
const α_press = 2.0
const m_p     = 2.6e-3         # 40 gr
const g_ms    = 9.8044         # g en m/s² (= 386 in/s²)

chamber_pressure(t) = (t <= 0 || t >= 6*t_peak) ? 0.0 :
    p_max * (t/t_peak)^α_press * exp(α_press*(1 - t/t_peak))

function consistent_point_load(F, x, ndof_a)
    (x <= 0 || x >= L) && return zeros(ndof_a)
    e = clamp(Int(floor(x / L_e)) + 1, 1, N_elements)
    ξ = (x - (e-1)*L_e) / L_e
    Ns = (1 - 3ξ^2 + 2ξ^3, L_e*(ξ - 2ξ^2 + ξ^3), 3ξ^2 - 2ξ^3, L_e*(-ξ^2 + ξ^3))
    Fa = zeros(ndof_a)
    for (k, idx_g) in enumerate((2*e-1, 2*e, 2*e+1, 2*e+2))
        ia = idx_g - 2
        1 <= ia <= ndof_a && (Fa[ia] += F * Ns[k])
    end
    return Fa
end

# Précharge de gravité : poids propre du canon (charge répartie cohérente sur
# chaque élément) + poids du tuner (charge ponctuelle au nœud de bouche).
# Constante dans le temps ⇒ calculée une fois, réutilisée à chaque pas.
function gravity_force(m_tuner, ndof_a)
    F = zeros(ndof_a)
    for e in 1:N_elements
        A, _ = element_section(e)
        q  = -ρ_steel * A * g_ms                  # N/m, vers le bas
        Fe = q * L_e * [0.5, L_e/12, 0.5, -L_e/12]
        for (k, idx_g) in enumerate((2*e-1, 2*e, 2*e+1, 2*e+2))
            ia = idx_g - 2                        # décalage d'encastrement
            1 <= ia <= ndof_a && (F[ia] += Fe[k])
        end
    end
    F[end-1] += -m_tuner * g_ms                   # poids du tuner à la bouche
    return F
end

function force_vector(t, x_p, ndof_a, h_offset, F_grav)
    F = copy(F_grav)                                  # gravité, en permanence
    F[2] += chamber_pressure(t) * A_bore * h_offset   # moment de recul → θ₂
    xp = x_p(t)
    0 < xp < L && (F .+= consistent_point_load(-m_p * g_ms, xp, ndof_a))
    return F
end

# -----------------------------------------------------------------------------
# 9. AMORTISSEMENT DE RAYLEIGH
# -----------------------------------------------------------------------------
function rayleigh_damping(Ma, Ka, ω1, ω2, ζ1, ζ2)
    A = [1/(2ω1) ω1/2; 1/(2ω2) ω2/2]
    α, β = A \ [ζ1; ζ2]
    return α*Ma + β*Ka
end

# -----------------------------------------------------------------------------
# 10. NEWMARK-β (γ=1/2, β=1/4)
# -----------------------------------------------------------------------------
function newmark_solve(Ma, Ca, Ka, F_of_t, t_end, Δt; U0 = nothing)
    γ, β = 0.5, 0.25
    n = size(Ma, 1)
    ts = collect(0:Δt:t_end); Nt = length(ts)
    U = zeros(n, Nt); V = zeros(n, Nt); Ac = zeros(n, Nt)
    U0 !== nothing && (U[:, 1] = U0)      # départ à l'équilibre statique
    Ac[:, 1] = Ma \ (F_of_t(ts[1]) - Ca*V[:,1] - Ka*U[:,1])
    Kf = factorize(Ma + γ*Δt*Ca + β*Δt^2*Ka)
    for i in 1:Nt-1
        up = U[:,i] + Δt*V[:,i] + Δt^2*(0.5-β)*Ac[:,i]
        vp = V[:,i] + Δt*(1-γ)*Ac[:,i]
        Ac[:, i+1] = Kf \ (F_of_t(ts[i+1]) - Ca*vp - Ka*up)
        U[:, i+1] = up + β*Δt^2*Ac[:,i+1]
        V[:, i+1] = vp + γ*Δt*Ac[:,i+1]
    end
    return ts, U, V
end

# -----------------------------------------------------------------------------
# 11. UN TIR : renvoie l'historique bouche (y, θ, ẏ) et les valeurs à t_b
# -----------------------------------------------------------------------------
function simulate_shot(m_tuner, v_muzzle; h_offset, ζ1 = 0.01, ζ2 = 0.06,
                       Δt = 2e-6, t_end = 8e-3)
    Ka, Ma = build_system(m_tuner)
    _, ωs  = modal_analysis(Ka, Ma; n_modes = 3)
    Ca     = rayleigh_damping(Ma, Ka, ωs[1], ωs[2], ζ1, ζ2)
    ndof_a = size(Ka, 1)
    x_p    = projectile_pos(v_muzzle)

    F_grav = gravity_force(m_tuner, ndof_a)
    U_stat = Ka \ F_grav                  # flèche statique sous gravité
    ts, U, V = newmark_solve(Ma, Ca, Ka,
                             t -> force_vector(t, x_p, ndof_a, h_offset, F_grav),
                             t_end, Δt; U0 = U_stat)

    y_L  = U[end-1, :]      # déflexion transverse bouche (m)
    θ_L  = U[end,   :]      # angle bouche (rad)
    ẏ_L  = V[end-1, :]      # vitesse transverse bouche (m/s)

    # Composante dynamique = écart à l'équilibre statique (gravité retirée).
    θ_dyn = θ_L .- U_stat[end]

    t_b   = exit_time(v_muzzle)
    ib    = argmin(abs.(ts .- t_b))
    return (ts=ts, y_L=y_L, θ_L=θ_L, ẏ_L=ẏ_L, θ_dyn=θ_dyn,
            y_stat=U_stat[end-1], θ_stat=U_stat[end],
            t_b=ts[ib], ib=ib,
            θ_tb=θ_L[ib], ẏ_tb=ẏ_L[ib], f1=ωs[1]/(2π))
end

# -----------------------------------------------------------------------------
# 12. POINT D'IMPACT (méthode de Harral)
#   proj  = θ(L,t_b)·D_cible          [pouces]   (D_cible = 1800 in = 50 yd)
#   vvel  = ẏ(L,t_b)                  [in/s]
#   chute et TOF : valeurs balistiques de Harral (BC 0.128, 40 gr, .224)
#     1035 fps → chute 4.345 in ;  1075 fps → chute 4.172 in
#   POI   = proj + vvel·TOF − chute
# -----------------------------------------------------------------------------
const D_TARGET_IN = 50 * 36.0                 # 50 yd en pouces
const DROP = Dict(1035 => 4.345, 1075 => 4.172)   # chute balistique (Harral), in
TOF(drop_in) = sqrt(2 * drop_in / 386.0)          # s  (g = 386 in/s²)

function point_of_impact(res, vfps)
    proj = res.θ_tb * D_TARGET_IN              # in
    vvel = res.ẏ_tb * M2IN                     # m/s → in/s
    drop = DROP[vfps]
    return proj + vvel * TOF(drop) - drop, proj, vvel
end

# -----------------------------------------------------------------------------
# 13. CALIBRATION de l'amplitude d'excitation (h_offset)
#   Système linéaire en h_offset ⇒ un tir suffit, puis on rééchelonne pour que
#   la projection de bouche du CANON NU à 1035 fps atteigne la cible de Harral.
# -----------------------------------------------------------------------------
const V_LO, V_HI = 315.47, 327.66              # 1035 et 1075 fps en m/s

# Calibration de l'amplitude d'excitation — PARAMÈTRE NON CONTRAINT.
#
# Historique : on calait le pic de |projection| du canon nu sur 1.5 in, en
# référence à la projection de bouche publiée par Harral (≈ −1.32 in). C'était
# une ERREUR : ce nombre est une FLÈCHE STATIQUE sous gravité (« Gravity is
# applied and the stock deforms as well as the barrel sag. That is why the
# muzzle starts out pointing approximately 1.35 inches below zero. »), pas une
# amplitude de vibration. Comparer l'une à l'autre gonflait l'excitation d'un
# ordre de grandeur.
#
# La flèche statique est désormais produite par la gravité (§ gravity_force),
# sans aucun calage. Il ne reste donc RIEN dans le tableau de Harral sur quoi
# caler l'amplitude DYNAMIQUE : h_offset est un paramètre libre. La valeur
# ci-dessous est un ordre de grandeur assumé, pas une mesure — le bras de
# levier du moment de recul restera mal défini tant que la culasse sera
# encastrée (il dépend de l'appui réel de la crosse : cf. rotation de corps
# rigide, chantier à venir).
const PEAK_DYN_PROJ = 0.15                # pic |projection| DYNAMIQUE visé (in)

function calibrate_h_offset()
    h_ref = 1e-3
    res = simulate_shot(0.0, V_LO; h_offset = h_ref)
    peak_ref = maximum(abs.(res.θ_dyn)) * D_TARGET_IN   # dynamique seule
    h = h_ref * (PEAK_DYN_PROJ / peak_ref)
    @printf("Amplitude dynamique (paramètre libre) : pic visé = %.2f in → h_offset = %.2f mm\n",
            PEAK_DYN_PROJ, h*1e3)
    return h
end

# =============================================================================
# EXÉCUTION
# =============================================================================
function main()
println("="^78)
println(" Recréation de l'analyse FEA de Harral — canon reverse-taper .22 LR d'Esten")
println("="^78)

# Fréquences propres du canon nu (repère)
let (fs, _) = modal_analysis(build_system(0.0)...; n_modes = 3)
    @printf("Canon nu — fréquences propres : f₁=%.1f Hz  f₂=%.1f Hz  f₃=%.1f Hz\n",
            fs[1], fs[2], fs[3])
end
@printf("Temps de sortie : t_b(1035)=%.3f ms   t_b(1075)=%.3f ms   (Δ=%.3f ms)\n",
        exit_time(V_LO)*1e3, exit_time(V_HI)*1e3, (exit_time(V_LO)-exit_time(V_HI))*1e3)
println()

h_cal = calibrate_h_offset()
println()

# Balayage des masses de tuner de Harral
configs = [("Reverse Taper nu", 0.0),
           ("+ 4.9 oz",  4.9*OZ),
           ("+ 8.6 oz",  8.6*OZ),
           ("+ 16.0 oz", 16.0*OZ)]

harral_spread = Dict("Reverse Taper nu"=>0.0910, "+ 4.9 oz"=>0.1218,
                     "+ 8.6 oz"=>0.1554, "+ 16.0 oz"=>0.1832)

println("Tableau 1 recréé (dispersion verticale à 50 yd, pouces)")
println("-"^78)
@printf("%-18s | %7s | %8s %8s | %7s %7s | %8s | %8s\n",
        "Config", "f₁(Hz)", "proj1035", "proj1075", "vv1035", "vv1075",
        "spread", "Harral")
println("-"^78)

results = Dict{String,Any}()
for (name, m) in configs
    r_lo = simulate_shot(m, V_LO; h_offset = h_cal)
    r_hi = simulate_shot(m, V_HI; h_offset = h_cal)
    poi_lo, proj_lo, vv_lo = point_of_impact(r_lo, 1035)
    poi_hi, proj_hi, vv_hi = point_of_impact(r_hi, 1075)
    spread = abs(poi_lo - poi_hi)
    results[name] = (r_lo=r_lo, r_hi=r_hi, spread=spread)
    @printf("%-18s | %7.1f | %+8.4f %+8.4f | %+7.4f %+7.4f | %8.4f | %8.4f\n",
            name, r_lo.f1, proj_lo, proj_hi, vv_lo, vv_hi, spread, harral_spread[name])
end
# Ligne bouche rigide (purement balistique)
let spread_rigid = abs((-DROP[1035]) - (-DROP[1075]))
    @printf("%-18s | %7s | %+8.4f %+8.4f | %+7.4f %+7.4f | %8.4f | %8.4f\n",
            "Bouche rigide", "—", 0.0, 0.0, 0.0, 0.0, spread_rigid, 0.1737)
end
println("-"^78)
println("proj = projection de bouche (in) ; vv = vitesse verticale de bouche (in/s).")
println()
println("proj : flèche statique sous gravité, reproduite SANS calage. La tendance")
println("       de Harral (−1.32 → −2.42) est retrouvée ; le reliquat = la crosse,")
println("       qu'il modélise et que nous ignorons.")
println("spread : NON reproduit. Culasse encastrée ⇒ pas de rotation de corps rigide")
println("       ⇒ pas de compensation : tout converge vers la dispersion non")
println("       compensée (≈ ligne « bouche rigide »). Voir l'en-tête du fichier.")

# -----------------------------------------------------------------------------
# Tracés : contour du canon + courbes de projection de bouche
# -----------------------------------------------------------------------------
if PLOTS_AVAILABLE
    # (a) Contour du canon reverse-taper
    xs = range(0, L, length = 400)
    Ds = [D_out_of_x(x) * M2IN for x in xs]        # pouces
    p_contour = plot(xs .* M2IN, Ds ./ 2, color=:steelblue, lw=2, legend=false,
                     xlabel="Position depuis la culasse (in)", ylabel="Rayon (in)",
                     title="Contour du canon reverse-taper d'Esten")
    plot!(p_contour, xs .* M2IN, -Ds ./ 2, color=:steelblue, lw=2)
    plot!(p_contour, xs .* M2IN, fill(0.224/2, length(xs)), color=:gray, ls=:dot)
    plot!(p_contour, xs .* M2IN, fill(-0.224/2, length(xs)), color=:gray, ls=:dot)
    savefig(p_contour, "harral_contour.png")
    println("\n→ harral_contour.png")

    # (b) Courbes de projection de bouche vs temps, FLÈCHE STATIQUE RETIRÉE.
    #     Même convention que Harral pour ses graphiques : « The sag due to
    #     gravity was in the calculations, but subtracted out for the chart so
    #     the curves can more easily be compared. » Sans cela, les courbes sont
    #     décalées de −1.13 à −2.02 in et l'oscillation devient illisible.
    #     ATTENTION : sa colonne `proj` du tableau 1, elle, GARDE la flèche.
    p_proj = plot(xlabel="Temps (ms)", ylabel="Projection de bouche à 50 yd (in)",
                  title="Projection de bouche — flèche statique retirée",
                  legend=:bottomright)   # en haut à droite, elle couvrait les courbes
    palette = [:black, :darkorange, :seagreen, :crimson]
    for (i, (name, m)) in enumerate(configs)
        r = results[name].r_lo
        mask = r.ts .<= 4e-3
        plot!(p_proj, r.ts[mask].*1e3, r.θ_dyn[mask].*D_TARGET_IN,
              color=palette[i], lw=2, label=name)
        scatter!(p_proj, [r.t_b*1e3], [(r.θ_tb - r.θ_stat)*D_TARGET_IN],
                 color=palette[i], ms=5, label="")
    end
    vline!(p_proj, [exit_time(V_LO)*1e3], color=:gray, ls=:dash, label="t_b(1035)")
    hline!(p_proj, [0.0], color=:black, ls=:dot, alpha=0.4, label="")
    savefig(p_proj, "harral_projection.png")
    println("→ harral_projection.png")
end

println("="^78)
end # main

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
