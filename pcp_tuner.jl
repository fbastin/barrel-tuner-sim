# =============================================================================
# ACCORD D'UN TUNER SUR ARME À AIR PRÉCOMPRIMÉE (PCP)
#
# Un PCP est quasi sans recul : le moment de recul qui excite le canon d'une
# arme à feu y est négligeable (cf. la section « Cas du PCP » du document). La
# vibration naît d'AUTRES sources — surtout la frappe marteau/soupape — qu'on
# ne sait pas mesurer sans banc instrumenté. Et la cible d'accord n'est plus la
# compensation positive (θ̇ = +6 MOA/ms) mais la BOUCHE STATIONNAIRE : θ̇ ≈ 0 à
# la sortie, où une gigue sur l'instant de sortie ne change plus l'angle de
# lancement.
#
# IDÉE-CLÉ QUI REND CE CALCUL POSSIBLE. Le système est linéaire :
#     θ̇(t) = (amplitude d'excitation) × (réponse unitaire).
# Le ZÉRO de θ̇ — donc la POSITION des sweet spots dans le balayage du tuner —
# ne dépend que de la réponse unitaire (fréquences modales + temps de sortie
# t_b), PAS de l'amplitude de l'excitation. On peut donc prédire À QUEL réglage
# la bouche devient stationnaire sans connaître l'excitation. Ce qui reste hors
# de portée : la PROFONDEUR du bénéfice (de combien le groupement se resserre),
# qui dépend de l'amplitude — donc d'une mesure.
#
# CE QUE LE CALCUL MONTRE (et corrige une idée reçue). Le whip du canon
# (~250-700 Hz) oscille vite DANS LE TEMPS, mais sa fréquence ne se décale que
# de ~20 % quand on charge le tuner de 100 g : la phase ω·t_b ne balaie qu'une
# fraction de cycle. Il y a donc, sur la plage réaliste d'un tuner (~0-80 g),
# UN sweet spot, pas un peigne rapproché — ce qui colle à la pratique (on
# scanne la course du tuner par tours pour trouver LE bon réglage). Les sweet
# spots se répètent, mais très espacés (~100 g), hors de portée pratique.
#
# LE RÉGLAGE RÉEL = LA POSITION. En pratique on FIXE la masse du tuner (le poids
# qu'on monte) et on accorde en le VISSANT en porte-à-faux devant la bouche —
# réglage continu et fin, compté en tours de filetage. Ce script centre donc la
# POSITION (§ main §1) : à masse fixe, θ̇(t_b) croise zéro à un porte-à-faux
# précis — LE sweet spot, trouvé en scannant la course, comme la procédure
# d'accord réelle. La masse, elle, déplace ce sweet spot le long de la course
# (§ main §2, la « courbe d'accord ») : masse et position forment UN espace
# d'accord, le poids fixe la courbe, la position accorde dessus.
#
# RÉSERVE : la position exacte des sweet spots suppose qu'UNE source
# d'excitation domine (ici la frappe marteau, modélisée par une impulsion
# brève à la culasse). Si plusieurs sources de poids comparable se mélangent,
# les positions se décalent — mais l'existence et l'espacement des sweet spots
# restent robustes, car ils procèdent du décalage de phase induit par le tuner.
#
# Usage :  julia pcp_tuner.jl
# =============================================================================

using LinearAlgebra
using Printf

# -----------------------------------------------------------------------------
# CONFIGURATION — bloc « munition + canon » (facile à changer)
# -----------------------------------------------------------------------------
const AMMO = (
    name   = ".22 PCP (plomb 16 gr, ~850 fps)",
    m_p    = 1.04e-3,     # 16 gr
    D_bore = 0.0055,      # .22 cal (5,5 mm)
    v_exit = 260.0,       # ~850 fps
    φ_burn = 0.35,        # fraction de canon en accélération (profil burnout)
)

const BARREL = (
    L     = 0.55,         # longueur (m)
    D_out = 0.016,        # diamètre extérieur (canon airgun mince, 16 mm)
    E     = 200e9,        # module d'Young acier (Pa)
    ρ     = 7850.0,       # masse volumique acier (kg/m³)
    N     = 24,           # éléments finis
)

const TUNER_MAX  = 0.200  # masse de tuner balayée jusqu'à 200 g (au-delà du réaliste,
                          # pour révéler la récurrence — le 2e sweet spot)
const TUNER_STEP = 0.001  # pas de 1 g
const ζ_STEEL    = 0.005  # amortissement matériau (acier)

# RÉGLAGE RÉEL = LA POSITION. En pratique la masse du tuner est FIXE (on possède
# un poids), et l'on accorde en le VISSANT plus ou moins loin en porte-à-faux
# devant la bouche — réglage continu et fin, repéré par tours de filetage.
const TUNER_MASS_FIX = 0.070  # masse fixe du tuner (70 g) — le poids qu'on visse
const OVERHANG_MAX   = 0.100  # porte-à-faux balayé jusqu'à 100 mm
const OVERHANG_STEP  = 0.0005 # pas de 0,5 mm
const THREAD_PITCH   = 0.0010 # pas du filetage : 1 mm/tour (pour exprimer en « tours »)

const MOA_per_rad = 180 * 60 / π
moa_per_ms(rad_s) = rad_s * MOA_per_rad * 1e-3

# Grandeurs dérivées
const L_e    = BARREL.L / BARREL.N
const A_sec  = π/4  * (BARREL.D_out^2 - AMMO.D_bore^2)
const I_sec  = π/64 * (BARREL.D_out^4 - AMMO.D_bore^4)
const EI     = BARREL.E * I_sec
const ρA     = BARREL.ρ * A_sec
const ndof   = 2 * (BARREL.N + 1)

# -----------------------------------------------------------------------------
# FEM : poutre d'Euler–Bernoulli (Hermite cubique), culasse encastrée + tuner
# -----------------------------------------------------------------------------
function element_matrices()
    Ke = (EI / L_e^3) * [
        12.0     6*L_e     -12.0    6*L_e   ;
        6*L_e    4*L_e^2   -6*L_e   2*L_e^2 ;
       -12.0    -6*L_e      12.0   -6*L_e   ;
        6*L_e    2*L_e^2   -6*L_e   4*L_e^2 ]
    Me = (ρA * L_e / 420.0) * [
        156.0    22*L_e     54.0   -13*L_e   ;
        22*L_e   4*L_e^2    13*L_e  -3*L_e^2 ;
        54.0     13*L_e    156.0   -22*L_e   ;
       -13*L_e  -3*L_e^2   -22*L_e   4*L_e^2 ]
    return Ke, Me
end

# Tuner = masse m en PORTE-À-FAUX à la distance d DEVANT la bouche (stem rigide).
# Le déplacement vertical de son centre de masse est y_bouche + d·θ_bouche, d'où
# une matrice de masse ajoutée au nœud de bouche couplant translation/rotation :
#   M_add = [ m      m·d      ]     (J_cm = inertie propre du tuner)
#           [ m·d   m·d²+J_cm ]
# À d = 0, on retrouve la masse ponctuelle à la bouche.
function build_system(m_tuner; d_overhang = 0.0)
    Ke, Me = element_matrices()
    K = zeros(ndof, ndof); M = zeros(ndof, ndof)
    for e in 1:BARREL.N
        idx = (2*e - 1):(2*e + 2)
        @views K[idx, idx] .+= Ke
        @views M[idx, idx] .+= Me
    end
    active = 3:ndof                       # encastrement culasse : y₁ = θ₁ = 0
    Ka = K[active, active]; Ma = copy(M[active, active])
    J_cm = m_tuner * 2.0e-3               # inertie propre (rayon giration ~4,5 cm)
    d = d_overhang
    Ma[end-1, end-1] += m_tuner
    Ma[end-1, end  ] += m_tuner * d
    Ma[end,   end-1] += m_tuner * d
    Ma[end,   end  ] += m_tuner * d^2 + J_cm
    return Ka, Ma
end

function modal(Ka, Ma; n = 4)
    d = eigen(Ka, Ma)
    λ = real.(d.values); keep = λ .> 1e-3
    ω = sqrt.(sort(λ[keep]))
    return ω[1:min(n, length(ω))]
end

function rayleigh(Ma, Ka, ω1, ω2, ζ)
    A = [1/(2ω1) ω1/2; 1/(2ω2) ω2/2]
    α, β = A \ [ζ, ζ]
    return α*Ma + β*Ka
end

# -----------------------------------------------------------------------------
# Balistique interne (profil burnout) et EXCITATION placeholder (marteau)
# -----------------------------------------------------------------------------
exit_time() = (1 + AMMO.φ_burn) * BARREL.L / AMMO.v_exit

# Impulsion de frappe marteau/soupape : moment bref à la culasse. Amplitude
# ARBITRAIRE (ne change pas la position des sweet spots — voir en-tête).
const HAMMER_TPK = 0.15e-3                # choc mécanique bref (~0,15 ms)
hammer_moment(t) = (t <= 0 || t >= 6*HAMMER_TPK) ? 0.0 :
    (t/HAMMER_TPK)^2 * exp(2*(1 - t/HAMMER_TPK))

function newmark(Ma, Ca, Ka, F_of_t, t_end; Δt = 2e-6)
    n = size(Ma,1); ts = collect(0:Δt:t_end); Nt = length(ts)
    U = zeros(n,Nt); V = zeros(n,Nt); Ac = zeros(n,Nt)
    Ac[:,1] = Ma \ (F_of_t(ts[1]) - Ca*V[:,1] - Ka*U[:,1])
    Kf = factorize(Ma + 0.5*Δt*Ca + 0.25*Δt^2*Ka)
    for i in 1:Nt-1
        up = U[:,i] + Δt*V[:,i] + Δt^2*0.25*Ac[:,i]
        vp = V[:,i] + Δt*0.5*Ac[:,i]
        Ac[:,i+1] = Kf \ (F_of_t(ts[i+1]) - Ca*vp - Ka*up)
        U[:,i+1] = up + 0.25*Δt^2*Ac[:,i+1]
        V[:,i+1] = vp + 0.5*Δt*Ac[:,i+1]
    end
    return ts, U, V
end

# θ̇(L, t_b) en MOA/ms pour un tuner (masse m_tuner, porte-à-faux d_overhang)
function theta_dot_tb(m_tuner; d_overhang = 0.0)
    Ka, Ma = build_system(m_tuner; d_overhang)
    ωs = modal(Ka, Ma; n = 4)
    Ca = rayleigh(Ma, Ka, ωs[1], ωs[2], ζ_STEEL)   # modes de flexion du canon
    ndof_a = size(Ka, 1)
    F_of_t = t -> (F = zeros(ndof_a); F[2] += hammer_moment(t); F)  # moment à la culasse
    t_b = exit_time()
    ts, U, V = newmark(Ma, Ca, Ka, F_of_t, 1.5*t_b)
    ib = argmin(abs.(ts .- t_b))
    return moa_per_ms(V[end, ib]), ωs[1]/(2π)      # θ̇(t_b), f₁
end

# Repère les passages par zéro (sweet spots) dans un balayage x → θ̇
function sweet_spots(xs, θdots)
    s = Float64[]
    for i in 1:length(xs)-1
        if θdots[i] == 0 || sign(θdots[i]) != sign(θdots[i+1])
            push!(s, xs[i] - θdots[i]*(xs[i+1]-xs[i])/(θdots[i+1]-θdots[i]))
        end
    end
    return s
end

# =============================================================================
function main()
    println("="^74)
    println(" Accord d'un tuner sur PCP — recherche des sweet spots (bouche stationnaire)")
    println("="^74)
    @printf("Munition : %s\n", AMMO.name)
    @printf("Canon : L=%.2f m, Ø ext %.0f mm, alésage %.1f mm\n",
            BARREL.L, BARREL.D_out*1e3, AMMO.D_bore*1e3)
    @printf("Temps de sortie t_b = %.2f ms ; cible d'accord : θ̇(t_b) = 0 (bouche stationnaire)\n\n",
            exit_time()*1e3)

    # ==================================================================
    # (1) LE RÉGLAGE RÉEL : LA POSITION (masse fixe, on visse le tuner)
    # ==================================================================
    m = TUNER_MASS_FIX
    ds = 0.0:OVERHANG_STEP:OVERHANG_MAX
    θd = [theta_dot_tb(m; d_overhang = d)[1] for d in ds]
    sweet_d = sweet_spots(collect(ds), θd)

    println("(1) RÉGLAGE EN POSITION — tuner $(Int(m*1e3)) g (masse FIXE), vissé en porte-à-faux")
    println("    On accorde en le vissant devant la bouche — réglage continu, fin,")
    @printf("    repéré par tours de filetage (ici %.1f mm/tour).\n\n", THREAD_PITCH*1e3)
    @printf("    %-14s | %s\n", "porte-à-faux", "θ̇(t_b) MOA/ms")
    println("    " * "-"^36)
    for (d, θ) in zip(ds, θd)
        (round(d*1e3*2) % 10 == 0) && @printf("    %-14.0f | %+8.3f%s\n",
            d*1e3, θ, abs(θ) < 0.04 ? "   ← θ̇≈0 (sweet spot)" : "")
    end
    println("\n" * "-"^74)
    if isempty(sweet_d)
        println("Aucun sweet spot en position pour cette masse — en essayer une autre.")
    else
        println("SWEET SPOT(S) EN POSITION (θ̇ = 0, bouche stationnaire) :")
        for (k, d) in enumerate(sweet_d)
            @printf("  #%d : porte-à-faux ≈ %.1f mm\n", k, d*1e3)
        end
        # Largeur de l'optimum : plage de d où |θ̇| < 0,04 MOA/ms
        band = [d for (d, θ) in zip(ds, θd) if abs(θ) < 0.04]
        if !isempty(band)
            w = (maximum(band) - minimum(band)) * 1e3
            @printf("  Optimum LARGE : |θ̇|<0,04 sur ≈ %.0f mm (%.0f tours) — tune tolérant,\n",
                    w, w / (THREAD_PITCH*1e3))
            println("  la résolution fine du filetage sert à se caler au centre, sans")
            println("  qu'un écart de quelques tours ne gâche le groupement.")
        end
        println("  → UN sweet spot sur la course : on le trouve en scannant, exactement")
        println("    comme la procédure d'accord réelle (séries à N tours d'intervalle).")
    end

    # ==================================================================
    # (2) LA COURBE D'ACCORD : où tombe le sweet spot selon la masse fixée
    # ==================================================================
    println("\n(2) COURBE D'ACCORD — position du sweet spot selon la masse (le poids qu'on monte)")
    @printf("    %-12s | %s\n", "masse (g)", "porte-à-faux du sweet spot")
    println("    " * "-"^44)
    for mg in 40:10:100
        θ = [theta_dot_tb(mg*1e-3; d_overhang = d)[1] for d in ds]
        ss = sweet_spots(collect(ds), θ)
        @printf("    %-12d | %s\n", mg,
                isempty(ss) ? "— (aucun sur la course)" :
                join([@sprintf("%.0f mm", s*1e3) for s in ss], ", "))
    end
    println("    Une masse plus lourde rapproche le sweet spot de la bouche. Masse et")
    println("    position forment donc UN espace d'accord : le poids fixe la courbe,")
    println("    la position accorde dessus.")

    # (contexte) balayage en masse pur, pour mémoire
    masses = 0.0:TUNER_STEP:TUNER_MAX
    θm = [theta_dot_tb(mm)[1] for mm in masses]
    sweet_m = sweet_spots(collect(masses), θm)
    println()
    @printf("Pour mémoire, en jouant sur la MASSE (tuner ponctuel, d=0) : sweet spots à %s g.\n",
            isempty(sweet_m) ? "—" : join([@sprintf("%.0f", s*1e3) for s in sweet_m], ", "))
    println("La masse est un knob comparable, mais ce n'est pas ainsi qu'on règle : on")
    println("fixe le poids et on visse. La position offre la résolution fine (par tours).")
    println()
    println("Rappel : la POSITION des sweet spots ne dépend PAS de l'amplitude d'excitation")
    println("(linéarité) — on prédit le « où » sans la mesurer ; seule la PROFONDEUR du")
    println("gain en dépend.")
    println("="^74)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
