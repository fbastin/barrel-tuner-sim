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
# CAVEAT « masse vs position ». Un vrai tuner airgun s'ajuste en POSITION (on
# dévisse le poids) : il change le bras de levier ET l'inertie, plus vite que
# la seule masse. Ce script balaie la MASSE — il montre le mécanisme et donne
# une borne BASSE de la sensibilité ; un balayage en position serait plus fin.
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

function build_system(m_tuner)
    Ke, Me = element_matrices()
    K = zeros(ndof, ndof); M = zeros(ndof, ndof)
    for e in 1:BARREL.N
        idx = (2*e - 1):(2*e + 2)
        @views K[idx, idx] .+= Ke
        @views M[idx, idx] .+= Me
    end
    active = 3:ndof                       # encastrement culasse : y₁ = θ₁ = 0
    Ka = K[active, active]; Ma = copy(M[active, active])
    Ma[end-1, end-1] += m_tuner           # masse ponctuelle à la bouche
    Ma[end,   end]   += m_tuner * 2.0e-3  # inertie du tuner (rayon giration ~4,5 cm)
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

# θ̇(L, t_b) en MOA/ms pour une masse de tuner donnée
function theta_dot_tb(m_tuner)
    Ka, Ma = build_system(m_tuner)
    ωs = modal(Ka, Ma; n = 4)
    Ca = rayleigh(Ma, Ka, ωs[1], ωs[2], ζ_STEEL)   # canon nu : modes de flexion
    ndof_a = size(Ka, 1)
    F_of_t = t -> (F = zeros(ndof_a); F[2] += hammer_moment(t); F)  # moment à la culasse
    t_b = exit_time()
    ts, U, V = newmark(Ma, Ca, Ka, F_of_t, 1.5*t_b)
    ib = argmin(abs.(ts .- t_b))
    return moa_per_ms(V[end, ib]), ωs[1]/(2π)      # θ̇(t_b), f₁
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

    masses = 0.0:TUNER_STEP:TUNER_MAX
    θdots  = Float64[]; f1s = Float64[]
    for m in masses
        θd, f1 = theta_dot_tb(m); push!(θdots, θd); push!(f1s, f1)
    end

    # Détection des passages par zéro (sweet spots) par interpolation linéaire
    sweet = Float64[]
    for i in 1:length(masses)-1
        if θdots[i] == 0 || sign(θdots[i]) != sign(θdots[i+1])
            m0 = masses[i] - θdots[i]*(masses[i+1]-masses[i])/(θdots[i+1]-θdots[i])
            push!(sweet, m0)
        end
    end

    @printf("f₁ : %.1f Hz (nu) → %.1f Hz (%.0f g)\n\n", f1s[1], f1s[end], TUNER_MAX*1e3)
    println("Balayage θ̇(t_b) vs masse de tuner (extrait tous les 10 g) :")
    @printf("  %-10s | %s\n", "tuner (g)", "θ̇(t_b) MOA/ms")
    println("  " * "-"^32)
    for (m, θd) in zip(masses, θdots)
        (round(m*1e3) % 10 == 0) && @printf("  %-10.0f | %+8.3f%s\n",
            m*1e3, θd, abs(θd) < 0.05 ? "   ← θ̇≈0 (sweet spot)" : "")
    end

    println("\n" * "-"^74)
    if isempty(sweet)
        println("Aucun sweet spot (θ̇ = 0) dans 0–$(Int(TUNER_MAX*1e3)) g :")
        println("la bouche ne repasse pas par un rebroussement sur cette plage.")
    else
        println("SWEET SPOTS (θ̇ = 0, bouche stationnaire) :")
        for (k, m) in enumerate(sweet)
            @printf("  #%d : tuner ≈ %5.1f g%s\n", k, m*1e3,
                    m <= 0.080 ? "   (dans la plage réaliste ≤ 80 g)" : "   (au-delà du réaliste)")
        end
        if length(sweet) >= 2
            gaps = diff(sweet) .* 1e3
            @printf("\n  Espacement des sweet spots : %s g  (moyenne %.0f g).\n",
                    join([@sprintf("%.0f", g) for g in gaps], ", "), sum(gaps)/length(gaps))
            println("  Très espacés : sur la course réaliste d'un tuner, on n'en croise")
            println("  qu'UN — on le trouve en scannant, comme en pratique. La récurrence")
            println("  est robuste (décalage de phase induit par le tuner), mais l'amplitude")
            println("  d'excitation, inconnue, n'en fixe pas la position — seulement la")
            println("  profondeur du gain.")
        end
    end
    println("="^74)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
