# =============================================================================
# LE TUNER À TUBE N'EST PAS UN CORPS RIGIDE DANS LA BANDE QUI COMPTE
#
# `simulation.jl` modélise le tuner en MASSE PONCTUELLE affectée d'une inertie,
# attachée à la bouche avec un porte-à-faux. Cette idéalisation suppose que
# l'ensemble tube reste rigide aux fréquences de travail. Elle était justifiée
# tant qu'on croyait l'accord gouverné par la fondamentale du canon (35 Hz) :
# un tube résonnant bien plus haut se comporte alors en bloc.
#
# La décomposition modale du 2026-07-19 a montré que l'accord se joue en réalité
# sur des modes à 2-3 kHz. Or l'ensemble tube + masse coulissante résonne, en
# porte-à-faux depuis la bouche :
#
#     50 mm / 100 g  →  5514 Hz   (au-dessus de la bande : rigide, hypothèse OK)
#    100 mm / 100 g  →  1838 Hz   (EN PLEIN DEDANS : hypothèse fausse)
#    150 mm / 200 g  →   728 Hz   (en dessous : le tube suit à peine le canon)
#
# Aux porte-à-faux que le modèle utilise (~100-120 mm), l'hypothèse de corps
# rigide est donc violée. Ce script remplace la masse ponctuelle par un TUBE
# ÉLASTIQUE maillé en poutre, avec la masse coulissante placée le long de ce
# tube, et mesure l'écart.
#
# Le tube acier Ø40 × 1,25 mm est PLUS RAIDE que le canon (EI = 5718 contre
# 3256 N·m²) et trois fois plus léger — ce n'est pas un détail : il ne se
# contente pas de suivre la bouche, il en modifie la dynamique.
#
# Usage :   julia flexible_tube.jl
# =============================================================================

include(joinpath(@__DIR__, "simulation.jl"))
using Printf, LinearAlgebra

# Géométrie du tube (mêmes cotes que la section TUBE_* de simulation.jl)
const TUBE_I    = π/4 * (TUBE_RO^4 - TUBE_RI^4)
const N_TUBE    = 8          # éléments de tube

# MATÉRIAU DU TUBE. simulation.jl suppose de l'ACIER, choix hérité du calage de
# k = 5 cm à 200 g. Ce choix a une conséquence qu'il n'avait pas vue : à 1,195
# kg/m, 200 g de tuner ne font que 167 mm de tube et il ne reste RIEN pour la
# masse coulissante. Le modèle décrivait donc un tube nu, tout en prétendant en
# régler la position. Les tubes Starik/Centra sont en ALUMINIUM (0,411 kg/m) :
# 200 mm de tube n'y pèsent que 82 g, laissant 118 g de curseur — un tuner
# réellement réglable.
const TUBE_MATERIALS = (acier = (ρ = ρ_steel, E = 200e9),
                        alu   = (ρ = 2700.0,  E = 70e9))
tube_props(mat) = (E = TUBE_MATERIALS[mat].E * TUBE_I,
                   ρA = TUBE_MATERIALS[mat].ρ * TUBE_A)
const TUBE_EI   = tube_props(:acier).E
const TUBE_RHOA = tube_props(:acier).ρA

# Assemblage canon + tube élastique. Le nœud de bouche (fin d'alésage) reste
# celui où θ est lu : c'est là que la balle sort, le tube ne la guide plus.
function build_flex(m_slider, L_tube, d_slider; mat = :acier)
    n_nodes = N_nodes + N_TUBE
    nd = 2 * n_nodes
    K = zeros(nd, nd); M = zeros(nd, nd)
    Kb, Mb = element_matrices(L_e, EI, ρA)
    for e in 1:N_elements
        idx = (2*e-1):(2*e+2)
        @views K[idx, idx] .+= Kb
        @views M[idx, idx] .+= Mb
    end
    Lt_e = L_tube / N_TUBE
    tp = tube_props(mat)
    Kt, Mt = element_matrices(Lt_e, tp.E, tp.ρA)
    for k in 1:N_TUBE
        n1 = N_nodes + k - 1
        idx = (2*n1-1):(2*n1+2)
        @views K[idx, idx] .+= Kt
        @views M[idx, idx] .+= Mt
    end
    # masse coulissante : au nœud de tube le plus proche de d_slider
    k_sl = clamp(round(Int, d_slider / Lt_e), 0, N_TUBE)
    n_sl = N_nodes + k_sl
    M[2*n_sl-1, 2*n_sl-1] += m_slider
    active = 3:nd
    return K[active, active], M[active, active], 2*N_nodes - 2   # ddl y de bouche
end

function shoot_flex(m_slider; L_tube = 0.12, d_slider = 0.10, h = H_OFFSET_EFF,
                    Δt = 2e-6, t_end = 6e-3, ζ1 = 0.005, ζ2 = 0.005, mat = :acier)
    Ka, Ma, dy = build_flex(m_slider, L_tube, d_slider; mat = mat)
    freqs, ωs, _ = modal_analysis(Ka, Ma; n_modes = 12)
    Ca, _, _ = rayleigh_damping(Ma, Ka, ωs[1], ωs[2], ζ1, ζ2)
    kin = projectile_kinematics(v_muzzle, L)
    nd = size(Ka, 1)
    F = t -> begin
        f = zeros(nd)
        f[2] += kin.p(t) * A_bore * h
        xp = kin.x(t)
        if 0 < xp < L
            e = clamp(Int(floor(xp / L_e)) + 1, 1, N_elements)
            ξ = (xp - (e-1)*L_e) / L_e
            Ns = (1-3ξ^2+2ξ^3, L_e*(ξ-2ξ^2+ξ^3), 3ξ^2-2ξ^3, L_e*(-ξ^2+ξ^3))
            for (k, ig) in enumerate((2*e-1, 2*e, 2*e+1, 2*e+2))
                ia = ig - 2
                1 <= ia <= nd && (f[ia] += -m_p * g_accel * Ns[k])
            end
        end
        f
    end
    ts, U, V, _ = newmark_solve(Ma, Ca, Ka, F, t_end, Δt)
    ib = argmin(abs.(ts .- kin.t_b))
    (freqs = freqs, θ = U[dy+1, ib], θdot = moa_per_ms(V[dy+1, ib]),
     pic = maximum(abs.(moa_per_ms.(V[dy+1, :]))))
end

function main()
    println("="^76)
    println(" Tuner à tube : masse ponctuelle rigide contre tube élastique")
    println("="^76)
    @printf("\nTube Ø%.0f × %.2f mm : EI = %.0f N·m² (canon : %.0f), %.3f kg/m\n\n",
            TUBE_OD*1e3, TUBE_WALL*1e3, TUBE_EI, EI, TUBE_RHOA)

    # COMPARAISON À MASSE ÉGALE. Le modèle rigide déduit son inertie d'un tube
    # de longueur ℓ = m/(ρ·A) : à 200 g cela décrit 167 mm de tube UNIFORME, sans
    # masse coulissante. Comparer à un tube de 120 mm PLUS 200 g de curseur
    # confronterait deux objets différents (343 g contre 200 g) et l'écart mesuré
    # ne voudrait rien dire. On impose donc le même budget de masse : le tube
    # consomme ρA·L_tube, le curseur reçoit le reste.
    println("  masse | L_tube | curseur | porte-à-faux |  RIGIDE θ̇  | ÉLASTIQUE θ̇ |  écart")
    println("  " * "-"^78)
    for m in (0.100, 0.200), d in (0.06, 0.10, 0.14)
        Lt = tuner_length(m)                    # longueur du tube équivalent
        m_tube = TUBE_RHOA * Lt
        m_sl = max(m - m_tube, 0.0)             # reste pour le curseur
        d > Lt && continue                      # porte-à-faux hors du tube
        rr = simulate_shot(m; d_overhang = d, h_offset = H_OFFSET_EFF, verbose = false)
        rf = shoot_flex(m_sl; L_tube = Lt, d_slider = d)
        @printf("  %4.0f g | %5.0f mm | %6.0f g |    %5.0f mm  | %+10.3f | %+11.3f | %+7.3f\n",
                m*1e3, Lt*1e3, m_sl*1e3, d*1e3, rr.θdot_MOAms, rf.θdot,
                rf.θdot - rr.θdot_MOAms)
    end

    # Tube ALUMINIUM : le seul qui laisse de la masse au curseur, donc le seul
    # où « régler la position » a un sens physique.
    println("\n  Tube ALUMINIUM, budget 200 g (tube 200 mm = 82 g, curseur 118 g) :")
    println("  position du curseur |  θ̇(t_b)  | θ(t_b) µrad")
    println("  " * "-"^48)
    for d in 0.02:0.03:0.20
        rf = shoot_flex(0.118; L_tube = 0.20, d_slider = d, mat = :alu)
        @printf("  %14.0f mm  | %+8.3f | %+11.1f\n", d*1e3, rf.θdot, rf.θ*1e6)
    end

    println("\nFréquences propres avec tube élastique (200 g à 100 mm) :")
    r = shoot_flex(0.0; L_tube = tuner_length(0.200), d_slider = 0.10)
    println("  " * join([@sprintf("%.0f", f) for f in r.freqs[1:8]], "  ") * " Hz")
    rr = simulate_shot(0.200; d_overhang = 0.10, h_offset = H_OFFSET_EFF, verbose = false)
    println("  rigide, pour comparaison :")
    println("  " * join([@sprintf("%.0f", f) for f in rr.freqs[1:5]], "  ") * " Hz")
    println("\n" * "-"^76)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
