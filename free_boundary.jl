# =============================================================================
# L'ENCASTREMENT EST-IL LE MÉCANISME MANQUANT ? — RÉPONSE : NON
#
# [C ter] désigne depuis le début la rotation d'ensemble supprimée par
# l'encastrement de culasse comme le mécanisme absent, et Kolbe l'écrit lui-même :
# « the recoil forces […] imparting a moment on the back of the barrel AS THE
# RIFLE ROTATES ABOUT ITS CENTRE OF GRAVITY ». Ce script rend au modèle les deux
# degrés de liberté de corps rigide et mesure ce qu'on y gagne.
#
# RÉSULTAT : on n'y gagne rien, et on y perd. Avec un bras de levier PHYSIQUE
# (h_cg = 25,4 mm, une vraie cote) le modèle libre donne θ̇(t_b) ≈ −0,18 MOA/ms et
# un pic de 0,49, contre −9,4 mesurés par Kolbe — soit ~20 fois trop peu. Le
# modèle encastré, avec son h_offset gonflé à 62,5 mm, en était PLUS proche.
# Libérer la culasse ne débloque donc pas l'amplitude.
#
# CE QUE L'EXPÉRIENCE ÉTABLIT AU PASSAGE :
#   • le décalage CG/culasse compte (θ̇ passe de −0,106 à −0,18 quand le CG recule
#     de 0 à 250 mm) mais d'un facteur 2, pas 20 : c'est la translation imposée à
#     la culasse quand l'arme pivote, et elle ne suffit pas ;
#   • la fréquence fondamentale ne bouge quasiment pas (34 → 54 Hz), donc
#     t_b/T₁ ≈ 0,15 : la balle sort toujours au tout début du premier cycle ;
#   • la rotation d'ensemble seule vaut ~0,2 MOA/ms, ce que kolbe_amplitude.jl
#     avait déjà obtenu analytiquement par une voie indépendante.
#
# OÙ CHERCHER ENSUITE. kolbe_amplitude.jl, qui modélise le fusil ENTIER sur deux
# sacs (crosse comprise, contact unilatéral), atteint un pic de 4,8 MOA/ms — dix
# fois ce script. La différence n'est donc pas la libération de la culasse mais
# la présence de la CROSSE et des APPUIS. C'est là qu'il faut regarder, non dans
# la condition aux limites du canon seul.
#
# Usage :   julia free_boundary.jl
# =============================================================================

include(joinpath(@__DIR__, "simulation.jl"))
using Printf, LinearAlgebra

# Fusil LIBRE : on rend les 2 ddl de corps rigide que l'encastrement supprime.
# Le corps de l'arme (masse et inertie) est attache au noeud de culasse, et
# l'excitation devient un MOMENT p*A_bore*h_cg autour de ce point -- h_cg est
# alors une COTE REELLE (25,4 mm), non un facteur de gain.
const M_BARREL = ρA * L
const M_BODY   = m_rifle_ref - M_BARREL        # arme moins canon
const K_GYR_RIFLE = 0.30                       # rayon de giration de l'arme (m)
const I_BODY   = m_rifle_ref * K_GYR_RIFLE^2
D_CG = -0.15                                   # CG en arriere de la culasse (m)

function build_free(m_tuner, J_tuner; d_overhang=0.0)
    Ke, Me = element_matrices(L_e, EI, ρA)
    K = zeros(ndof, ndof); M = zeros(ndof, ndof)
    for e in 1:N_elements
        idx = (2*e-1):(2*e+2)
        @views K[idx,idx] .+= Ke
        @views M[idx,idx] .+= Me
    end
    # Corps rigide DEPORTE : son CG est a D_CG de la culasse (negatif = en
    # arriere). Meme couplage masse/rotation que pour le tuner deporte. C'est ce
    # decalage qui fait TRANSLATER la culasse quand l'arme pivote, et donc qui
    # secoue le canon -- sans lui, la rotation d'ensemble ne fait que pivoter
    # l'encastrement et n'excite presque rien.
    M[1,1] += M_BODY
    M[1,2] += M_BODY*D_CG
    M[2,1] += M_BODY*D_CG
    M[2,2] += M_BODY*D_CG^2 + I_BODY
    d = d_overhang
    M[end-1,end-1] += m_tuner
    M[end-1,end  ] += m_tuner*d
    M[end,  end-1] += m_tuner*d
    M[end,  end  ] += m_tuner*d^2 + J_tuner
    return K, M
end

function shoot_free(m_tuner; d_overhang=0.0, h=h_cg_ref, Δt=5e-6, t_end=30e-3,
                    ζ1=0.005, ζ2=0.01, J_tuner=tuner_inertia(m_tuner))
    K, M = build_free(m_tuner, J_tuner; d_overhang)
    freqs, ωs, Φ = modal_analysis(K, M; n_modes=5)
    Ca, _, _ = rayleigh_damping(M, K, ωs[1], ωs[2], ζ1, ζ2)
    kin = projectile_kinematics(v_muzzle, L)
    nd = size(K,1)
    F = t -> begin
        f = zeros(nd)
        f[2] += kin.p(t) * A_bore * h        # moment de recul sur le CORPS
        xp = kin.x(t)
        if 0 < xp < L
            e = clamp(Int(floor(xp/L_e))+1, 1, N_elements)
            ξ = (xp - (e-1)*L_e)/L_e
            Ns = (1-3ξ^2+2ξ^3, L_e*(ξ-2ξ^2+ξ^3), 3ξ^2-2ξ^3, L_e*(-ξ^2+ξ^3))
            for (k,ig) in enumerate((2*e-1, 2*e, 2*e+1, 2*e+2))
                f[ig] += -m_p*g_accel*Ns[k]
            end
        end
        f
    end
    ts, U, V, _ = newmark_solve(M, Ca, K, F, t_end, Δt)
    ib = argmin(abs.(ts .- kin.t_b))
    # angle de bouche RELATIF au corps : c'est ce que voit la balle
    pic = maximum(abs.(moa_per_ms.(V[end,:])))
    (freqs=freqs, pic=pic, θ=U[end,ib], θdot=moa_per_ms(V[end,ib]),
     θ_body=U[2,ib], θdot_body=moa_per_ms(V[2,ib]), t_b=kin.t_b)
end

println("Effet du decalage CG/culasse (tuner 200 g, h = h_cg reel) :")
println("  D_CG (mm) | f1 (Hz) | theta bouche | theta_dot bouche")
for dc in (0.0, -0.05, -0.10, -0.15, -0.25, -0.40)
    global D_CG = dc
    q = shoot_free(0.200)
    @printf("  %8.0f  | %6.1f  | %+9.1f    | %+8.3f\n", dc*1e3, q.freqs[1], q.θ*1e6, q.θdot)
end
println()
global D_CG = -0.15
r = shoot_free(0.200)
@printf("Frequences libres (Hz) : %s\n", join([@sprintf("%.1f", f) for f in r.freqs[1:5]], "  "))
@printf("(encastre, pour memoire : 34,2  217,7  602,4 ...)\n\n")
@printf("h = h_cg REEL = %.1f mm (et non un facteur de gain)\n", h_cg_ref*1e3)
@printf("  bouche : θ = %+.1f urad, θdot = %+.3f MOA/ms\n", r.θ*1e6, r.θdot)
@printf("  corps  : θ = %+.1f urad, θdot = %+.3f MOA/ms\n", r.θ_body*1e6, r.θdot_body)
println()
println("Balayage en position, tuner 100 g :")
for d in 0.0:0.02:0.16
    q = shoot_free(0.100; d_overhang=d)
    @printf("  %5.0f mm : θ = %+8.1f urad, θdot = %+7.3f MOA/ms\n", d*1e3, q.θ*1e6, q.θdot)
end
q0 = shoot_free(0.0)
@printf("\ncanon nu : θdot = %+.3f MOA/ms   (Kolbe mesure -9,4)\n", q0.θdot)

println()
println("PIC de |theta_dot| contre valeur A LA SORTIE (canon nu, h reel) :")
z = shoot_free(0.0)
@printf("  pic          = %6.2f MOA/ms   (Kolbe mesure -9,4 A LA SORTIE)\n", z.pic)
@printf("  a t_b        = %+6.3f MOA/ms\n", z.θdot)
@printf("  rapport      = %.0f\n", z.pic/abs(z.θdot))
@printf("  periode mode 1 = %.2f ms, t_b = %.2f ms -> t_b/T1 = %.2f\n",
        1000/z.freqs[1], z.t_b*1e3, z.t_b*z.freqs[1])
