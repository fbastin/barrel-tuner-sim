# =============================================================================
# Étape 2/3 du chantier « modèle de Harral » : LE FUSIL ENTIER SUR SES DEUX SACS
#
# But : tester si la ROTATION DE CORPS RIGIDE de l'arme restaure la
# compensation positive — c.-à-d. si la dispersion CROÎT avec la masse de
# bouche, comme chez Harral (0.091 → 0.183 in), au lieu de rester plate à la
# valeur non compensée (≈ 0.173) comme dans `harral_a22lr.jl` (culasse
# encastrée).
#
# ⚠️ CE N'EST PAS UNE REPRODUCTION. Harral ne publie PAS les données requises :
# sur l'arme entière il ne donne qu'un chiffre, 10.5 lb. Ni crosse (matériau,
# raideur, masse), ni position des sacs, ni hauteur d'âme, ni amortissement.
# Caler ces inconnues pour retomber sur ses chiffres serait du curve-fitting.
# On les BALAIE donc, et on ne regarde qu'une chose : LA STRUCTURE (spread
# croissante avec la masse) tient-elle sur toute la plage ?
# Le succès NE se mesure PAS à l'écart aux chiffres de Harral.
#
# Ancrage réel (gratuit) : le budget de masse. Canon = 1.711 kg ; fusil total
# = 10.5 lb = 4.763 kg ⇒ boîte + crosse + avant-bras = 3.052 kg (64 % de l'arme).
#
# -----------------------------------------------------------------------------
# HISTORIQUE — une version antérieure de ce script était FAUSSE. Elle reportait
# le sac avant EN ARRIÈRE de la culasse, sur la crosse, faute de modéliser
# l'avant-bras. Le sac avant se retrouvait alors DERRIÈRE le centre de gravité
# (CG ≈ 16.8-19.5 in du talon selon la culasse ; sacs testés : 11.7-20.2 in) :
# l'arme BASCULAIT vers l'avant et n'était retenue que par un sac arrière
# TIRANT VERS LE BAS. Avec des appuis bilatéraux, c'est silencieux mais
# physiquement impossible ; le contact unilatéral l'a exposé en faisant
# décoller l'appui, laissant l'arme sur un seul point ⇒ raideur singulière ⇒
# divergence (spreads à 1e8 in). Le « 74 % de configs conformes à Harral » issu
# de cette version est CONTAMINÉ et ne doit pas être cité.
#
# CORRECTIF (version actuelle) : l'avant-bras est modélisé pour ce qu'il est —
# un membre PARALLÈLE au canon (canon flottant), partant de la boîte vers
# l'avant, avec le sac avant sous son extrémité. Structure RAMIFIÉE : deux
# chaînes d'éléments (canon et avant-bras) partagent le nœud de culasse. Le CG
# se retrouve ainsi ENTRE les deux sacs, comme sur une vraie arme.
# -----------------------------------------------------------------------------
#
# AUTRES CHOIX DE MODÉLISATION (et leurs limites) :
#  * Poutre plane. Harral maille en briques solides 3D à 8 nœuds.
#  * Le canon reste en PORTE-À-FAUX au-delà de l'appui avant : c'est cela qui
#    produit l'affaissement de bouche (les ≈ −1.32 in de Harral).
#  * Appuis = ressorts verticaux (sacs), en COMPRESSION SEULE par défaut
#    (`unilateral=true`) : le talon peut décoller. L'option `unilateral=false`
#    permet de mesurer l'effet de cette hypothèse — c'est elle qui masquait le
#    défaut ci-dessus.
#  * Le moment de recul est ANCRÉ PHYSIQUEMENT : la force axiale p(t)·A_âme
#    s'applique à la culasse, à la hauteur h_bore au-dessus de la ligne des
#    sacs ⇒ moment de tangage p·A·h_bore. Avec la culasse encastrée, ce bras de
#    levier n'avait aucune référence et flottait librement.
#  * Crosse + boîte + avant-bras : masse et EI UNIFORMES (l'acier de la boîte
#    est en réalité bien plus dense que le bois). Simplification assumée.
#
# Usage :  julia harral_rifle_sweep.jl
# =============================================================================

include(joinpath(@__DIR__, "harral_a22lr.jl"))   # géométrie canon, pression,
                                                 # Newmark, POI, exit_time…
using Printf
using LinearAlgebra

# -----------------------------------------------------------------------------
# Budget de masse (ancrage sur le seul chiffre publié)
# -----------------------------------------------------------------------------
const M_RIFLE_TOTAL = 10.5 * 0.45359237          # kg (Harral : « 10.5# »)
const M_BARREL      = sum(begin
                              A, _ = element_section(e)
                              ρ_steel * A * L_e
                          end for e in 1:N_elements)
const M_STOCK       = M_RIFLE_TOTAL - M_BARREL   # boîte + crosse + avant-bras

# -----------------------------------------------------------------------------
# EXCITATION PHYSIQUE : le moment de recul est calé sur l'impulsion RÉELLE de
# recul, J = (1+β_gaz)·m_p·v_muzzle, et non sur le p_max = 200 MPa de
# `chamber_pressure` (un simple gabarit de forme, ~2,6× le recul réel). On
# garde la FORME temporelle de chamber_pressure et on la rééchelonne pour que
# ∫ A_bore·p dt = J. Ainsi l'amplitude devient une PRÉDICTION (via h_bore,
# hauteur d'âme), plus un paramètre libre.
# -----------------------------------------------------------------------------
const β_GAZ = 0.3        # part des gaz dans l'impulsion de recul (~+30 %, rimfire)
const PRESSURE_SHAPE_INT = let s = 0.0, dt = 1e-7, t = 0.0
    while t < 6 * t_peak; s += chamber_pressure(t) * dt; t += dt; end
    s                                    # ∫ chamber_pressure dt (à p_max), Pa·s
end

const N_STOCK = 16      # talon → culasse
const N_BAR   = N_elements   # culasse → bouche (48)
const N_FORE  = 8       # culasse → pointe d'avant-bras

# -----------------------------------------------------------------------------
# TOPOLOGIE RAMIFIÉE
#   chaîne crosse  : nœuds 1 … N_STOCK+1        (x = 0 … x_breech)
#   chaîne canon   : repart du nœud de culasse  (x = x_breech … x_breech+L)
#   chaîne av.-bras: repart du nœud de culasse  (x = x_breech … x_breech+L_fore)
# Les deux dernières sont PARALLÈLES : elles ne se touchent qu'à la culasse.
# -----------------------------------------------------------------------------
# BRANCHE TUBE (ajoutée le 2026-07-19). Un tuner à tube n'est pas une masse
# ponctuelle : c'est une poutre chaînée en avant de la bouche, dans laquelle une
# masse coulisse. L'idéaliser en masse + inertie gonfle son effet d'un ordre de
# grandeur (cf. flexible_tube.jl), parce qu'une masse rigide déportée transmet à
# la bouche un moment `m·d·ÿ + (m d² + J)·θ̈` qu'un tube fléchissant ne transmet
# pas — et dans la bande qui gouverne l'accord (2-3 kHz) il fléchit.
# `L_tube = 0` restitue exactement le comportement antérieur.
const N_TUBE_EL = 6

function rifle_topology(x_breech, L_fore; L_tube = 0.0)
    n_tube = L_tube > 0 ? N_TUBE_EL : 0
    nn = (N_STOCK + 1) + N_BAR + N_FORE + n_tube
    xs = zeros(nn)
    for i in 1:(N_STOCK+1)
        xs[i] = (i-1) * x_breech / N_STOCK
    end
    i_breech = N_STOCK + 1
    for k in 1:N_BAR
        xs[i_breech + k] = x_breech + k * L / N_BAR
    end
    off = i_breech + N_BAR
    for k in 1:N_FORE
        xs[off + k] = x_breech + k * L_fore / N_FORE
    end

    els = Tuple{Int,Int,Float64,Symbol}[]
    for e in 1:N_STOCK
        push!(els, (e, e+1, x_breech/N_STOCK, :stock))
    end
    prev = i_breech
    for k in 1:N_BAR
        push!(els, (prev, i_breech+k, L/N_BAR, :barrel)); prev = i_breech+k
    end
    prev = i_breech
    for k in 1:N_FORE
        push!(els, (prev, off+k, L_fore/N_FORE, :fore)); prev = off+k
    end
    i_muzzle = i_breech + N_BAR
    i_tubetip = i_muzzle
    if n_tube > 0
        off_t = off + N_FORE
        for k in 1:n_tube
            xs[off_t + k] = x_breech + L + k * L_tube / n_tube
        end
        prev = i_muzzle
        for k in 1:n_tube
            push!(els, (prev, off_t+k, L_tube/n_tube, :tube)); prev = off_t+k
        end
        i_tubetip = off_t + n_tube
    end
    return (xs=xs, els=els, nn=nn, i_breech=i_breech,
            i_muzzle=i_muzzle, i_foretip=off+N_FORE, i_tubetip=i_tubetip,
            n_tube=n_tube)
end

dofs_of(n1, n2) = (2*n1-1, 2*n1, 2*n2-1, 2*n2)

# `d_overhang` : porte-à-faux du centre de masse du tuner DEVANT la bouche.
# AJOUTÉ le 2026-07-19. Jusqu'ici le tuner était une simple masse ponctuelle à la
# bouche : le modèle d'arme complète ne disposait donc que du levier FAIBLE — la
# masse — et ignorait le réglage de terrain, qui est la POSITION. C'est le
# couplage masse/rotation d'une masse déportée, identique à celui de
# simulation.jl, où le balayage en position s'était révélé bien plus autoritaire
# que le balayage en masse.
# Tube du tuner : aluminium par défaut, comme les Starik/Centra. L'acier, hérité
# du calage de k = 5 cm dans simulation.jl, est intenable : à 1,195 kg/m un tuner
# de 200 g ne fait que 167 mm de tube et ne laisse RIEN au curseur — on décrivait
# un tube nu tout en prétendant en régler la position.
const TUBE_OD_R, TUBE_WALL_R = 0.040, 0.00125
const TUBE_RO_R = TUBE_OD_R/2
const TUBE_RI_R = TUBE_RO_R - TUBE_WALL_R
const TUBE_A_R  = π * (TUBE_RO_R^2 - TUBE_RI_R^2)
const TUBE_I_R  = π/4 * (TUBE_RO_R^4 - TUBE_RI_R^4)
const TUBE_EI_R   = 70e9   * TUBE_I_R      # aluminium
const TUBE_RHOA_R = 2700.0 * TUBE_A_R

function build_rifle(; m_tuner, x_breech, L_fore, x_rear, EI_stock, k_rest,
                     d_overhang = 0.0, L_tube = 0.0, m_slider = 0.0, d_slider = 0.0)
    T  = rifle_topology(x_breech, L_fore; L_tube = L_tube)
    nd = 2 * T.nn
    K  = zeros(nd, nd); M = zeros(nd, nd); Fg = zeros(nd)
    ρA_stock = M_STOCK / (x_breech + L_fore)      # masse répartie sur crosse+av.-bras

    for (n1, n2, le, kind) in T.els
        if kind === :barrel
            x_mid = 0.5*(T.xs[n1] + T.xs[n2]) - x_breech   # abscisse LOCALE canon
            D = D_out_of_x(x_mid)
            A = π/4  * (D^2 - D_bore^2)
            I = π/64 * (D^4 - D_bore^4)
            EI, ρA = E * I, ρ_steel * A
        elseif kind === :tube
            EI, ρA = TUBE_EI_R, TUBE_RHOA_R
        else
            EI, ρA = EI_stock, ρA_stock
        end
        Ke, Me = element_matrices(le, EI, ρA)
        idx = collect(dofs_of(n1, n2))
        @views K[idx, idx] .+= Ke
        @views M[idx, idx] .+= Me
        q  = -ρA * g_ms
        @views Fg[idx] .+= q * le * [0.5, le/12, 0.5, -le/12]
    end

    # Masse coulissante DANS le tube (modèle élastique), au nœud le plus proche.
    if T.n_tube > 0 && m_slider > 0
        k_sl = clamp(round(Int, d_slider / (L_tube / T.n_tube)), 0, T.n_tube)
        n_sl = k_sl == 0 ? T.i_muzzle : (T.i_tubetip - T.n_tube + k_sl)
        M[2*n_sl-1, 2*n_sl-1] += m_slider
        Fg[2*n_sl-1]          += -m_slider * g_ms
    end

    dm = 2*T.i_muzzle - 1
    d  = d_overhang
    M[dm,   dm  ] += m_tuner
    M[dm,   dm+1] += m_tuner * d
    M[dm+1, dm  ] += m_tuner * d
    M[dm+1, dm+1] += m_tuner * d^2 + tuner_inertia(m_tuner)
    Fg[dm]        += -m_tuner * g_ms
    Fg[dm+1]      += -m_tuner * g_ms * d

    i_rear = argmin(abs.(T.xs[1:(N_STOCK+1)] .- x_rear))   # sac arrière : sur la crosse
    spring_dofs = [2*i_rear - 1, 2*T.i_foretip - 1]        # sac avant : pointe d'avant-bras

    return (K=K, M=M, Fg=Fg, T=T, nd=nd, spring_dofs=spring_dofs)
end

# Centre de gravité (contrôle de stabilité : doit tomber ENTRE les deux sacs)
function rifle_cg(x_breech, L_fore)
    T = rifle_topology(x_breech, L_fore)
    ρA_stock = M_STOCK / (x_breech + L_fore)
    num = 0.0; den = 0.0
    for (n1, n2, le, kind) in T.els
        m = (kind === :barrel) ? begin
                x_mid = 0.5*(T.xs[n1]+T.xs[n2]) - x_breech
                D = D_out_of_x(x_mid); ρ_steel * π/4 * (D^2 - D_bore^2) * le
            end : ρA_stock * le
        num += m * 0.5*(T.xs[n1]+T.xs[n2]); den += m
    end
    return num/den
end

# -----------------------------------------------------------------------------
# CONTACT UNILATÉRAL : la surface du sac est fixe en y = 0 ; l'arme repose
# DESSUS. Un ressort ne travaille qu'en compression (nœud pénétrant, y < 0).
# Si y > 0, l'appui a décollé et ne rend rien. 2 appuis ⇒ 4 états (masque).
# -----------------------------------------------------------------------------
function contact_mask(u, spring_dofs; unilateral = true)
    unilateral || return (1 << length(spring_dofs)) - 1
    m = 0
    for (j, d) in enumerate(spring_dofs)
        u[d] < 0 && (m |= (1 << (j-1)))
    end
    return m
end

function K_with_mask(K0, spring_dofs, k_rest, mask)
    Kk = copy(K0)
    for (j, d) in enumerate(spring_dofs)
        ((mask >> (j-1)) & 1) == 1 && (Kk[d, d] += k_rest)
    end
    return Kk
end

# Équilibre statique sous gravité. Renvoie aussi `stable` : si un appui décolle
# sous la SEULE gravité, l'arme bascule — la raideur devient singulière (mode de
# rotation libre) et toute valeur calculée serait du bruit. Config à REJETER.
function static_contact(K0, spring_dofs, k_rest, Fg; unilateral = true)
    full = (1 << length(spring_dofs)) - 1
    U = K_with_mask(K0, spring_dofs, k_rest, full) \ Fg
    stable = contact_mask(U, spring_dofs; unilateral) == full
    return U, full, stable
end

function newmark_contact(M, C, K0, spring_dofs, k_rest, F_of_t, t_end, Δt;
                         U0, unilateral = true)
    γ, β = 0.5, 0.25
    n  = size(M, 1)
    ts = collect(0:Δt:t_end); Nt = length(ts)
    U  = zeros(n, Nt); V = zeros(n, Nt); Ac = zeros(n, Nt)
    U[:, 1] = U0

    nmask = 1 << length(spring_dofs)
    Kmat = Dict(m => K_with_mask(K0, spring_dofs, k_rest, m) for m in 0:(nmask-1))
    # (M + γΔtC + βΔt²K) reste défini positif même sans appui : M régularise.
    Kfac = Dict(m => factorize(M + γ*Δt*C + β*Δt^2*Kmat[m]) for m in 0:(nmask-1))

    m0 = contact_mask(U0, spring_dofs; unilateral)
    Ac[:, 1] = M \ (F_of_t(ts[1]) - C*V[:,1] - Kmat[m0]*U[:,1])

    n_lift = 0
    for i in 1:Nt-1
        up = U[:,i] + Δt*V[:,i] + Δt^2*(0.5-β)*Ac[:,i]
        vp = V[:,i] + Δt*(1-γ)*Ac[:,i]
        Fi = F_of_t(ts[i+1])
        mask = contact_mask(U[:,i], spring_dofs; unilateral)
        a = zeros(n); u = zeros(n)
        for _ in 1:6                        # itération d'ensemble actif
            a = Kfac[mask] \ (Fi - C*vp - Kmat[mask]*up)
            u = up + β*Δt^2*a
            nm = contact_mask(u, spring_dofs; unilateral)
            nm == mask && break
            mask = nm
        end
        mask != (nmask-1) && (n_lift += 1)
        Ac[:, i+1] = a; U[:, i+1] = u; V[:, i+1] = vp + γ*Δt*a
    end
    return ts, U, V, n_lift / (Nt-1)
end

# -----------------------------------------------------------------------------
# Un tir sur le fusil entier
# -----------------------------------------------------------------------------
function point_load_barrel(Fval, x_global, T, x_breech, nd)
    le = L / N_BAR
    s  = x_global - x_breech
    (s <= 0 || s >= L) && return zeros(nd)
    k  = clamp(ceil(Int, s / le), 1, N_BAR)
    n1 = (k == 1) ? T.i_breech : T.i_breech + k - 1
    n2 = T.i_breech + k
    ξ  = (s - (k-1)*le) / le
    Ns = (1 - 3ξ^2 + 2ξ^3, le*(ξ - 2ξ^2 + ξ^3), 3ξ^2 - 2ξ^3, le*(-ξ^2 + ξ^3))
    F = zeros(nd)
    for (j, d) in enumerate(dofs_of(n1, n2))
        F[d] += Fval * Ns[j]
    end
    return F
end

# `p_of_t`, `x_of_t` et `t_b_override` permettent d'INJECTER une balistique
# intérieure externe (typiquement la version couplée de simulation.jl, où la
# cinématique est intégrée depuis la pression). Par défaut, comportement
# inchangé : gabarit `chamber_pressure` rééchelonné sur l'impulsion de recul et
# cinématique burnout locale. Ce point d'entrée existe pour confronter la
# STRUCTURE (crosse + appuis) à la meilleure EXCITATION disponible, les deux
# n'ayant jamais été combinées.
function shoot_rifle(v_muzzle; m_tuner, x_breech, L_fore, x_rear, EI_stock,
                     k_rest, h_bore, unilateral = true, ζ1 = 0.005, ζ2 = 0.005,
                     Δt = 2e-6, t_end = 3.0e-3,
                     p_of_t = nothing, x_of_t = nothing, t_b_override = nothing,
                     d_overhang = 0.0, return_state = false,
                     L_tube = 0.0, m_slider = 0.0, d_slider = 0.0)
    S = build_rifle(; m_tuner, x_breech, L_fore, x_rear, EI_stock, k_rest, d_overhang,
                    L_tube, m_slider, d_slider)
    xp = x_of_t === nothing ? projectile_pos(v_muzzle) : x_of_t

    U_stat, mask0, stable = static_contact(S.K, S.spring_dofs, k_rest, S.Fg; unilateral)
    stable || return nothing                 # l'arme bascule ⇒ config rejetée

    K_lin = K_with_mask(S.K, S.spring_dofs, k_rest, mask0)
    # ⚠️ AMORTISSEMENT : caler Rayleigh sur les modes de WHIP du canon
    # (flexion, > ~120 Hz), PAS sur les modes de corps rigide sur les sacs
    # (~20-80 Hz). Caler sur ces derniers imposait ζ ≈ 0,8 au whip (267 Hz) et
    # l'étouffait d'un facteur ~10 (bug corrigé le 2026-07-17). Les modes de
    # sacs, lents (< 0,1 cycle sur la fenêtre t_b), sont négligeables ici.
    fs, ωs = modal_analysis(K_lin, S.M; n_modes = 8)
    whip = findall(f -> f > 120.0, fs)
    ωa, ωb = length(whip) >= 2 ? (ωs[whip[1]], ωs[whip[2]]) : (ωs[end-1], ωs[end])
    f_whip = length(whip) >= 1 ? fs[whip[1]] : fs[end]
    C = rayleigh_damping(S.M, K_lin, ωa, ωb, ζ1, ζ2)

    d_breech_θ = 2 * S.T.i_breech
    # Rééchelonnement physique : ∫ A_bore·(p·p_scale) dt = (1+β)·m_p·v_muzzle.
    p_scale = (1 + β_GAZ) * m_p * v_muzzle / (A_bore * PRESSURE_SHAPE_INT)
    function F_of_t(t)
        F = copy(S.Fg)
        # Moment de recul physique. Avec une pression injectée, elle est déjà
        # cohérente en impulsion par construction : pas de rééchelonnement.
        pt = p_of_t === nothing ? chamber_pressure(t) * p_scale : p_of_t(t)
        F[d_breech_θ] += pt * A_bore * h_bore
        F .+= point_load_barrel(-m_p * g_ms, x_breech + xp(t), S.T, x_breech, S.nd)
        return F
    end

    ts, U, V, f_lift = newmark_contact(S.M, C, S.K, S.spring_dofs, k_rest,
                                       F_of_t, t_end, Δt; U0 = U_stat, unilateral)

    dm  = 2 * S.T.i_muzzle - 1
    t_b = t_b_override === nothing ? exit_time(v_muzzle) : t_b_override
    ib  = argmin(abs.(ts .- t_b))
    const_MOAms = (180*60/π) * 1e-3        # rad/s → MOA/ms
    if return_state
        # Pour la décomposition modale : on rend l'état complet, les modes et la
        # matrice de masse. Sert à savoir QUEL mode porte θ̇ à l'instant de sortie.
        fs_all, ωs_all, Φ_all = modal_analysis(K_lin, S.M; n_modes = 40, want_modes = true)
        return (ts = ts, U = U, V = V, ib = ib, dm = dm, M = S.M,
                fs = fs_all, Φ = Φ_all, t_b = t_b, f1 = f_whip)
    end
    return (θ_tb = U[dm+1, ib], ẏ_tb = V[dm, ib],
            θdot_tb   = V[dm+1, ib] * const_MOAms,          # taux angulaire à t_b, MOA/ms
            θdot_peak = maximum(abs.(V[dm+1, :])) * const_MOAms,
            θ_stat = U_stat[dm+1], f1 = f_whip, f_lift = f_lift)
end

function spread_for(m_tuner; kw...)
    r_lo = shoot_rifle(V_LO; m_tuner = m_tuner, kw...)
    r_hi = shoot_rifle(V_HI; m_tuner = m_tuner, kw...)
    (r_lo === nothing || r_hi === nothing) && return nothing, nothing
    poi(r, v) = r.θ_tb * D_TARGET_IN + (r.ẏ_tb * M2IN) * TOF(DROP[v]) - DROP[v]
    return abs(poi(r_lo, 1035) - poi(r_hi, 1075)), r_lo
end

# =============================================================================
# BALAYAGE
# =============================================================================
const TUNERS = [0.0, 4.9*OZ, 8.6*OZ, 16.0*OZ]
const HARRAL = [0.0910, 0.1218, 0.1554, 0.1832]

# Plages plausibles (pratique benchrest) — AUCUNE n'est publiée par Harral.
const X_BREECH = [18.0, 20.0, 22.0] .* IN   # culasse depuis le talon
const L_FORE   = [6.0, 10.0, 14.0] .* IN    # avant-bras EN AVANT de la culasse
const EI_STOCK = [3.0e3, 1.0e4, 4.0e4]      # raideur crosse/boîte/avant-bras (N·m²)
const H_BORE   = [1.0, 2.0, 3.0] .* IN      # âme au-dessus de la ligne des sacs
const K_REST   = [1.0e4, 1.0e5, 1.0e6]      # raideur d'un sac (N/m)
const X_REAR   = 1.5 * IN                   # sac arrière, près du talon

function run_sweep(; unilateral)
    rows = NamedTuple[]; n_reject = 0
    for xb in X_BREECH, lf in L_FORE, eis in EI_STOCK, hb in H_BORE, kr in K_REST
        kw = (x_breech=xb, L_fore=lf, x_rear=X_REAR, EI_stock=eis,
              k_rest=kr, h_bore=hb, unilateral=unilateral)
        sp = Float64[]; local r0; bad = false
        for m in TUNERS
            s, r = spread_for(m; kw...)
            if s === nothing || !isfinite(s); bad = true; break; end
            push!(sp, s); m == 0.0 && (r0 = r)
        end
        if bad; n_reject += 1; continue; end
        push!(rows, (xb=xb, lf=lf, eis=eis, hb=hb, kr=kr, sp=copy(sp),
                     Δ=sp[end]-sp[1], mono=all(diff(sp) .> 0),
                     proj=r0.θ_stat*D_TARGET_IN, f1=r0.f1, f_lift=r0.f_lift))
    end
    return rows, n_reject
end

function report(rows, n_reject, title)
    n = length(rows)
    if n == 0
        println("-"^78); println(title)
        @printf("  AUCUNE configuration valide (%d rejetées).\n", n_reject)
        return nothing
    end
    println("-"^78); println(title)
    @printf("Configurations valides : %d   (rejetées — l'arme bascule : %d)\n", n, n_reject)
    @printf("  tendance CROISSANTE (spread(16oz) > spread(nu)) : %d / %d  (%.0f %%)\n",
            count(r->r.Δ>0, rows), n, 100*count(r->r.Δ>0, rows)/n)
    @printf("  croissance STRICTEMENT MONOTONE sur les 4 masses : %d / %d  (%.0f %%)\n",
            count(r->r.mono, rows), n, 100*count(r->r.mono, rows)/n)
    @printf("  configs où un appui DÉCOLLE (>1 %% du temps) : %d / %d  (%.0f %%)\n",
            count(r->r.f_lift>0.01, rows), n, 100*count(r->r.f_lift>0.01, rows)/n)
    @printf("  flèche statique de bouche : %.2f … %.2f in   (Harral : −1.32 in)\n",
            minimum(r->r.proj, rows), maximum(r->r.proj, rows))
    println("-"^78)
    return rows[argmax([r.Δ for r in rows])]
end

function main_sweep()
    println("="^78)
    println(" Étape 2/3 — fusil entier sur deux sacs, avant-bras modélisé")
    println("="^78)
    @printf("Budget de masse : canon %.3f kg + (boîte+crosse+av.-bras) %.3f kg = %.3f kg (%.1f lb)\n",
            M_BARREL, M_STOCK, M_RIFLE_TOTAL, M_RIFLE_TOTAL/0.45359237)
    @printf("Harral (référence) : spread %.3f → %.3f in (CROISSANT)\n", HARRAL[1], HARRAL[end])
    @printf("Culasse encastrée (étape 1) : ≈ 0.171 → 0.164 in (plat, non compensé)\n\n")

    println("Contrôle de stabilité — le CG doit tomber ENTRE les deux sacs :")
    for xb in X_BREECH, lf in L_FORE
        cg = rifle_cg(xb, lf); xf = xb + lf
        @printf("  culasse %.0f in, avant-bras %.0f in → sac avant à %.1f in, CG à %.1f in : %s\n",
                xb*M2IN, lf*M2IN, xf*M2IN, cg*M2IN, cg < xf ? "OK" : "⚠ BASCULE")
    end
    println()

    rows_bi, rej_bi = run_sweep(unilateral = false)
    report(rows_bi, rej_bi, "APPUIS BILATÉRAUX (le sac peut TIRER — hypothèse à tester)")
    rows, rej = run_sweep(unilateral = true)
    best = report(rows, rej, "CONTACT UNILATÉRAL (compression seule — le talon décolle)")

    if best !== nothing && !isempty(rows_bi)
        println("\nEffet du contact unilatéral, par raideur de sac (fraction de Δ>0) :")
        for kr in K_REST
            b = filter(r->r.kr==kr, rows_bi); u = filter(r->r.kr==kr, rows)
            (isempty(b) || isempty(u)) && continue
            @printf("  k_sac=%-8.0e : bilatéral %3.0f %%  →  unilatéral %3.0f %%   (décollement : %3.0f %% des configs)\n",
                    kr, 100*count(r->r.Δ>0,b)/length(b), 100*count(r->r.Δ>0,u)/length(u),
                    100*count(r->r.f_lift>0.01,u)/length(u))
        end

        println("\nMEILLEURE configuration (Δ maximal, contact unilatéral) :")
        @printf("  culasse=%.0f in  avant-bras=%.0f in  EI=%.0e  h_âme=%.1f in  k_sac=%.0e\n",
                best.xb*M2IN, best.lf*M2IN, best.eis, best.hb*M2IN, best.kr)
        @printf("  proj statique = %+.3f in (Harral −1.32)   f_whip=%.0f Hz   décollement %.0f %% du temps\n",
                best.proj, best.f1, 100*best.f_lift)
        @printf("  spread : %s\n", join([@sprintf("%.4f", s) for s in best.sp], "  "))
        @printf("  Harral : %s\n", join([@sprintf("%.4f", s) for s in HARRAL], "  "))

        # ---------------------------------------------------------------------
        # FILTRE PHYSIQUE — et ce n'est PAS du curve-fitting.
        # Harral publie DEUX colonnes indépendantes : `proj` (flèche statique de
        # bouche, −1.32 in) et `spread` (dispersion). On se sert de la PREMIÈRE
        # pour éliminer les configurations absurdes — le balayage produit des
        # crosses qui fléchissent de 9 pouces ! — puis on TESTE la seconde, que
        # rien n'a calée. Caler sur `spread` serait tricher ; caler sur `proj`
        # et prédire `spread`, c'est la démarche normale.
        # ---------------------------------------------------------------------
        adm = filter(r -> abs(r.proj - (-1.32)) <= 0.25, rows)
        println("\n" * "-"^78)
        println("CONFIGS PHYSIQUEMENT ADMISSIBLES : flèche statique = −1.32 ± 0.25 in")
        println("(contrainte tirée de la colonne `proj` de Harral — INDÉPENDANTE de `spread`)")
        if isempty(adm)
            println("  aucune — les plages balayées ne contiennent pas l'arme de Harral.")
        else
            @printf("  %d / %d configurations retenues\n", length(adm), length(rows))
            @printf("  tendance CROISSANTE : %d / %d  (%.0f %%)\n",
                    count(r->r.Δ>0, adm), length(adm), 100*count(r->r.Δ>0, adm)/length(adm))
            @printf("  STRICTEMENT MONOTONE : %d / %d  (%.0f %%)\n",
                    count(r->r.mono, adm), length(adm), 100*count(r->r.mono, adm)/length(adm))
            sp1 = [r.sp[1] for r in adm]; sp4 = [r.sp[end] for r in adm]
            @printf("  spread canon nu  : %.3f … %.3f in   (Harral : %.3f)\n",
                    minimum(sp1), maximum(sp1), HARRAL[1])
            @printf("  spread + 16.0 oz : %.3f … %.3f in   (Harral : %.3f)\n",
                    minimum(sp4), maximum(sp4), HARRAL[end])
            println("  EI de crosse retenus : ",
                    join(sort(unique([@sprintf("%.0e", r.eis) for r in adm])), ", "))
        end

        println("\n" * "-"^78)
        println("Sensibilité du signe de Δ à chaque inconnue (contact unilatéral) :")
        for (nm, get, vals) in (("culasse (in)",  r->r.xb*M2IN, X_BREECH .* M2IN),
                                ("av.-bras (in)", r->r.lf*M2IN, L_FORE .* M2IN),
                                ("EI crosse",     r->r.eis,     EI_STOCK),
                                ("h_âme (in)",    r->r.hb*M2IN, H_BORE .* M2IN),
                                ("k_sac (N/m)",   r->r.kr,      K_REST))
            parts = String[]
            for v in vals
                sel = filter(r -> isapprox(get(r), v; rtol=1e-6), rows)
                isempty(sel) && continue
                push!(parts, @sprintf("%.6g→%3.0f%%", v,
                                      100*count(r->r.Δ>0, sel)/length(sel)))
            end
            @printf("  %-14s : %s\n", nm, join(parts, "   "))
        end
    end
    println("="^78)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_sweep()
end
