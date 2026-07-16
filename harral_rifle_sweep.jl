# =============================================================================
# Étape 2 du chantier « modèle de Harral » : LE FUSIL ENTIER SUR SES DEUX SACS
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
# On les BALAIE donc sur des plages plausibles, et on ne regarde qu'une chose :
# LA STRUCTURE (spread croissante avec la masse) tient-elle sur toute la plage ?
#   - sur toute la plage  → mécanisme confirmé, indépendant des valeurs choisies
#   - dans un coin étroit → signal d'alarme, à dire
#   - jamais              → le diagnostic « il manque la rotation » est FAUX
# Le succès NE se mesure PAS à l'écart aux chiffres de Harral.
#
# Ancrage réel (gratuit) : le budget de masse. Canon = 1.711 kg ; fusil total
# = 10.5 lb = 4.763 kg ⇒ boîte + crosse = 3.052 kg (64 % de l'arme).
#
# CHOIX DE MODÉLISATION (et leurs limites) :
#  * Poutre plane unique, talon → bouche. Harral maille en briques solides 3D.
#  * L'avant-bras est un membre PARALLÈLE au canon (canon flottant), pas
#    modélisable dans une poutre unique : on reporte donc le sac avant sur la
#    crosse, en arrière de la culasse. Le canon reste en PORTE-À-FAUX au-delà
#    de l'appui avant — c'est cela qui produit l'affaissement de bouche.
#  * Appuis = ressorts verticaux (sacs de sable), BILATÉRAUX. Justification :
#    un sac est MOU. À k = 1e4..1e6 N/m sur 4.76 kg, la fréquence d'appui est
#    de ~10 à ~100 Hz, soit une période très supérieure aux 2.4 ms de séjour de
#    la balle : le sac n'a pas le temps de réagir, et bilatéral ≈ unilatéral
#    sur cette fenêtre. (k est balayé : s'il s'avère déterminant, il faudra le
#    contact unilatéral — le talon décolle — et ce sera à traiter.)
#  * Le moment de recul est enfin ANCRÉ PHYSIQUEMENT : la force axiale
#    p(t)·A_âme s'applique à la culasse, à la hauteur h_bore au-dessus de la
#    ligne des sacs ⇒ moment de tangage p·A·h_bore. Avec la culasse encastrée
#    ce bras de levier n'avait aucune référence et flottait librement.
#  * Crosse + boîte : masse uniforme (l'acier de la boîte est en réalité bien
#    plus dense que le bois du talon) et EI uniforme. Simplification assumée.
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
const M_STOCK       = M_RIFLE_TOTAL - M_BARREL   # boîte + crosse

# -----------------------------------------------------------------------------
# Maillage du fusil entier : x = 0 au talon, culasse en x_breech, bouche au bout
# -----------------------------------------------------------------------------
const N_STOCK = 16
const N_BAR   = N_elements                        # 48, comme le canon nu

function rifle_mesh(x_breech)
    xs_stock = collect(range(0.0, x_breech, length = N_STOCK + 1))
    xs_bar   = collect(range(x_breech, x_breech + L, length = N_BAR + 1))
    return vcat(xs_stock, xs_bar[2:end])          # nœuds, m
end

nearest_node(xs, x) = argmin(abs.(xs .- x))

# Propriétés d'un élément : (EI, ρA)
function rifle_element(e, xs, x_breech, EI_stock)
    if e <= N_STOCK
        return EI_stock, M_STOCK / x_breech
    else
        x_mid = 0.5 * (xs[e] + xs[e+1]) - x_breech     # abscisse LOCALE canon
        D = D_out_of_x(x_mid)
        A = π/4  * (D^2 - D_bore^2)
        I = π/64 * (D^4 - D_bore^4)
        return E * I, ρ_steel * A
    end
end

# -----------------------------------------------------------------------------
# Assemblage : K, M, et le vecteur de gravité (aucun d.d.l. supprimé — l'arme
# n'est PAS encastrée ; ce sont les ressorts d'appui qui la tiennent)
# -----------------------------------------------------------------------------
function build_rifle(; m_tuner, x_breech, x_front, x_rear, EI_stock, k_rest)
    xs = rifle_mesh(x_breech)
    nn = length(xs)
    nd = 2 * nn
    K = zeros(nd, nd); M = zeros(nd, nd); Fg = zeros(nd)

    for e in 1:(nn - 1)
        le = xs[e+1] - xs[e]
        EI, ρA = rifle_element(e, xs, x_breech, EI_stock)
        Ke, Me = element_matrices(le, EI, ρA)
        idx = (2*e - 1):(2*e + 2)
        @views K[idx, idx] .+= Ke
        @views M[idx, idx] .+= Me
        q  = -ρA * g_ms                                  # N/m, vers le bas
        Fe = q * le * [0.5, le/12, 0.5, -le/12]
        @views Fg[idx] .+= Fe
    end

    # Tuner à la bouche (dernier nœud)
    M[nd-1, nd-1] += m_tuner
    M[nd,   nd]   += tuner_inertia(m_tuner)
    Fg[nd-1]      += -m_tuner * g_ms

    # Sacs de sable : ressorts VERTICAUX seuls (contact sans frottement)
    i_rear  = nearest_node(xs, x_rear)
    i_front = nearest_node(xs, x_front)
    K[2*i_rear  - 1, 2*i_rear  - 1] += k_rest
    K[2*i_front - 1, 2*i_front - 1] += k_rest

    i_breech = N_STOCK + 1
    return (K=K, M=M, Fg=Fg, xs=xs, nd=nd, i_breech=i_breech,
            i_rear=i_rear, i_front=i_front)
end

# Charge ponctuelle cohérente sur un maillage quelconque
function point_load_rifle(Fval, x, xs, nd)
    (x <= xs[1] || x >= xs[end]) && return zeros(nd)
    e = searchsortedlast(xs, x); e = clamp(e, 1, length(xs) - 1)
    le = xs[e+1] - xs[e]
    ξ  = (x - xs[e]) / le
    Ns = (1 - 3ξ^2 + 2ξ^3, le*(ξ - 2ξ^2 + ξ^3), 3ξ^2 - 2ξ^3, le*(-ξ^2 + ξ^3))
    F = zeros(nd)
    for (k, ig) in enumerate((2*e-1, 2*e, 2*e+1, 2*e+2))
        F[ig] += Fval * Ns[k]
    end
    return F
end

# -----------------------------------------------------------------------------
# Un tir sur le fusil entier
# -----------------------------------------------------------------------------
function shoot_rifle(v_muzzle; m_tuner, x_breech, x_front, x_rear, EI_stock,
                     k_rest, h_bore, ζ1 = 0.01, ζ2 = 0.06,
                     Δt = 2e-6, t_end = 3.0e-3)
    S = build_rifle(; m_tuner, x_breech, x_front, x_rear, EI_stock, k_rest)
    _, ωs = modal_analysis(S.K, S.M; n_modes = 3)
    C  = rayleigh_damping(S.M, S.K, ωs[1], ωs[2], ζ1, ζ2)
    xp = projectile_pos(v_muzzle)

    U_stat = S.K \ S.Fg                       # équilibre statique sur les sacs

    function F_of_t(t)
        F = copy(S.Fg)
        # Moment de tangage : force axiale × hauteur d'âme au-dessus des sacs
        F[2 * S.i_breech] += chamber_pressure(t) * A_bore * h_bore
        xg = x_breech + xp(t)                 # position ABSOLUE du projectile
        if x_breech < xg < x_breech + L
            F .+= point_load_rifle(-m_p * g_ms, xg, S.xs, S.nd)
        end
        return F
    end

    ts, U, V = newmark_solve(S.M, C, S.K, F_of_t, t_end, Δt; U0 = U_stat)

    t_b = exit_time(v_muzzle)
    ib  = argmin(abs.(ts .- t_b))
    return (θ_tb = U[S.nd, ib], ẏ_tb = V[S.nd - 1, ib],
            θ_stat = U_stat[S.nd], y_stat = U_stat[S.nd - 1],
            f1 = ωs[1] / (2π))
end

function spread_for(m_tuner; kw...)
    r_lo = shoot_rifle(V_LO; m_tuner = m_tuner, kw...)
    r_hi = shoot_rifle(V_HI; m_tuner = m_tuner, kw...)
    poi(r, vfps) = r.θ_tb * D_TARGET_IN + (r.ẏ_tb * M2IN) * TOF(DROP[vfps]) - DROP[vfps]
    return abs(poi(r_lo, 1035) - poi(r_hi, 1075)), r_lo
end

# =============================================================================
# BALAYAGE
# =============================================================================
const TUNERS = [0.0, 4.9*OZ, 8.6*OZ, 16.0*OZ]
const HARRAL = [0.0910, 0.1218, 0.1554, 0.1832]

# Plages plausibles (pratique benchrest) — AUCUNE n'est publiée par Harral.
const X_BREECH  = [18.0, 20.0, 22.0] .* IN      # culasse depuis le talon
const F_FRONT   = [0.65, 0.80, 0.92]            # sac avant, en fraction de x_breech
const EI_STOCK  = [3.0e3, 1.0e4, 4.0e4]         # raideur crosse+boîte (N·m²)
const H_BORE    = [1.0, 2.0, 3.0] .* IN         # âme au-dessus de la ligne des sacs
const K_REST    = [1.0e4, 1.0e5, 1.0e6]         # raideur d'un sac (N/m)
const X_REAR    = 1.5 * IN                      # sac arrière, près du talon

function main_sweep()
    println("="^78)
    println(" Étape 2 — fusil entier sur deux sacs : LA STRUCTURE APPARAÎT-ELLE ?")
    println("="^78)
    @printf("Budget de masse : canon %.3f kg + (boîte+crosse) %.3f kg = %.3f kg (%.1f lb)\n",
            M_BARREL, M_STOCK, M_RIFLE_TOTAL, M_RIFLE_TOTAL/0.45359237)
    @printf("Harral (référence) : spread %.3f → %.3f in (CROISSANT)\n", HARRAL[1], HARRAL[end])
    @printf("Culasse encastrée (étape 1) : ≈ 0.171 → 0.164 in (plat, non compensé)\n\n")

    n_ok = 0; n_mono = 0; n_tot = 0
    best = nothing; worst = nothing
    rows = NamedTuple[]

    for xb in X_BREECH, ff in F_FRONT, eis in EI_STOCK, hb in H_BORE, kr in K_REST
        kw = (x_breech = xb, x_front = ff*xb, x_rear = X_REAR,
              EI_stock = eis, k_rest = kr, h_bore = hb)
        sp = Float64[]
        local r0
        ok = true
        try
            for m in TUNERS
                s, r = spread_for(m; kw...)
                push!(sp, s)
                m == 0.0 && (r0 = r)
            end
        catch err
            ok = false
        end
        ok || continue

        n_tot += 1
        Δ    = sp[end] - sp[1]                      # > 0 ⇒ tendance de Harral
        mono = all(diff(sp) .> 0)
        Δ > 0 && (n_ok += 1)
        mono && (n_mono += 1)
        push!(rows, (xb=xb, ff=ff, eis=eis, hb=hb, kr=kr, sp=copy(sp), Δ=Δ,
                     mono=mono, proj=r0.θ_stat*D_TARGET_IN, f1=r0.f1))
        if best === nothing || Δ > best.Δ;  best  = rows[end]; end
        if worst === nothing || Δ < worst.Δ; worst = rows[end]; end
    end

    println("-"^78)
    @printf("Configurations balayées : %d\n", n_tot)
    @printf("  tendance CROISSANTE (spread(16oz) > spread(nu)) : %d / %d  (%.0f %%)\n",
            n_ok, n_tot, 100*n_ok/n_tot)
    @printf("  croissance STRICTEMENT MONOTONE sur les 4 masses : %d / %d  (%.0f %%)\n",
            n_mono, n_tot, 100*n_mono/n_tot)
    println("-"^78)

    for (lbl, r) in (("Δ MAXIMAL", best), ("Δ MINIMAL", worst))
        println("\n$lbl :")
        @printf("  x_breech=%.1f in  sac_avant=%.2f·x_b  EI_crosse=%.0e  h_âme=%.1f in  k_sac=%.0e\n",
                r.xb*M2IN, r.ff, r.eis, r.hb*M2IN, r.kr)
        @printf("  proj statique bouche = %+.3f in   (Harral : −1.32 in)   f₁=%.1f Hz\n",
                r.proj, r.f1)
        @printf("  spread : %s\n", join([@sprintf("%.4f", s) for s in r.sp], "  "))
        @printf("  Harral : %s\n", join([@sprintf("%.4f", s) for s in HARRAL], "  "))
        @printf("  Δ(16oz − nu) = %+.4f in   %s\n", r.Δ,
                r.Δ > 0 ? "(sens de Harral)" : "(sens OPPOSÉ)")
    end

    # Sensibilité : quel paramètre gouverne le signe de Δ ?
    println("\n" * "-"^78)
    println("Sensibilité du signe de Δ à chaque inconnue (fraction de Δ>0) :")
    for (nm, get, vals) in (("x_breech (in)", r->r.xb*M2IN, X_BREECH .* M2IN),
                            ("sac_avant/x_b", r->r.ff,      F_FRONT),
                            ("EI_crosse",     r->r.eis,     EI_STOCK),
                            ("h_âme (in)",    r->r.hb*M2IN, H_BORE .* M2IN),
                            ("k_sac (N/m)",   r->r.kr,      K_REST))
        parts = String[]
        for v in vals
            sel = filter(r -> isapprox(get(r), v; rtol=1e-6), rows)
            isempty(sel) && continue
            frac = 100 * count(r -> r.Δ > 0, sel) / length(sel)
            push!(parts, @sprintf("%.6g→%3.0f%%", v, frac))
        end
        @printf("  %-14s : %s\n", nm, join(parts, "   "))
    end
    println("="^78)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_sweep()
end
