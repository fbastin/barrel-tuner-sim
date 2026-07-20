# =============================================================================
# Simulation des vibrations transversales d'un canon de carabine
# Éléments finis (Euler-Bernoulli + Hermite cubique) + Newmark-β
# Application : analyse de la compensation positive et accord par tuner
#
# Conventions :
#   - 2 d.d.l. par nœud : déflexion transverse y et rotation θ = ∂y/∂x
#   - Encastrement à la culasse (nœud 1), bouche libre (nœud N+1)
#   - Une masse ponctuelle (tuner) est ajoutée à la bouche
#   - Excitation : (a) moment à la culasse via pression de chambre × bras de
#     levier h_offset, (b) poids du projectile traité comme charge mobile
#   - Le bras de levier h_offset est ancré physiquement sur deux grandeurs
#     mesurables de l'arme (poids m_rifle, hauteur âme↔CG h_cg) : amplitude
#     ∝ h_cg et ∝ 1/m_rifle (cf. physical_h_offset)
#   - Sortie : θ(L,t), y(L,t), θ̇(L,t), évalués à t_b
#
# Usage :   julia simulation.jl
# =============================================================================

using LinearAlgebra
using Printf
using Pkg

# Auto-installation de Plots.jl si absent (utilisé en fin de script)
const PLOTS_AVAILABLE = try
    using Plots
    true
catch
    println("→ Installation de Plots.jl (1ʳᵉ exécution)…")
    Pkg.add("Plots")
    using Plots
    true
end

# -----------------------------------------------------------------------------
# 1. PARAMÈTRES PHYSIQUES (SI)
# -----------------------------------------------------------------------------
const L         = 0.66          # Longueur canon (26 in)
const D_out     = 0.024         # Diamètre extérieur (profil match)
const D_in      = 0.0056        # Diamètre intérieur (.22 LR)
const E         = 200e9         # Module d'Young acier
const ρ_steel   = 7850.0        # Masse volumique acier

# Tuner (Kolbe : 200 g à la bouche)
const m_tuner_0 = 0.200

# Inertie PROPRE du tuner autour de son centre de masse. L'architecture retenue
# est un ENSEMBLE TUBE façon Starik/Centra — l'objet étalé (tubes de 19, 32 ou
# 36 cm) dans lequel une masse coulisse — et ce n'est pas un choix neutre :
# c'est la SEULE qui autorise le porte-à-faux de ~10 cm auquel aboutit l'accord
# au nœud (§6 du wiki), un tuner à corps vissé (Harrell's, Ezell, PMA) n'offrant
# que quelques millimètres de filetage. Une bague compacte donnerait k ≈ 1,7 cm
# à 200 g et un modèle qui ne compense plus du tout en masse pure — cohérent,
# mais décrivant un autre produit que celui dont il est question.
#
# POURQUOI UNE FONCTION ET NON UNE CONSTANTE. Jusqu'au 2026-07-18 J_tuner_0
# valait 5,0e-4, FIXE — donc indépendant de la masse. Le balayage en masse
# faisait varier m de 0 à 400 g en gardant J = 5,0e-4, si bien que le « canon
# nu » (m = 0) portait encore 5,0e-4 kg·m² d'inertie de rotation à la bouche.
# Le canon nu n'était pas nu.
#
# k VARIABLE AVEC LA MASSE (2026-07-19). k = 5 cm FIXE restait une approximation
# fautive dans un balayage en masse : elle prêtait au tuner de 25 g l'étalement
# d'un tube de 17 cm, et à celui de 400 g un tube deux fois trop court. J est
# désormais DÉRIVÉ de la géométrie, comme dans harral_a22lr.jl — même formule
# d'anneau, section de TUBE (paroi mince) au lieu d'une bague pleine :
#
#     J = m·(3(Ro² + Ri²) + ℓ²)/12,   ℓ = m/(ρ·A)     [autour du CdM du tuner]
#
# La section Ø40 × 1,25 mm n'est pas libre : elle est choisie pour redonner
# k = 5,02 cm à 200 g, soit la valeur nominale antérieure — le cas de référence
# est donc préservé, seule la DÉPENDANCE en masse change. Contrôle de réalisme
# non imposé par ce calage : les longueurs dérivées (4,2 cm à 50 g → 33,5 cm à
# 400 g) couvrent exactement la gamme des tubes Starik réels (19, 32, 36 cm).
#
# CONVENTION ET SA LIMITE. Faire croître la LONGUEUR avec la masse est la
# convention de harral_a22lr.jl, retenue ici pour que les deux modèles restent
# structurellement comparables. Sur un tube réel on n'allonge pas le tube : on
# fait COULISSER une masse dans un tube de longueur fixe. Les deux coïncident
# tant que la masse ajoutée reste répartie, et divergent pour un poids compact
# placé loin du centre. Un modèle fidèle demanderait deux corps (tube fixe +
# masse mobile) et coupleraient J à d_overhang, aujourd'hui réglage indépendant
# — refonte non faite, et non requise par les questions traitées.
#
# DIVERGENCE D'ARCHITECTURE ASSUMÉE AVEC LA FAMILLE harral_*. Les deux modèles
# ne décrivent PAS le même produit, et leurs inerties ne sont pas comparables.
# Depuis le 2026-07-19 tous deux DÉRIVENT J de la géométrie (même formule
# d'anneau), mais sur des sections différentes :
#
#   ici (simulation.jl)   TUBE paroi mince   Ø40 × 1,25 mm
#                         → k = 1,8 cm (50 g) à 9,8 cm (400 g)
#                         → k = 5,02 cm et J = 5,05e-4 kg·m² à 200 g
#   harral_a22lr.jl       BAGUE pleine       OD 1,4" / ID 0,915"
#                         → k = 1,1 cm (50 g) à 2,8 cm (400 g)
#                         → k = 1,67 cm et J = 5,6e-5 kg·m² à 200 g
#
# Soit un facteur ~9 sur J à 200 g. Ce n'est pas une incohérence à résoudre mais
# deux points de la gamme réelle : une bague vissée courte d'un côté, un tube
# type Starik/Centra de l'autre. Conséquence pratique : ne PAS confronter
# directement les θ̇ des deux familles ; leurs résultats ne se recoupent que sur
# les grandeurs cinématiques (τ_v, taux requis), qui ne dépendent d'aucune
# inertie de tuner.
const TUBE_OD   = 0.040                       # diamètre extérieur du tube (m)
const TUBE_WALL = 0.00125                     # épaisseur de paroi (m)
const TUBE_RO   = TUBE_OD / 2
const TUBE_RI   = TUBE_RO - TUBE_WALL
const TUBE_A    = π * (TUBE_RO^2 - TUBE_RI^2)

tuner_length(m_tuner) = m_tuner / (ρ_steel * TUBE_A)

function tuner_inertia(m_tuner)
    m_tuner <= 0 && return 0.0                # canon nu : aucune inertie ajoutée
    ℓ = tuner_length(m_tuner)
    return m_tuner * (3 * (TUBE_RO^2 + TUBE_RI^2) + ℓ^2) / 12
end

const J_tuner_0 = tuner_inertia(m_tuner_0)

# Projectile et balistique interne (.22 LR Match, Eley Tenex)
const m_p       = 2.6e-3        # 40 grains
const v_muzzle  = 318.0         # 1043 ft/s (Eley Tenex : cohérent avec le τ_v de Kolbe)
# NI t_b NI LA PRESSION NE SONT PLUS POSÉS. Les deux découlent de la balistique
# intérieure couplée (section 6). Les anciennes constantes p_max / t_peak /
# α_press / PHI_BURN ont été retirées le 2026-07-19 : elles décrivaient un profil
# autonome dont l'impulsion valait 5,19× le recul physique, incompatible avec la
# cinématique posée à côté. Le profil fautif survit dans pcp_vs_firearm.jl, qui
# documente le défaut, et nulle part ailleurs.
const h_offset_default = 0.005  # Bras de levier (calibrable)

# Arme complète : grandeurs mesurables pilotant l'AMPLITUDE des vibrations.
# Le moment de recul vaut M0(t) = p(t)·A_bore·h_cg ; l'amplitude angulaire
# résultante décroît comme 1/m_rifle (une arme lourde recule moins). Les valeurs
# de référence sont celles par défaut de l'outil de Kolbe.
const m_rifle_ref = 5.0         # Masse totale de l'arme de référence (kg)
const h_cg_ref    = 0.0254      # Hauteur âme ↔ centre de gravité de réf. (1 in)

# Cible balistique et constantes
const D_target  = 50.0
const g_accel   = 9.81
const θdot_optimum_MOAms = 6.0  # Kolbe : optimum à 50 m

# -----------------------------------------------------------------------------
# 2. MAILLAGE
# -----------------------------------------------------------------------------
const N_elements = 20
const N_nodes    = N_elements + 1
const ndof       = 2 * N_nodes
const L_e        = L / N_elements

const A_sec  = π/4  * (D_out^2 - D_in^2)
const I_sec  = π/64 * (D_out^4 - D_in^4)
const A_bore = π/4  * D_in^2
const EI     = E * I_sec

# -----------------------------------------------------------------------------
# ÉCHAUFFEMENT DU CANON (ajouté le 2026-07-19)
#
# Dai, Fu, Cao, Lyu et Xu (Mechanics of Solids, 2024) mesurent sur fusil
# automatique une chute de la fréquence propre entre canon froid et canon chaud :
# −12,7 % en calcul (Timoshenko), −26,4 % en mesure. La cause est la baisse du
# module d'Young de l'acier avec la température.
#
# POURQUOI NE PAS APPLIQUER LEUR CHIFFRE TEL QUEL. Leur « canon chaud » est
# celui d'une arme automatique après rafales ; une carabine de match sur une
# série de 25 coups ne monte que de 30 à 60 °C. On paramètre donc par la
# PHYSIQUE plutôt que par leur écart : E décroît d'environ 3,6e-4 par degré
# (≈ −25 % à 700 °C), et f ∝ √E. Contrôle : ce coefficient reproduit leurs
# −12,7 % pour ΔT ≈ 660 °C, valeur plausible pour leur essai.
#
# ORDRE DE GRANDEUR POUR NOUS. +30 à +60 °C ⇒ chute de fréquence de 0,5 à 1,1 %,
# soit 5 à 11 % de la course qu'offre un tuner de 200 g (~10 %). L'effet est donc
# réel mais modeste en tir de match — et déterminant en tir soutenu. Il pourrait
# expliquer la dérive d'un accord sur une longue série et la singularité du coup
# à canon froid, que rien dans ce modèle ne reproduisait jusqu'ici.
# AMORTISSEMENT (révisé le 2026-07-19 d'après le dépouillement de la littérature)
#
# Les valeurs antérieures — ζ₁ = 0,5 %, ζ₂ = 1 % — étaient posées sans référence.
# Le dépouillement a produit quatre estimations publiées, et elles couvrent plus
# d'un ORDRE DE GRANDEUR :
#
#   Sava (2015), 5,56 automatique, non bridé ......... 0,7 à 1,4 %
#   Étude fusil de précision, DANS UN ÉTAU ........... 1,7 %
#   Benet ARCCB, canon de 25 mm, mesure directe ...... 3,6 à 4,6 %
#   Dai (2024), fusil automatique .................... ~8 à 10 %
#
# On retient **1 %**, milieu de la fourchette de Sava — seule mesure sur arme de
# type carabine dans une configuration non bridée, donc la plus proche de notre
# cas. Les valeurs plus élevées concernent un canon de 25 mm, une arme dans un
# étau (qui ajoute ses pertes de contact) ou une arme automatique.
#
# ⚠️ LES RÉSULTATS DE CE MODÈLE EN DÉPENDENT CRITIQUEMENT, et c'est à dire :
#
#   ζ = 1 %  → θ̇ atteint +5,98, la cible de compensation est accessible
#   ζ = 2 %  → θ̇ plafonne à +3,94, la cible devient INATTEIGNABLE
#   ζ = 5 %  → θ̇ plafonne à +2,26 et le canon nu passe POSITIF : la structure
#              de signe mesurée par Kolbe s'inverse
#
# Autrement dit, tout ce que ce modèle produit d'exploitable vit sur le BAS de
# la fourchette publiée. Si l'amortissement réel d'une carabine est de 2 % ou
# plus, ce modèle ne décrit pas la compensation positive. Deux des quatre
# sources ci-dessus sont au-dessus de ce seuil. C'est la réserve la plus lourde
# de tout le modèle, et elle ne se lèvera que par une mesure (cf. la page wiki
# « Mesurer l'amortissement d'un canon »).
const ZETA_DEFAULT = 0.01

const K_E_TEMP = 3.6e-4                       # baisse relative de E par °C
young_at(ΔT) = E * (1 - K_E_TEMP * ΔT)
EI_at(ΔT) = young_at(ΔT) * I_sec
const ρA     = ρ_steel * A_sec

# Conversions
const MOA_per_rad = 1 / (π / (180 * 60))
moa_per_ms(rad_per_s) = rad_per_s * MOA_per_rad * 1e-3

# -----------------------------------------------------------------------------
# 3. MATRICES ÉLÉMENTAIRES (poutre de Bernoulli, Hermite cubique)
# -----------------------------------------------------------------------------
function element_matrices(L_e, EI, ρA)
    Ke = (EI / L_e^3) * [
        12.0     6*L_e       -12.0    6*L_e   ;
        6*L_e    4*L_e^2     -6*L_e   2*L_e^2 ;
       -12.0    -6*L_e        12.0   -6*L_e   ;
        6*L_e    2*L_e^2     -6*L_e   4*L_e^2
    ]
    Me = (ρA * L_e / 420.0) * [
        156.0    22*L_e       54.0   -13*L_e   ;
        22*L_e   4*L_e^2      13*L_e  -3*L_e^2 ;
        54.0     13*L_e      156.0   -22*L_e   ;
       -13*L_e  -3*L_e^2     -22*L_e   4*L_e^2
    ]
    return Ke, Me
end

# -----------------------------------------------------------------------------
# 4. ASSEMBLAGE + ENCASTREMENT + TUNER
# -----------------------------------------------------------------------------
# d_overhang : porte-à-faux du centre de masse du tuner DEVANT la bouche (m).
# C'est le vrai réglage de terrain — on visse le tuner plus ou moins loin, la
# masse restant fixe. Le couplage masse/rotation d'une masse déportée s'écrit
#     M_add = [m      m·d    ;
#              m·d    m·d² + J]
# sur le couple (y, θ) du nœud de bouche : à d = 0 on retrouve le cas classique
# de la masse ponctuelle.
# ENCASTREMENT ÉLASTIQUE (ajouté le 2026-07-19). `K_root = Inf` conserve
# l'encastrement rigide, comportement par défaut et inchangé.
#
# MOTIF. O'Neil (2022) mesure sur Tikka T3X un rapport f₂/f₁ = 8,16, quand tous
# les modèles — le nôtre, son ANSYS, son analytique — donnent 6,27, valeur
# théorique de la poutre encastrée-libre. Un canon réel, vissé dans une boîte
# posée dans une crosse, n'est pas parfaitement encastré. Un ressort de torsion
# de 1,09e4 N·m/rad à la culasse reproduit le rapport mesuré.
#
# POURQUOI CE N'EST PAS LE DÉFAUT. Les fréquences donnant 8,16 viennent d'un
# essai que l'auteur juge lui-même insuffisamment résolu, et le pic à 263,7 Hz
# pourrait appartenir à la crosse plutôt qu'au canon. On offre la capacité, on
# ne la présume pas.
function build_system(m_tuner, J_tuner; d_overhang = 0.0, ΔT = 0.0, K_root = Inf)
    Ke, Me = element_matrices(L_e, EI_at(ΔT), ρA)
    K = zeros(ndof, ndof)
    M = zeros(ndof, ndof)
    for e in 1:N_elements
        idx = (2*e - 1):(2*e + 2)
        @views K[idx, idx] .+= Ke
        @views M[idx, idx] .+= Me
    end
    if isfinite(K_root)
        K[2, 2] += K_root          # ressort de torsion à la culasse
        active = 2:ndof            # θ₁ libéré, y₁ toujours bloqué
    else
        active = 3:ndof            # encastrement rigide (défaut)
    end
    Ka = K[active, active]
    Ma = copy(M[active, active])
    d = d_overhang
    Ma[end-1, end-1] += m_tuner
    Ma[end-1, end  ] += m_tuner * d
    Ma[end,   end-1] += m_tuner * d
    Ma[end,   end  ] += m_tuner * d^2 + J_tuner
    return Ka, Ma
end

# -----------------------------------------------------------------------------
# 5. ANALYSE MODALE
# -----------------------------------------------------------------------------
function modal_analysis(Ka, Ma; n_modes = 5)
    decomp = eigen(Ka, Ma)
    λ = real.(decomp.values)
    Φ = real.(decomp.vectors)
    keep = λ .> 1e-3
    λ_ok = λ[keep]
    Φ_ok = Φ[:, keep]
    order = sortperm(λ_ok)
    λ_ok = λ_ok[order]
    Φ_ok = Φ_ok[:, order]
    ω = sqrt.(λ_ok)
    f = ω ./ (2π)
    for k in 1:size(Φ_ok, 2)
        Φ_ok[:, k] ./= sqrt(Φ_ok[:, k]' * Ma * Φ_ok[:, k])
    end
    nm = min(n_modes, length(f))
    return f[1:nm], ω[1:nm], Φ_ok[:, 1:nm]
end

# -----------------------------------------------------------------------------
# 6. BALISTIQUE INTÉRIEURE COUPLÉE
#
# REFONTE DU 2026-07-19. Jusqu'ici le modèle portait DEUX représentations
# incompatibles du même coup de feu : une cinématique « burnout » posée a priori
# (accélération constante sur φ=0,35 du canon, calée pour reproduire le
# τ_v = 8,8 µs/(m/s) de Kolbe) et, à côté, un profil de pression servant à
# l'excitation dont l'impulsion valait 5,19× le recul physique. La première
# était juste, la seconde fausse, et rien ne les reliait.
#
# Désormais la cinématique est INTÉGRÉE depuis la pression :
#
#     m_eff(t)·ẍ = p(t,x)·A_bore ,     p = f·ω·z(t) / (V₀ + A_bore·x)
#
# avec z(t) la fraction de poudre brûlée. La conservation de la quantité de
# mouvement devient alors EXACTE PAR CONSTRUCTION — ∫p·A dt = m_eff·v_bouche —
# et le défaut corrigé ici ne peut plus réapparaître par dérive d'un paramètre.
#
# CE QUE CELA CHANGE DE STATUT. τ_v cesse d'être un paramètre ajusté pour
# devenir une PRÉDICTION du modèle, donc un test. Prédit : 8,3 µs/(m/s) contre
# 8,8 mesuré par Kolbe, soit −6 % SANS calage — là où l'ancien φ=0,35 était
# choisi pour tomber juste. On perd 6 % d'accord nominal et on gagne un test.
#
# TENSION NON RÉSOLUE. Aucune paramétrisation essayée ne réconcilie le τ_v
# mesuré avec le pic de pression SAAMI de la .22 LR (165 MPa) : viser 8,3
# pousse le pic à ~34 MPa. Trois formulations donnent le même compromis (loi
# de puissance en temps, détente adiabatique en volume, combustion progressive).
# On retient ici le τ_v, qui est mesuré et qui gouverne la compensation ; le pic
# est une conséquence, et il est bas. À reprendre avec un vrai modèle de
# combustion si la question devient critique.
# -----------------------------------------------------------------------------
# PARAMÈTRES CONTRAINTS PAR LA LITTÉRATURE (2026-07-19), non ajustés librement.
# Source : G. Kolbe, « The Ballistics Handbook » (2000), données .22 LR reprises
# par R. Kenchington (Long Range Rimfire Club) — même auteur que le τ_v de
# référence, donc cohérent avec le reste de la validation :
#   • pic de pression à 0,25 ms, après 9,4 mm de trajet
#   • poudre entièrement brûlée vers 1 pouce de trajet
#   • vitesse maximale à 19 pouces (~1,5 ms), sortie du 28" à 2,3 ms
#   • pic ~15 000 psi = 103 MPa  (et NON les 165 MPa de la limite SAAMI, qui
#     est un maximum admissible et non une valeur de travail — erreur commise
#     lors du premier calage)
# τ_b = 150 µs place la fin de combustion au bon endroit ; le modèle rend alors
# un pic de 99 MPa (−3 % sur la cible), et surtout retrouve la STRUCTURE DE
# SIGNE de Kolbe : canon nu négatif, accord ramenant positif. À τ_b = 600 µs,
# valeur choisie arbitrairement lors de la refonte, cette structure était perdue
# (canon nu déjà au-dessus de la cible, tuner purement dégradant).
const V0_CHAMBER = 1.0e-7       # volume de chambre efficace (m³)
const TAU_BURN   = 150e-6       # temps caractéristique de combustion (s)
const N_BURN     = 2.0          # exposant de la loi de combustion
const M_POWDER   = 0.10e-3      # charge de poudre (kg)
const DT_IB      = 2e-8         # pas d'intégration de la balistique intérieure

burnt_fraction(t) = t <= 0 ? 0.0 : 1 - exp(-(t / TAU_BURN)^N_BURN)

# Intègre la trajectoire jusqu'à la bouche. `fω` est l'impétus total (J), seul
# paramètre d'échelle : il est calé une fois pour rendre v(L) = v_muzzle.
function integrate_bore(fω, v_target, L; dt = DT_IB, tmax = 20e-3)
    x, v, t = 0.0, 0.0, 0.0
    ts, xs, vs, ps = Float64[], Float64[], Float64[], Float64[]
    while t < tmax
        z    = burnt_fraction(t)
        pr   = fω * z / (V0_CHAMBER + A_bore * x)
        m_ef = m_p + M_POWDER * z / 3          # Lagrange : gaz entraîné
        push!(ts, t); push!(xs, x); push!(vs, v); push!(ps, pr)
        v += pr * A_bore / m_ef * dt
        x += v * dt
        t += dt
        if x >= L
            push!(ts, t); push!(xs, x); push!(vs, v); push!(ps, pr)
            return (ts = ts, xs = xs, vs = vs, ps = ps, t_b = t, v_b = v, ok = true)
        end
    end
    return (ts = ts, xs = xs, vs = vs, ps = ps, t_b = NaN, v_b = v, ok = false)
end

# Cale fω pour que la vitesse de bouche soit celle visée (bissection).
function calibrate_impetus(v_target, L)
    lo, hi = 1e-2, 1e4
    for _ in 1:80
        mid = 0.5 * (lo + hi)
        r = integrate_bore(mid, v_target, L)
        (!r.ok || r.v_b < v_target) ? (lo = mid) : (hi = mid)
    end
    return 0.5 * (lo + hi)
end

# Interpolation linéaire sur une trajectoire tabulée.
function _interp(ts, ys, t)
    (t <= ts[1]) && return ys[1]
    (t >= ts[end]) && return 0.0
    i = searchsortedfirst(ts, t)
    i <= 1 && return ys[1]
    i > length(ys) && return 0.0
    w = (t - ts[i-1]) / (ts[i] - ts[i-1])
    return (1 - w) * ys[i-1] + w * ys[i]
end

# Interface inchangée pour l'aval : x(t), v(t), t_b, τ_v — mais tout est
# désormais dérivé, plus rien n'est posé. τ_v est obtenu par différences
# finies sur l'impétus (varier la charge, relire (v_bouche, t_b)).
# MÉMOÏSATION INDISPENSABLE. La cinématique ne dépend que de (v, L), mais son
# calcul coûte 85 intégrations de ~150 000 pas (calage de fω par bissection, plus
# cinq tirs pour τ_v). simulate_shot est appelé des dizaines de fois par balayage
# et des dizaines de milliers de fois par le Monte-Carlo : sans cache, le coût
# devient prohibitif alors que le résultat est rigoureusement identique.
const _KIN_CACHE = Dict{Tuple{Float64,Float64},Any}()

projectile_kinematics(v_muzzle, L) =
    get!(_KIN_CACHE, (v_muzzle, L)) do
        _projectile_kinematics(v_muzzle, L)
    end

function _projectile_kinematics(v_muzzle, L)
    fω = calibrate_impetus(v_muzzle, L)
    tr = integrate_bore(fω, v_muzzle, L)
    pts = Tuple{Float64,Float64}[]
    for k in (0.92, 0.96, 1.0, 1.04, 1.08)
        r = integrate_bore(fω * k, v_muzzle, L)
        r.ok && push!(pts, (r.v_b, r.t_b))
    end
    n  = length(pts)
    mv = sum(q[1] for q in pts) / n
    mt = sum(q[2] for q in pts) / n
    τv = abs(sum((q[1]-mv)*(q[2]-mt) for q in pts) / sum((q[1]-mv)^2 for q in pts))
    return (
        x   = t -> _interp(tr.ts, tr.xs, t),
        v   = t -> _interp(tr.ts, tr.vs, t),
        p   = t -> _interp(tr.ts, tr.ps, t),
        t_b = tr.t_b,
        τ_v = τv,
        fω  = fω,
        p_peak = maximum(tr.ps),
    )
end

# -----------------------------------------------------------------------------
# 7. EXCITATION
# -----------------------------------------------------------------------------

function consistent_point_load(F, x, ndof_a)
    if x <= 0 || x >= L
        return zeros(ndof_a)
    end
    e = clamp(Int(floor(x / L_e)) + 1, 1, N_elements)
    ξ̄ = (x - (e - 1) * L_e) / L_e
    Ns = (1 - 3ξ̄^2 + 2ξ̄^3,
          L_e * (ξ̄ - 2ξ̄^2 + ξ̄^3),
          3ξ̄^2 - 2ξ̄^3,
          L_e * (-ξ̄^2 + ξ̄^3))
    Fa = zeros(ndof_a)
    for (k, idx_g) in enumerate((2*e - 1, 2*e, 2*e + 1, 2*e + 2))
        idx_a = idx_g - 2
        if 1 <= idx_a <= ndof_a
            Fa[idx_a] += F * Ns[k]
        end
    end
    return Fa
end

# Vecteur force global à l'instant t (h_offset paramétrable)
#
# moment_of_t : permet d'INJECTER un autre historique de moment à la culasse, en
# N·m, au lieu du p(t)·A_bore·h par défaut. Sert à comparer des architectures
# dont le profil de pression n'a pas la même forme — typiquement un PCP, dont la
# détente est plate et longue là où la combustion est brève et pointue. Sans ce
# point d'entrée on ne peut comparer que des amplitudes, jamais des FORMES, ce
# qui est précisément la question quand on demande si un PCP excite peu.
# m_proj : masse du projectile pour la charge mobile (b), indépendante de m_p.
function force_vector(t, x_p_of_t, ndof_a, h_offset;
                      moment_of_t = nothing, m_proj = m_p, p_of_t = nothing)
    F = zeros(ndof_a)
    # (a) Moment de recul à la culasse → d.d.l. de rotation du nœud 2.
    #     La pression vient de la balistique intérieure COUPLÉE (section 6) :
    #     c'est la même p(t) qui accélère la balle et qui pousse la culasse, ce
    #     qui rend ∫p·A dt = m_eff·v exact et interdit la dérive corrigée le
    #     2026-07-19 (l'ancien profil autonome valait 5,19× le recul physique).
    F[2] += moment_of_t !== nothing ? moment_of_t(t) :
            (p_of_t === nothing ? 0.0 : p_of_t(t)) * A_bore * h_offset
    # (b) Poids du projectile (charge mobile)
    xp = x_p_of_t(t)
    if 0 < xp < L
        F .+= consistent_point_load(-m_proj * g_accel, xp, ndof_a)
    end
    return F
end

# -----------------------------------------------------------------------------
# 8. AMORTISSEMENT DE RAYLEIGH
# -----------------------------------------------------------------------------
function rayleigh_damping(Ma, Ka, ω1, ω2, ζ1, ζ2)
    Amat = [1/(2ω1)  ω1/2;
            1/(2ω2)  ω2/2]
    α_R, β_R = Amat \ [ζ1; ζ2]
    return α_R * Ma + β_R * Ka, α_R, β_R
end

# -----------------------------------------------------------------------------
# 9. NEWMARK-β (γ=1/2, β=1/4)
# -----------------------------------------------------------------------------
function newmark_solve(Ma, Ca, Ka, F_of_t, t_end, Δt; γ = 0.5, β = 0.25)
    n  = size(Ma, 1)
    ts = collect(0:Δt:t_end)
    Nt = length(ts)
    U  = zeros(n, Nt)
    V  = zeros(n, Nt)
    Ac = zeros(n, Nt)
    Ac[:, 1] = Ma \ (F_of_t(ts[1]) - Ca * V[:, 1] - Ka * U[:, 1])
    K_eff = Ma + γ * Δt * Ca + β * Δt^2 * Ka
    K_fact = factorize(K_eff)
    for i in 1:Nt-1
        u_pred = U[:, i] + Δt * V[:, i] + Δt^2 * (0.5 - β) * Ac[:, i]
        v_pred = V[:, i] + Δt * (1 - γ) * Ac[:, i]
        rhs    = F_of_t(ts[i+1]) - Ca * v_pred - Ka * u_pred
        Ac[:, i+1] = K_fact \ rhs
        U[:, i+1]  = u_pred + β * Δt^2 * Ac[:, i+1]
        V[:, i+1]  = v_pred + γ * Δt   * Ac[:, i+1]
    end
    return ts, U, V, Ac
end

# -----------------------------------------------------------------------------
# 10. SIMULATION D'UN TIR (h_offset paramétrable)
# -----------------------------------------------------------------------------
function simulate_shot(m_tuner; J_tuner = tuner_inertia(m_tuner),
                       d_overhang = 0.0,
                       Δt = 5e-6, t_end = 30e-3,
                       ζ1 = ZETA_DEFAULT, ζ2 = ZETA_DEFAULT,
                       h_offset = h_offset_default,
                       moment_of_t = nothing, v_p = v_muzzle, m_proj = m_p,
                       ΔT = 0.0, K_root = Inf, verbose = true)
    Ka, Ma  = build_system(m_tuner, J_tuner; d_overhang = d_overhang, ΔT = ΔT,
                           K_root = K_root)
    freqs, ωs, Φ = modal_analysis(Ka, Ma; n_modes = 5)
    Ca, _, _ = rayleigh_damping(Ma, Ka, ωs[1], ωs[2], ζ1, ζ2)
    kin     = projectile_kinematics(v_p, L)
    t_b     = kin.t_b                    # découle du profil, plus imposé
    ndof_a  = size(Ka, 1)
    F_of_t  = t -> force_vector(t, kin.x, ndof_a, h_offset;
                                moment_of_t = moment_of_t, m_proj = m_proj,
                                p_of_t = kin.p)

    ts, U, V, _ = newmark_solve(Ma, Ca, Ka, F_of_t, t_end, Δt)

    y_L     = U[end-1, :]
    θ_L     = U[end,   :]
    θdot_L  = V[end,   :]

    idx_b      = argmin(abs.(ts .- t_b))
    θ_at_tb    = θ_L[idx_b]
    θdot_at_tb = θdot_L[idx_b]
    θdot_MOAms = moa_per_ms(θdot_at_tb)

    if verbose
        @printf("=== Tir simulé : tuner = %.0f g, h_offset = %.2f mm ===\n",
                m_tuner * 1e3, h_offset * 1e3)
        print("Fréquences propres (Hz) : ")
        for f in freqs; @printf("%8.2f  ", f); end
        println()
        @printf("τ_v (sensibilité) = %.2f µs/(m/s)   [Kolbe : 8,8]\n", kin.τ_v * 1e6)
        @printf("À t = t_b = %.3f ms :\n", ts[idx_b] * 1e3)
        @printf("    y(L)  = %+.3e m\n",      y_L[idx_b])
        @printf("    θ(L)  = %+.3e rad\n",    θ_at_tb)
        @printf("    θ̇(L)  = %+.3f rad/s  =  %+.3f MOA/ms\n",
                θdot_at_tb, θdot_MOAms)
    end

    return (
        m_tuner    = m_tuner,
        d_overhang = d_overhang,
        h_offset   = h_offset,
        freqs      = freqs,
        ωs         = ωs,
        Φ          = Φ,
        ts         = ts,
        y_L        = y_L,
        θ_L        = θ_L,
        θdot_L     = θdot_L,
        idx_b      = idx_b,
        t_b        = ts[idx_b],
        θ_at_tb    = θ_at_tb,
        θdot_at_tb = θdot_at_tb,
        θdot_MOAms = θdot_MOAms,
    )
end

# -----------------------------------------------------------------------------
# 11. CALIBRATION AUTOMATIQUE de h_offset
#
# Le système est linéaire en h_offset (l'excitation est proportionnelle, et la
# structure est linéaire), donc UN SEUL tir suffit : on multiplie h_offset par
# le ratio de la cible sur la valeur mesurée. La cible est le pic absolu de
# |θ̇(L,t)| sur l'historique temporel, qui correspond à l'enveloppe vibratoire.
# -----------------------------------------------------------------------------
function calibrate_h_offset(target_peak_MOAms; m_tuner = m_tuner_0, kwargs...)
    h_ref = 1e-3                    # 1 mm de référence
    res = simulate_shot(m_tuner;
                        h_offset = h_ref, verbose = false, kwargs...)
    peak_rad = maximum(abs.(res.θdot_L))
    peak_MOAms = moa_per_ms(peak_rad)
    h_cal = h_ref * (target_peak_MOAms / peak_MOAms)
    @printf("Calibration : pic |θ̇| visé = %.2f MOA/ms  →  h_offset = %.2f mm\n",
            target_peak_MOAms, h_cal * 1e3)
    return h_cal
end

# -----------------------------------------------------------------------------
# 11 bis. BRAS DE LEVIER PHYSIQUE (ancrage sur les paramètres de l'arme)
#
# Le modèle encastré ne contient pas le mouvement de corps rigide du recul ;
# on relie donc le bras de levier h_offset au moment de recul physique
# M0 = p·A_bore·h_cg en calant UNE constante (h_cal_ref, obtenue par
# calibrate_h_offset sur la mesure de Kolbe pour l'arme de référence), puis en
# appliquant la dépendance documentée : ∝ h_cg (bras de levier du recul) et
# ∝ 1/m_rifle (une arme lourde recule, et donc vibre, moins). On retombe sur
# h_cal_ref pour l'arme de référence (m_rifle_ref, h_cg_ref).
# -----------------------------------------------------------------------------
function physical_h_offset(h_cal_ref; m_rifle = m_rifle_ref, h_cg = h_cg_ref)
    @assert m_rifle > 0 "La masse de l'arme doit être strictement positive."
    return h_cal_ref * (h_cg / h_cg_ref) * (m_rifle_ref / m_rifle)
end

# -----------------------------------------------------------------------------
# 12. BALAYAGE PARAMÉTRIQUE
# -----------------------------------------------------------------------------
function tuner_sweep(m_range; θdot_target = θdot_optimum_MOAms, kwargs...)
    println("\n========== Balayage paramétrique sur m_tuner ==========")
    @printf("%10s | %9s | %14s | %14s | %10s\n",
            "m (g)", "f₁ (Hz)", "θ(L,t_b)(µrad)", "θ̇(L,t_b) MOA/ms", "|écart|")
    println("-"^72)
    results = Any[]
    for m in m_range
        res = simulate_shot(m; verbose = false, kwargs...)
        ecart = abs(res.θdot_MOAms - θdot_target)
        @printf("%10.1f | %9.2f | %+14.2f | %+14.3f | %10.3f\n",
                m * 1e3, res.freqs[1], res.θ_at_tb * 1e6,
                res.θdot_MOAms, ecart)
        push!(results, res)
    end
    println("-"^72)
    ecarts = [abs(r.θdot_MOAms - θdot_target) for r in results]
    idx = argmin(ecarts)
    @printf("Optimum : m_tuner ≈ %.0f g  →  f₁ = %.2f Hz, θ̇(t_b) = %+.3f MOA/ms (cible %.1f)\n",
            results[idx].m_tuner * 1e3, results[idx].freqs[1],
            results[idx].θdot_MOAms, θdot_target)
    return results
end

# -----------------------------------------------------------------------------
# 12 bis. BALAYAGE EN POSITION — LE VRAI RÉGLAGE DE TERRAIN
#
# Sur le terrain la masse du tuner est FIXE (on monte un poids) et l'accord se
# fait en VISSANT le tuner plus ou moins loin en porte-à-faux devant la bouche.
# Ce balayage reproduit donc la procédure réelle : à masse fixée, on scanne la
# course et on cherche le porte-à-faux qui amène θ̇(t_b) sur la cible de
# compensation positive (+6 MOA/ms, Kolbe), et non le zéro visé en PCP.
# -----------------------------------------------------------------------------
function position_sweep(d_range; m_tuner = m_tuner_0,
                        θdot_target = θdot_optimum_MOAms, kwargs...)
    @printf("\n========== Balayage en position (masse fixe %.0f g) ==========\n",
            m_tuner * 1e3)
    @printf("%12s | %9s | %14s | %14s | %10s\n",
            "porte-à-faux", "f₁ (Hz)", "θ(L,t_b)(µrad)", "θ̇(L,t_b) MOA/ms", "|écart|")
    println("-"^74)
    results = Any[]
    for d in d_range
        res = simulate_shot(m_tuner; d_overhang = d, verbose = false, kwargs...)
        @printf("%9.1f mm | %9.2f | %+14.2f | %+14.3f | %10.3f\n",
                d * 1e3, res.freqs[1], res.θ_at_tb * 1e6,
                res.θdot_MOAms, abs(res.θdot_MOAms - θdot_target))
        push!(results, res)
    end
    println("-"^74)
    ecarts = [abs(r.θdot_MOAms - θdot_target) for r in results]
    idx = argmin(ecarts)
    @printf("Optimum : porte-à-faux ≈ %.1f mm  →  f₁ = %.2f Hz, θ̇(t_b) = %+.3f MOA/ms (cible %.1f)\n",
            results[idx].d_overhang * 1e3, results[idx].freqs[1],
            results[idx].θdot_MOAms, θdot_target)
    # Tolérance : plage de porte-à-faux restant à moins de 1 MOA/ms de la cible.
    ok = [r.d_overhang * 1e3 for r in results if abs(r.θdot_MOAms - θdot_target) <= 1.0]
    isempty(ok) || @printf("Optimum large : %.0f–%.0f mm à moins de 1 MOA/ms de la cible (%.0f mm de tolérance)\n",
                           minimum(ok), maximum(ok), maximum(ok) - minimum(ok))
    return results
end

# -----------------------------------------------------------------------------
# 13. TRACÉS (Plots.jl)
# -----------------------------------------------------------------------------
function plot_shot(res; save_path = "")
    t_ms = res.ts .* 1e3

    p_y = plot(t_ms, res.y_L .* 1e6,
               xlabel = "Temps (ms)", ylabel = "y(L) (µm)",
               title  = "Déflexion verticale à la bouche",
               lw = 1.5, color = :steelblue, legend = false)
    vline!(p_y, [res.t_b * 1e3], color = :red, ls = :dash, label = "t_b")

    p_θ = plot(t_ms, res.θ_L .* 1e6,
               xlabel = "Temps (ms)", ylabel = "θ(L) (µrad)",
               title  = "Angle de bouche",
               lw = 1.5, color = :darkgreen, legend = false)
    vline!(p_θ, [res.t_b * 1e3], color = :red, ls = :dash)
    hline!(p_θ, [0.0], color = :black, ls = :dot, alpha = 0.5)

    p_dθ = plot(t_ms, moa_per_ms.(res.θdot_L),
                xlabel = "Temps (ms)", ylabel = "θ̇(L) (MOA/ms)",
                title  = "Vitesse angulaire de bouche",
                lw = 1.5, color = :darkorange, legend = false)
    vline!(p_dθ, [res.t_b * 1e3], color = :red, ls = :dash)
    hline!(p_dθ, [θdot_optimum_MOAms], color = :green, ls = :dot,
           alpha = 0.7)
    hline!(p_dθ, [0.0], color = :black, ls = :dot, alpha = 0.5)

    fig = plot(p_y, p_θ, p_dθ, layout = (3, 1), size = (800, 900),
               plot_title = @sprintf("Tir nominal — tuner %.0f g, h_offset %.2f mm",
                                     res.m_tuner * 1e3, res.h_offset * 1e3))
    if !isempty(save_path)
        savefig(fig, save_path)
        println("→ Tracé du tir enregistré : $save_path")
    end
    return fig
end

# `sweeps` : un ou plusieurs balayages en position (un par masse de tuner).
# Superposer deux masses encadre la fourchette réelle du .22 LR — cf. les poids
# fabricants relevés dans la page wiki (~115 à 230 g).
function plot_position_sweep(sweeps...; save_path = "", θdot_target = θdot_optimum_MOAms)
    couleurs = [:darkorange, :purple, :teal]

    # NB : le point de θ̇ ne s'affiche pas dans la police par défaut de GR ;
    # on écrit « taux angulaire » en toutes lettres dans titres et libellés.
    p1 = plot(xlabel = "Porte-à-faux du tuner (mm)",
              ylabel = "Taux angulaire à t_b (MOA/ms)",
              title  = "Le réglage réel : taux angulaire de bouche vs position",
              legend = :outertop)
    # Bande de tolérance : ±1 MOA/ms autour de la cible (tracée en premier
    # pour rester sous les courbes).
    hspan!(p1, [θdot_target - 1, θdot_target + 1], color = :green, alpha = 0.10,
           label = "Cible Kolbe ±1 MOA/ms")
    hline!(p1, [θdot_target], color = :green, ls = :dot, label = "")

    p2 = plot(xlabel = "Porte-à-faux du tuner (mm)", ylabel = "f₁ (Hz)",
              title = "Fréquence fondamentale vs position", legend = false)

    for (k, results) in enumerate(sweeps)
        ds  = [r.d_overhang * 1e3 for r in results]
        m_g = results[1].m_tuner * 1e3
        c   = couleurs[mod1(k, length(couleurs))]
        plot!(p1, ds, [r.θdot_MOAms for r in results],
              label = @sprintf("Tuner %.0f g", m_g),
              lw = 2, marker = :circle, ms = 3, color = c)
        plot!(p2, ds, [r.freqs[1] for r in results],
              lw = 2, marker = :circle, ms = 3, color = c)
    end

    masses = join([@sprintf("%.0f", s[1].m_tuner * 1e3) for s in sweeps], " et ")
    fig = plot(p1, p2, layout = (2, 1), size = (800, 660),
               plot_title = "Accord en position, à masse fixe — $masses g (.22 LR)")
    if save_path != ""
        savefig(fig, save_path)
        println("Figure sauvegardée : $save_path")
    end
    return fig
end

function plot_sweep(results; save_path = "", θdot_target = θdot_optimum_MOAms)
    ms        = [r.m_tuner * 1e3        for r in results]
    f1s       = [r.freqs[1]             for r in results]
    θdot_arr  = [r.θdot_MOAms           for r in results]
    θ_arr     = [r.θ_at_tb * 1e6        for r in results]

    p1 = plot(ms, θdot_arr,
              xlabel = "m_tuner (g)", ylabel = "θ̇(L,t_b) (MOA/ms)",
              title  = "Vitesse angulaire à t_b vs masse du tuner",
              lw = 2, marker = :circle, color = :darkorange, legend = :outertop)
    hline!(p1, [θdot_target], color = :green, ls = :dot, label = "Cible Kolbe ($θdot_target)")
    hline!(p1, [0.0], color = :black, ls = :dot, alpha = 0.5, label = "")

    p2 = plot(ms, θ_arr,
              xlabel = "m_tuner (g)", ylabel = "θ(L,t_b) (µrad)",
              title  = "Angle de bouche à t_b vs masse du tuner",
              lw = 2, marker = :circle, color = :darkgreen, legend = false)
    hline!(p2, [0.0], color = :black, ls = :dot, alpha = 0.5)

    p3 = plot(ms, f1s,
              xlabel = "m_tuner (g)", ylabel = "f₁ (Hz)",
              title  = "Fréquence fondamentale vs masse du tuner",
              lw = 2, marker = :circle, color = :steelblue, legend = false)

    fig = plot(p1, p2, p3, layout = (3, 1), size = (800, 900),
               plot_title = "Balayage paramétrique du tuner")
    if !isempty(save_path)
        savefig(fig, save_path)
        println("→ Tracé du balayage enregistré : $save_path")
    end
    return fig
end

function plot_modes(res; n_modes = 4, save_path = "")
    # Reconstruction des modes : on récupère les composantes y(x) sur chaque nœud
    # y₁ = 0 (encastrement), puis y_n correspond au d.d.l. actif (2n-3)
    xs_nodes = collect(0:L_e:L)
    Φ = res.Φ
    fig = plot(layout = (n_modes, 1), size = (800, 700),
               plot_title = "Modes propres du canon + tuner")
    for k in 1:min(n_modes, size(Φ, 2))
        # extraire y aux nœuds (composantes translation)
        ys = [0.0]
        for n in 2:N_nodes
            push!(ys, Φ[2*(n-1) - 1, k])  # d.d.l. actif 2(n-1)-1 = y_n
        end
        plot!(fig[k], xs_nodes .* 1e3, ys,
              xlabel = "Position le long du canon (mm)",
              ylabel = "φ($k)(x)",
              title  = @sprintf("Mode %d — f = %.2f Hz", k, res.freqs[k]),
              lw = 2, marker = :circle, legend = false, color = :purple)
        hline!(fig[k], [0.0], color = :black, ls = :dot, alpha = 0.5)
    end
    if !isempty(save_path)
        savefig(fig, save_path)
        println("→ Tracé des modes enregistré : $save_path")
    end
    return fig
end

# -----------------------------------------------------------------------------
# 14. EXÉCUTION — sous garde PROGRAM_FILE, comme les autres scripts du
#     dépôt : simulation.jl s'inclut ainsi comme BIBLIOTHÈQUE (cf.
#     variability.jl) sans rejouer tous ses balayages.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 13 bis. PARAMÈTRES DE MODÈLE ISSUS DES SOURCES EXTERNES
# -----------------------------------------------------------------------------
# Étape A — Bras de levier h_offset, ancré sur une MESURE (Vaughn 1998)
#
# HISTORIQUE. Jusqu'au 2026-07-18, h_offset était fixé en calant le PIC absolu
# de |θ̇(L,t)| sur 10 MOA/ms — un chiffre choisi par nous pour son « amplitude
# plausible », que personne n'a mesuré. Toutes les amplitudes du modèle
# reposaient donc sur un paramètre libre. On l'ancre désormais sur la seule
# mesure publiée de l'EXCITATION elle-même.
#
# SOURCE. H. R. Vaughn, « Rifle Accuracy Facts » (Precision Shooting, 1998),
# ch. 4 : jauges de contrainte sur l'anneau de culasse d'une .270 Win à
# 53 000 psi. Piège que Vaughn signale lui-même — les jauges lisent un moment
# de RÉPONSE (~450 in-lb), bien plus petit que le moment APPLIQUÉ (~1500 in-lb,
# fig. 4-33) parce que le canon « can't respond quickly enough ». C'est le
# moment APPLIQUÉ qui correspond à notre excitation :
#     1500 in-lb / (53 000 psi × A_âme(.277)) → bras de levier ≈ 11,9 mm.
#
# TRANSPORT .270 → .22 LR. La différence de poussée est déjà portée par le
# p·A_bore du modèle ; h ne transporte que la géométrie et la dynamique de
# l'arme, via la loi documentée en physical_h_offset (h ∝ h_cg, h ∝ 1/m_rifle) :
#     h(.22 LR) = 11,9 mm × (h_cg_22 / h_cg_270) × (m_270 / m_22)
# Les deux h_cg sont du même ordre (~25 mm, axe d'âme → CG). Vaughn ne publie
# PAS la masse totale de son arme (seulement le canon, 2,8 lb) ; en retenant
# 3,9 kg pour une .270 sporter avec lunette, h(.22 LR) ≈ 9,3 mm — à lire comme
# un ORDRE DE GRANDEUR (~10 mm, ni 5 ni 40), l'estimation héritant d'une masse
# supposée, d'un ansatz non validé et d'un report centerfire → rimfire.
#
# POURQUOI ON NE RETIENT PAS 9,3 mm. Essayé : le modèle plafonne alors à
# 3,5 MOA/ms et n'atteint plus JAMAIS l'optimum de compensation de 6,0, en
# masse comme en position — plus aucun sweet spot. C'est cohérent avec l'étape
# [C ter] : il manque à l'ossature encastrée la rotation d'ensemble de l'arme.
# Un h « physiquement juste » dans un modèle amputé donne un modèle
# physiquement faux ET inutilisable.
#
# CE QU'ON RETIENT. h_offset n'est PLUS lisible comme un bras de levier. Avec
# l'ancien profil de pression (impulsion 5,19× le recul physique) il valait
# 16,5 mm, soit ~1,6× le bras physique borné par Vaughn (~9-12 mm), et l'on
# écrivait que ce rapport « chiffrait ce que l'encastrement omet ». C'ÉTAIT FAUX :
# le facteur était contaminé par l'erreur de profil, qui poussait dans le même
# sens. Le modèle sur-excitait, tout en sous-prédisant Kolbe d'un facteur ~100 ;
# h_offset absorbait silencieusement les deux défauts à la fois, et le calage
# « marchait » pour cette seule raison.
#
# Avec la balistique intérieure couplée (section 6), où la conservation de la
# quantité de mouvement est exacte par construction, la valeur requise pour
# atteindre la compensation (+6,0 MOA/ms) monte à 62,5 mm — soit ~6× le bras
# physique. À ce niveau la grandeur excède toute cote de l'arme : ce n'est plus
# un bras de levier mais un FACTEUR DE GAIN, dont la seule lecture honnête est
# qu'il mesure l'ampleur du mécanisme absent (rotation d'ensemble supprimée par
# l'encastrement). Le chiffre est laid, et c'est sa vertu : il ne se déguise plus
# en grandeur géométrique plausible.
const H_OFFSET_EFF = 62.5e-3

# Mesures de Kolbe (2015) sur .22 LR, taux d'angle de bouche à la sortie.
# Servent au diagnostic sans dimension de l'étape [C ter].
const KOLBE_BARE, KOLBE_TUNED = -9.4, 6.0


if abspath(PROGRAM_FILE) == @__FILE__

    println("="^72)
    println(" Simulation MEF + Newmark-β des vibrations transversales d'un canon")
    println(" Canon 26\" / .22 LR / Tuner réglable / Application à Kolbe (2015)")
    println("="^72)
    println()


    println("[A] Bras de levier h_offset (effectif ; référent physique = Vaughn 1998)")
    println("-"^72)
    h_cal = H_OFFSET_EFF
    @printf("h_offset = %.2f mm (EFFECTIF) — physique d'après Vaughn ≈ 9-12 mm,\n", h_cal * 1e3)
    @printf("           soit ×%.1f-%.1f : la part que l'encastrement de culasse omet.\n",
            h_cal * 1e3 / 12, h_cal * 1e3 / 9.3)
    # Pour mémoire, l'ancienne voie — un pic arbitraire de 10 MOA/ms — reste
    # disponible mais n'a AUCUN référent mesuré :
    #     h_cal = calibrate_h_offset(10.0; m_tuner = m_tuner_0)
    println()

    # Étape B — Tir nominal avec le bras de levier effectif
    println("[B] Tir nominal (tuner $(Int(m_tuner_0 * 1000)) g, h_offset effectif)")
    println("-"^72)
    result_nom = simulate_shot(m_tuner_0; h_offset = h_cal)
    println()

    # Étape B bis — Sensibilité à l'arme (poids + hauteur âme↔CG)
    # h_cal est l'étalon obtenu pour l'arme de référence (m_rifle_ref, h_cg_ref) ;
    # on prédit ici l'effet de monter le même canon sur des armes différentes.
    println("[B bis] Sensibilité à l'arme : amplitude ∝ h_cg, ∝ 1/m_rifle")
    println("-"^72)
    @printf("%12s | %14s | %14s | %16s\n",
            "m_rifle (kg)", "h_cg (mm)", "h_offset (mm)", "pic |θ̇| (MOA/ms)")
    println("-"^72)
    for (m_rifle, h_cg) in [(m_rifle_ref, h_cg_ref), (3.5, h_cg_ref),
                            (7.0, h_cg_ref), (m_rifle_ref, 0.0381),
                            (m_rifle_ref, 0.0127)]
        h_eff = physical_h_offset(h_cal; m_rifle = m_rifle, h_cg = h_cg)
        res   = simulate_shot(m_tuner_0; h_offset = h_eff, verbose = false)
        peak  = moa_per_ms(maximum(abs.(res.θdot_L)))
        @printf("%12.1f | %14.1f | %14.3f | %16.2f\n",
                m_rifle, h_cg * 1e3, h_eff * 1e3, peak)
    end
    println("-"^72)
    println("→ Arme plus légère ou âme plus haute ⇒ vibrations plus amples (retuning).")
    println()

    # Étape C — Balayage paramétrique
    result_sweep = tuner_sweep(0.0:0.025:0.4; h_offset = h_cal)
    println()

    # Étape C bis — Balayage en position, à masse fixe : LE réglage de terrain.
    # Deux masses encadrant la fourchette réelle des tuners .22 LR, d'après les
    # poids fabricants : PMA 4,2 oz (119 g) et EC v2 ~4 oz en bas, Harrell's rimfire
    # 8 oz (227 g) en haut. Le tube carbone Starik/Centra — l'un des plus répandus —
    # pèse ~200-220 g ENSEMBLE COMPLET, sa bague mobile seule étant bien plus
    # légère : il relève donc du bas de la fourchette.
    const M_TUNER_LOW  = 0.100
    const M_TUNER_HIGH = 0.200
    # Plage étendue à 200 mm le 2026-07-19. Avec k dérivé de la masse, l'optimum
    # à 100 g tombe à 100 mm, soit exactement le dernier point de l'ancienne plage
    # (0-100 mm) : un optimum au bord n'en est pas un. Un premier élargissement à
    # 150 mm rendait l'optimum intérieur mais tronquait encore la TOLÉRANCE (la
    # bande à moins de 1 MOA/ms se refermait au-delà) — on aurait publié une
    # largeur de plage fixée par la fenêtre de calcul, non par la physique. À
    # 200 mm les deux bornes sont atteintes pour les deux masses.
    result_pos_low  = position_sweep(0.0:0.005:0.20; m_tuner = M_TUNER_LOW,  h_offset = h_cal)
    result_pos_high = position_sweep(0.0:0.005:0.20; m_tuner = M_TUNER_HIGH, h_offset = h_cal)
    println()

    # Étape C ter — ÉCART À COMBLER PAR L'ACCORD, rapporté à ce que l'accord fournit.
    #
    # POURQUOI CETTE FORME, ET PAS UN RAPPORT nu/accordé. La version antérieure
    # de ce diagnostic divisait par θ̇(canon nu). C'était une faute de
    # construction : θ̇(t_b) est la valeur d'une OSCILLATION à un instant donné,
    # elle peut tomber n'importe où, y compris tout près de zéro — et c'est
    # précisément ce qui se produit une fois l'inertie du tuner corrigée
    # (θ̇(nu) = −0,09). Un dénominateur qui traverse zéro ne fait pas un
    # indicateur : le rapport explosait (×3,6 hier, ×100 aujourd'hui) sans que
    # rien de physique n'ait changé entre les deux.
    #
    # On garde l'idée — une grandeur SANS DIMENSION, donc insensible au choix
    # de h_offset, qui multiplie numérateur et dénominateur à l'identique —
    # mais on la rapporte à une référence STABLE : le maximum de θ̇ que
    # l'accord sait produire. D'où
    #
    #     R = θ̇(canon non accordé) / max θ̇(accordé)
    #
    # qui se lit « de combien le canon non accordé est-il en dessous de la
    # cible, rapporté à ce que l'accord peut fournir ». C'est exactement le
    # travail demandé au tuner.
    #
    # CE QUE LE RÉSULTAT LOCALISE. Le modèle donne R ≈ −1,5 % : son canon non
    # accordé est déjà quasiment sur la cible, le tuner n'a presque rien à
    # rattraper. Kolbe mesure R = −157 % : son canon part très en dessous et
    # l'accord doit remonter plus que toute sa course. L'écart n'est donc PAS
    # dans le mécanisme d'accord — la plage 0 → +6,2 MOA/ms que le modèle
    # produit est du bon ordre — mais dans l'ÉTAT NON ACCORDÉ, c'est-à-dire
    # dans l'excitation de base. C'est la signature attendue de la rotation
    # d'ensemble absente : elle tire fortement la bouche vers le bas sur un
    # canon nu, sans changer grand-chose à ce qu'un tuner peut ensuite ajouter.
    #
    # RÉSERVE MAINTENUE. Le −9,4 de Kolbe dépend de la souplesse non publiée de
    # son étau (cf. ROADMAP.md) — soit justement le mécanisme manquant. Ce
    # diagnostic ne prouve donc pas que le modèle a tort : il chiffre, sur une
    # grandeur que la calibration ne peut pas flatter, l'ampleur de ce qu'il
    # laisse de côté.
    println("[C ter] Écart à combler par l'accord (sans dimension, insensible à h_offset)")
    println("-"^72)
    θdot_bare = result_sweep[1].θdot_MOAms           # m_tuner = 0, canon non accordé
    R_kolbe   = KOLBE_BARE / KOLBE_TUNED
    for (lbl, θdot_max) in (("accord en masse",    maximum(r.θdot_MOAms for r in result_sweep)),
                            ("accord en position", maximum(r.θdot_MOAms for r in result_pos_low)))
        R = θdot_bare / θdot_max
        @printf("  %-19s θ̇(nu) = %+6.2f, max atteignable %+5.2f  →  R = %+7.1f %%  (Kolbe %+.0f %%)\n",
                lbl, θdot_bare, θdot_max, 100R, 100R_kolbe)
    end
    @printf("  ⇒ le modèle demande au tuner ~%.0f fois moins de travail que la mesure.\n",
            abs(R_kolbe / (θdot_bare / maximum(r.θdot_MOAms for r in result_pos_low))))
    println("    L'écart est dans l'état NON ACCORDÉ, non dans le mécanisme d'accord :")
    println("    signature de la rotation d'ensemble que l'encastrement supprime.")
    println()

    # Étape D bis — Dérive thermique (ajoutée le 2026-07-19 d'après Dai et al.)
    println("[D bis] Dérive thermique de l'accord")
    println("-"^72)
    println("Le canon s'échauffe en série ; E baisse, la fréquence propre avec lui.")
    # Le réglage de référence doit être l'OPTIMUM COURANT, non une cote figée :
    # coder 110 mm en dur laissait un écart de 0,45 mm dès ΔT = 0 après la
    # révision de l'amortissement, qui a déplacé l'optimum à 135 mm.
    d_ref = result_pos_low[argmin([abs(r.θdot_MOAms - θdot_optimum_MOAms)
                                   for r in result_pos_low])].d_overhang
    @printf("Effet à RÉGLAGE FIXÉ (tuner %d g à %.0f mm, l'optimum courant), cible 6,0 :\n",
            round(Int, M_TUNER_LOW*1e3), d_ref*1e3)
    println()
    @printf("  %6s | %8s | %10s | %s\n", "ΔT", "f₁ (Hz)", "θ̇(t_b)", "traînée résiduelle à 50 m")
    println("  " * "-"^62)
    let trainee = g_accel * D_target^2 * 10.0 / v_muzzle^3 * 1e3
        for dT in (0.0, 30.0, 60.0, 100.0)
            r = simulate_shot(M_TUNER_LOW; d_overhang = d_ref, h_offset = h_cal,
                              ΔT = dT, verbose = false)
            res = trainee * abs(r.θdot_MOAms - θdot_optimum_MOAms) / θdot_optimum_MOAms
            @printf("  %+5.0f°C | %8.2f | %+10.3f | %.2f mm\n",
                    dT, r.freqs[1], r.θdot_MOAms, res)
        end
    end
    println()
    println("  EXISTE-T-IL DES RÉGLAGES PLUS TOLÉRANTS ? Question explorée, réponse")
    println("  NÉGATIVE. Un balayage minimax — minimiser l'écart MAXIMAL sur toute")
    println("  la plage 0-60 °C plutôt que l'écart à froid — redonne exactement")
    println("  l'optimum classique, aux deux masses testées. Aucun compromis ne se")
    println("  justifie.")
    println()
    println("  La raison est une question d'échelles : la sensibilité vaut")
    println("  0,027 MOA/ms par mm de position contre 0,0012 par degré, si bien que")
    println("  60 °C d'échauffement équivalent à 2,7 mm de tuner (un degré ≈ 44 µm).")
    println("  L'erreur de position domine, et viser juste reste la bonne stratégie.")
    println()
    println("  Il existe bien un point thermiquement NEUTRE (dérive nulle), à ~5 mm")
    println("  de l'optimum, mais s'y placer coûte trois fois plus en écart à la")
    println("  cible que d'accepter la dérive à l'optimum.")
    println()
    println("  Ordre de grandeur : quelques dixièmes de millimètre à 50 m sur une")
    println("  série de match (+30 à +60 °C). En tir soutenu (plusieurs centaines de")
    println("  degrés), la dérive dépasserait la course même du tuner.")
    println()

    # Étape D — Tracés
    println("[D] Génération des tracés")
    println("-"^72)
    plot_shot(result_nom; save_path = "plot_tir_nominal.png")
    plot_position_sweep(result_pos_low, result_pos_high;
                        save_path = "plot_balayage_position.png")
    plot_sweep(result_sweep; save_path = "plot_balayage_tuner.png")
    plot_modes(result_nom; save_path = "plot_modes_propres.png", n_modes = 4)
    println()

    println("="^72)
    println("Remarques :")
    @printf(" • h_offset = FACTEUR DE GAIN (%.0f mm), non un bras de levier : ~6× le bras\n", H_OFFSET_EFF*1e3)
    println("   physique (Vaughn ≈ 9-12 mm). Il mesure l'ampleur du mécanisme absent.")
    println(" • Balistique intérieure COUPLÉE : ∫p·A dt = m_eff·v exact par construction.")
    @printf("   τ_v est désormais PRÉDIT (%.2f µs/(m/s)) et non calé — Kolbe mesure 8,8.\n",
            projectile_kinematics(v_muzzle, L).τ_v * 1e6)
    println(" • Les tracés sont enregistrés en PNG dans le répertoire courant.")
    println(" • Amortissement de Rayleigh : ζ₁ = 0.5 %, ζ₂ = 1 %.")
    println(" • Schéma de Newmark : (γ, β) = (1/2, 1/4), Δt = 5 µs.")
    println("="^72)


end