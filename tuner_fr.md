---
title: |
  Dynamique et Optimisation des Vibrations de Canon\
  De la théorie balistique à la simulation par éléments finis
---

# Introduction : principes de la compensation positive {#sec:intro}

Cette première section présente le *pourquoi* et le *comment* de la
compensation positive sans recours aux outils mathématiques. Un lecteur
curieux qui souhaite seulement comprendre l'idée et son intérêt pratique
peut s'arrêter à la fin de la
section [1.8](#sec:plan){reference-type="ref" reference="sec:plan"} ;
les parties suivantes développent ensuite la modélisation rigoureuse et
la simulation numérique.

## Le défi du tir de précision

Dans les disciplines de tir de précision (couché à 50 mètres ISSF par
exemple), la recherche du groupement parfait se heurte aux limites
physiques de la munition. Même avec des lots de cartouches de qualité
*match*, il subsiste une variation inévitable de la vitesse de sortie du
canon --- la *vitesse initiale* --- d'une cartouche à l'autre. Cette
variation, souvent mesurée en pieds par seconde (fps), est typiquement
de l'ordre de quelques pour cent.

## Pourquoi une variation de vitesse fait varier le point d'impact

Une balle plus lente passe plus de temps en vol pour atteindre la cible
et subit donc l'accélération de la pesanteur plus longtemps : elle tombe
davantage. Mécaniquement, en cible, cela se traduit par un impact plus
bas qu'une balle rapide tirée avec le même réglage. Si le canon était un
tube parfaitement rigide et immobile, la dispersion des vitesses des
cartouches se traduirait inévitablement par une *traînée verticale* sur
la cible : les balles rapides en haut, les lentes en bas.

À titre d'ordre de grandeur, la chute vaut
$h = \tfrac12\,g\,(D/\bar v)^2$ pour une distance $D$ et une vitesse
initiale $\bar v$, si bien qu'un écart de vitesse $\Delta v$ décale
l'impact de $|\Delta h| \approx g\,D^2\,\Delta v / \bar v^3$. Pour une
.22 LR de match ($\bar v \approx 330$ m/s) tirée à $D = 50$ m, une
dispersion de $30$ fps ($\Delta v \approx 9$ m/s) produit ainsi, par le
seul temps de vol, une traînée verticale d'environ $6$ mm --- largement
de quoi ruiner un groupement de match. C'est précisément cette
dispersion que la compensation positive vise à annuler.

## Le canon vibre comme un diapason

En réalité, lorsque le coup part, le canon ne reste pas immobile : il se
met à vibrer. L'explosion de la poudre, le recul de l'arme, la course du
projectile dans le tube et son frottement contre les rayures excitent
des oscillations transversales --- c'est-à-dire de bas en haut et de
gauche à droite --- qui font bouger la bouche du canon de quelques
centièmes de degré, à la manière d'un diapason qui sonne après un choc.
Ces vibrations durent plusieurs dizaines de millisecondes, bien au-delà
du temps que met la balle pour parcourir le canon (typiquement 2 à 3
millisecondes).

Pendant ces 2 à 3 millisecondes, la bouche du canon est donc en
mouvement : elle monte, redescend, repasse par sa position d'équilibre,
etc. L'angle exact avec lequel la balle sort dépend non seulement de la
visée du tireur, mais aussi de *l'instant précis* de la sortie.

## Ce qui détermine l'ampleur des vibrations

D'où vient cette mise en branle ? La source dominante, soulignée par
G. Kolbe `\cite{KolbeSim}`{=latex}, est le *recul* : sous la poussée des
gaz, l'arme entière part en arrière et *pivote autour de son centre de
gravité*. Or l'axe du canon (la ligne de l'âme, par où sort la balle)
passe *au-dessus* de ce centre de gravité ; la poussée de recul s'exerce
donc avec un *bras de levier* et imprime un couple à l'arrière du canon,
exactement comme une chiquenaude donnée hors du centre fait basculer une
règle posée sur la table. Ce couple est le « coup de marteau » qui fait
sonner le diapason.

Deux grandeurs concrètes commandent alors l'*amplitude* du mouvement de
bouche :

-   le **poids total de l'arme** : une arme lourde encaisse le recul en
    bougeant moins, donc vibre moins ;

-   la **hauteur de l'âme au-dessus du centre de gravité** : plus ce
    bras de levier est grand, plus le couple --- et donc l'amplitude ---
    est important.

Ces deux paramètres n'apparaissent pas dans le réglage du tuner (qui
agit sur le *rythme* des vibrations, pas sur leur cause), mais ils
expliquent pourquoi deux carabines de canon identique mais de montage
différent ne se règlent pas pareil. Ils seront réintroduits formellement
à la section [3.4](#sec:excitation){reference-type="ref"
reference="sec:excitation"} comme le moment d'excitation $M_0(t)$.

## Le principe de la compensation positive

L'idée géniale, proposée dès 1901 par
A. Mallock `\cite{Mallock}`{=latex} et formalisée pour le tir moderne
par G. Kolbe (Border Barrels) `\cite{Kolbe,KolbeWeb}`{=latex}, est
d'*exploiter* ces vibrations plutôt que de les subir. Si l'on parvient à
régler le canon de sorte que la bouche soit en train de *monter* au
moment précis où la balle en sort, alors :

-   une balle **rapide** sort un peu *en avance* : à ce moment, le canon
    n'a pas encore eu le temps de beaucoup remonter, et la balle est
    lancée avec un angle plus *bas*
    (figure [1](#fig:compensation_positive){reference-type="ref"
    reference="fig:compensation_positive"}, trajectoire rouge) ;

-   une balle **lente** sort un peu *en retard* : le canon a continué à
    se relever, et la balle est lancée avec un angle plus *haut*
    (trajectoire bleue).

La balle lente est donc projetée *plus haut* initialement ; en cible,
son surcroît de chute (parce qu'elle a volé plus longtemps) est
exactement compensé par son surcroît de hauteur de départ. Toutes les
balles, rapides ou lentes, peuvent ainsi se retrouver au *même* point
d'impact. C'est ce que Kolbe appelle la *positive compensation*.

<figure id="fig:compensation_positive">

<figcaption>Principe géométrique de la compensation positive. Si la
bouche du canon est en train de monter lorsque les balles en sortent,
les balles plus lentes (qui sortent plus tard) bénéficient d’un angle de
lancement supérieur. Leur surcroît de hauteur de départ compense le
surcroît de chute dû à leur vol prolongé, et toutes les balles peuvent
se retrouver au même point d’impact en cible.</figcaption>
</figure>

## Le rôle du tuner

Pour que ce phénomène se produise, encore faut-il *accorder* les
vibrations du canon avec le temps de sortie de la balle --- d'où le
terme *tuner*, qui désigne une masse mobile (typiquement 100 à 400
grammes) fixée à la bouche du canon. En vissant ou en dévissant cette
masse, le tireur la rapproche ou l'éloigne de la bouche, ce qui modifie
subtilement la fréquence de vibration du canon. Cela revient à ajuster
le rythme du « diapason » : on rallonge ou raccourcit le temps qu'il
faut à la bouche pour repasser par sa position d'équilibre, afin que cet
instant tombe pile au moment où la balle moyenne sort.

Tout changement de munition ou même de lot demande en principe un
nouveau réglage, puisque la vitesse moyenne --- et donc le temps de
parcours dans le canon --- change. Les tireurs procèdent typiquement par
essais successifs (méthode dite du *ladder* ou
*Audette* `\cite{Audette}`{=latex}) : on tire plusieurs séries en
changeant la position du tuner d'un cran à chaque fois, et on retient la
position pour laquelle les groupements se resserrent le plus.

## Synthèse intuitive

En résumé, la compensation positive demande trois ingrédients :

1.  un canon qui vibre suffisamment pour bouger sensiblement entre la
    sortie d'une balle rapide et celle d'une balle lente ;

2.  un *accord* de ces vibrations qui place le moment de sortie de la
    balle moyenne au bon point de l'oscillation (la phase *ascendante*)
    ;

3.  un compromis quantitatif entre vitesse du mouvement angulaire et
    chute balistique : un canon qui se relève trop lentement ne
    compensera pas assez la dispersion, et un canon qui se relève trop
    vite la sur-compensera.

Les sections suivantes formalisent chacun de ces trois ingrédients et
montrent comment les paramètres physiques du canon, du tuner et de la
munition se combinent pour atteindre l'optimum.

## Plan du document {#sec:plan}

Le reste du document développe la théorie et la simulation :

-   Section [2](#sec:notations){reference-type="ref"
    reference="sec:notations"} : liste des notations et conventions
    utilisées.

-   Section [3](#sec:math){reference-type="ref" reference="sec:math"} :
    modélisation mathématique du canon (poutre d'Euler-Bernoulli,
    méthode des éléments finis, intégration par Newmark-$\beta$).

-   Section [4](#sec:modale){reference-type="ref"
    reference="sec:modale"} : analyse des fréquences propres et
    sensibilité au tuner.

-   Section [5](#sec:cinematique){reference-type="ref"
    reference="sec:cinematique"} : critère mathématique de la
    compensation positive, distinction nœud / ventre temporel,
    confirmation expérimentale de Kolbe et validation par la simulation
    numérique.

-   Section [6](#sec:conclusions){reference-type="ref"
    reference="sec:conclusions"} : conclusions pratiques pour le tireur.

# Notations et conventions {#sec:notations}

On adopte le système d'unités SI
($\text{m}, \text{kg}, \text{s}, \text{Pa}, \text{rad}$). L'axe $x$ est
aligné avec l'âme du canon, orienté de la culasse ($x=0$) vers la bouche
($x=L$). L'axe $y$ est transverse, orienté vers le haut. Sauf mention
contraire, les angles sont en radians ; les conversions usuelles sont
$1\ \text{MOA} = \pi/(180 \times 60)\ \text{rad} \approx 2{,}909 \times 10^{-4}\ \text{rad}$
et $1\ \text{mrad} \approx 3{,}438\ \text{MOA}.$

::: center
  Symbole                                           Signification                                                                  Unité
  ------------------------------------------------- ------------------------------------------------------------------------------ ----------------
  *Géométrie et matériau du canon*                                                                                                 
  $L$                                               Longueur du canon                                                              m
  $D_\text{ext}, D_\text{int}$                      Diamètres extérieur (profil) et intérieur (âme)                                m
  $A$                                               Aire de section transverse $\tfrac{\pi}{4}(D_\text{ext}^2 - D_\text{int}^2)$   m$^2$
  $I$                                               Moment quadratique $\tfrac{\pi}{64}(D_\text{ext}^4 - D_\text{int}^4)$          m$^4$
  $E$                                               Module d'Young (acier $\approx 200$ GPa)                                       Pa
  $\rho$                                            Masse volumique (acier $\approx 7850$ kg/m$^3$)                                kg/m$^3$
  $EI$                                              Rigidité en flexion                                                            N$\cdot$m$^2$
  $\rho A$                                          Masse linéique                                                                 kg/m
  *Cinématique de la poutre*                                                                                                       
  $y(x,t)$                                          Déflexion transverse au point $x$ et à l'instant $t$                           m
  $\theta(x,t) \equiv \partial y/\partial x$        Pente / angle de fibre neutre                                                  rad
  $\theta(L,t)$                                     Angle à la bouche du canon (angle de tir)                                      rad
  $\dot\theta = \partial \theta/\partial t$         Vitesse angulaire de bouche                                                    rad/s
  *Tuner*                                                                                                                          
  $m_t$                                             Masse du tuner (typiquement 100--400 g)                                        kg
  $J_t$                                             Moment d'inertie transverse du tuner                                           kg$\cdot$m$^2$
  *Modes propres*                                                                                                                  
  $\omega_n$                                        Pulsation propre du mode $n$                                                   rad/s
  $f_n = \omega_n/(2\pi)$                           Fréquence propre du mode $n$                                                   Hz
  $T_n = 1/f_n$                                     Période propre du mode $n$                                                     s
  $\phi_n(x)$                                       Forme spatiale du mode $n$                                                     ---
  $\zeta_n$                                         Taux d'amortissement modal                                                     ---
  *Projectile et balistique interne*                                                                                               
  $m_p$                                             Masse du projectile (.22 LR : $\approx 2{,}6$ g)                               kg
  $v(x,t)$                                          Vitesse du projectile à l'abscisse $x$                                         m/s
  $v_0, v_\text{muzzle}$                            Vitesse initiale (à la sortie de la bouche)                                    m/s
  $\bar v$                                          Vitesse moyenne effective dans le canon                                        m/s
  $t_b$                                             Temps de parcours du projectile dans le canon                                  s
  $p(t)$                                            Pression dans la chambre / l'âme                                               Pa
  $A_\text{bore} = \tfrac{\pi}{4} D_\text{int}^2$   Aire de la section utile de l'âme                                              m$^2$
  $\tau_v \equiv -\partial t_b/\partial v_0$        Sensibilité du temps de sortie à la vitesse initiale                           s/(m/s)
  *Cible et compensation*                                                                                                          
  $D$                                               Distance de tir (50 m pour la discipline ISSF couché)                          m
  $g$                                               Accélération de la pesanteur ($\approx 9{,}81$ m/s$^2$)                        m/s$^2$
  $h_\text{bal}(v_0)$                               Chute balistique à la distance $D$                                             m
  $\theta_\text{out}$                               Angle de bouche au moment $t_b$ : $\theta(L, t_b)$                             rad
  $\dot\theta_\text{out}^\star$                     Vitesse angulaire *cible* pour compensation complète                           rad/s
  *Simulation numérique*                                                                                                           
  $N$                                               Nombre d'éléments finis                                                        ---
  $L_e = L/N$                                       Longueur d'un élément                                                          m
  $[K], [M], [C]$                                   Matrices globales de raideur, masse, amortissement                             ---
  $\{u(t)\}$                                        Vecteur des d.d.l. nodaux (translations + rotations)                           ---
  $\Delta t$                                        Pas de temps d'intégration Newmark-$\beta$                                     s
  $\gamma, \beta$                                   Paramètres du schéma de Newmark                                                ---
  $h_\text{offset}$                                 Bras de levier du moment de recul à la culasse (excitation)                    m
:::

# Développement Mathématique {#sec:math}

## Modèle continu : poutre d'Euler--Bernoulli

On modélise le canon comme une poutre prismatique encastrée à la culasse
($x=0$) et libre à la bouche ($x=L$), de section $A(x)$ et de moment
quadratique $I(x)$ (admis constants par morceaux ; les profils *bull*,
*sporter* et coniques se traitent par variation par élément). En notant
$y(x,t)$ la déflexion transversale, l'équation aux dérivées partielles
régissant les petites oscillations vaut : $$\label{eq:EB}
        \rho A \,\frac{\partial^2 y}{\partial t^2} \;+\; \frac{\partial^2}{\partial x^2}\!\left(EI \,\frac{\partial^2 y}{\partial x^2}\right) \;=\; q(x,t),$$
où $q(x,t)$ représente le chargement transverse linéique. Les conditions
aux limites sont :
$$y(0,t) = 0,\qquad \frac{\partial y}{\partial x}(0,t) = 0,\qquad EI\,\frac{\partial^2 y}{\partial x^2}\bigg|_{x=L}\!\! = M_t,\qquad \frac{\partial}{\partial x}\!\left(EI\,\frac{\partial^2 y}{\partial x^2}\right)\bigg|_{x=L}\!\! = -F_t,$$
où $M_t$ et $F_t$ représentent respectivement le moment et la force
appliqués par le tuner (incluant ses effets inertiels). Pour un canon de
profil match à parois épaisses, l'hypothèse d'Euler--Bernoulli (sections
planes restant planes et normales à la fibre neutre) reste raisonnable
tant que $L/D_\text{ext} \gtrsim 20$ ; au-delà, un modèle de Timoshenko
intégrant l'inertie de rotation et la déformation par cisaillement est
préférable.

<figure id="fig:schema_canon">

<figcaption>Modèle du canon : poutre d’Euler-Bernoulli encastrée à la
culasse (<span class="math inline"><em>x</em> = 0</span>), libre à la
bouche (<span class="math inline"><em>x</em> = <em>L</em></span>),
portant à son extrémité une masse ponctuelle et son inertie de rotation
<span
class="math inline">(<em>m</em><sub><em>t</em></sub>,<em>J</em><sub><em>t</em></sub>)</span>
représentant le tuner. La courbe pointillée illustre une déflexion
transverse <span
class="math inline"><em>y</em>(<em>x</em>,<em>t</em>)</span>.</figcaption>
</figure>

## Discrétisation par éléments finis

On cherche une solution faible de
[\[eq:EB\]](#eq:EB){reference-type="eqref" reference="eq:EB"}. En
multipliant par une fonction test $w(x)$ de classe $\mathcal{C}^1$
vérifiant les conditions essentielles et en intégrant deux fois par
parties sur $[0,L]$, on obtient la *formulation variationnelle* :
$$\label{eq:weak}
        \int_0^L \rho A \,\ddot y\, w \,\mathrm{d}x \;+\; \int_0^L EI\, y'' w'' \,\mathrm{d}x \;=\; \int_0^L q\,w\,\mathrm{d}x \;+\; [\text{termes de bord}].$$

Le canon est partitionné en $N$ éléments de longueur $L_e$. Sur chaque
élément, on attribue deux degrés de liberté par nœud : la déflexion $y$
et la rotation $\theta = y'$. On exprime le champ local $y_e(\xi)$, où
$\xi \in [0,L_e]$ est la coordonnée locale, à l'aide des **fonctions de
forme d'Hermite cubiques** : $$\label{eq:hermite}
        y_e(\xi) = \mathbf{N}(\xi)\, \mathbf{u}_e, \qquad \mathbf{u}_e = \begin{pmatrix} y_1 \\ \theta_1 \\ y_2 \\ \theta_2 \end{pmatrix},$$
avec $$\begin{aligned}
        N_1(\xi) &= 1 - 3\bar\xi^2 + 2\bar\xi^3, & N_2(\xi) &= L_e\,(\bar\xi - 2\bar\xi^2 + \bar\xi^3), \\
        N_3(\xi) &= 3\bar\xi^2 - 2\bar\xi^3, & N_4(\xi) &= L_e\,(-\bar\xi^2 + \bar\xi^3),
    
\end{aligned}$$ où $\bar\xi = \xi/L_e$. Ces polynômes garantissent la
continuité $\mathcal{C}^1$ aux nœuds, exigée par l'opérateur
biharmonique de [\[eq:weak\]](#eq:weak){reference-type="eqref"
reference="eq:weak"}.

En reportant [\[eq:hermite\]](#eq:hermite){reference-type="eqref"
reference="eq:hermite"} dans
[\[eq:weak\]](#eq:weak){reference-type="eqref" reference="eq:weak"}, on
obtient pour chaque élément la matrice de raideur élémentaire
$$^e \;=\; \int_0^{L_e} EI\, (\mathbf{N}'')^{\!\top} \mathbf{N}'' \,\mathrm{d}\xi
        \;=\; \frac{EI}{L_e^3}
        \begin{bmatrix}
            12     & 6L_e   & -12    & 6L_e   \\
            6L_e   & 4L_e^2 & -6L_e  & 2L_e^2 \\
            -12    & -6L_e  & 12     & -6L_e  \\
            6L_e   & 2L_e^2 & -6L_e  & 4L_e^2
        \end{bmatrix},$$ et la matrice de masse cohérente
$$^e \;=\; \int_0^{L_e} \rho A\, \mathbf{N}^{\!\top} \mathbf{N} \,\mathrm{d}\xi
        \;=\; \frac{\rho A L_e}{420}
        \begin{bmatrix}
            156    & 22L_e  & 54     & -13L_e \\
            22L_e  & 4L_e^2 & 13L_e  & -3L_e^2 \\
            54     & 13L_e  & 156    & -22L_e \\
            -13L_e & -3L_e^2 & -22L_e & 4L_e^2
        \end{bmatrix}.$$ Le vecteur de chargement nodal s'obtient
analoguement par
$\mathbf{f}_e = \int_0^{L_e} q(\xi,t)\, \mathbf{N}^{\!\top}\,\mathrm{d}\xi$.

## Assemblage, conditions aux limites et intégration du tuner

Les matrices globales $[K]$ et $[M]$ de taille $2(N+1) \times 2(N+1)$
sont assemblées par sommation des contributions élémentaires sur les
degrés de liberté partagés. L'encastrement à la culasse impose
$y(0,t)=0$ et $\theta(0,t)=0$ : on supprime les deux premières lignes et
colonnes du système, ce qui définit les sous-matrices actives $[K]_a$ et
$[M]_a$.

L'ajout du tuner à la bouche (dernier nœud d'indice $N+1$) se traduit
par l'addition de sa masse $m_t$ sur le degré de liberté de translation
et de son moment d'inertie $J_t$ (calculé autour de l'axe
perpendiculaire à l'axe du canon, dans le plan de vibration) sur le
degré de rotation : $$_a \;\longleftarrow\; [M]_a \;+\;
        \underbrace{\begin{bmatrix}
                & \mathbf{0} & \\[2pt]
                \mathbf{0} & m_t & 0 \\[2pt]
                & 0 & J_t
        \end{bmatrix}}_{\text{aux d.d.l. de la bouche}}.$$ Si le tuner
est modélisé comme un cylindre de masse $m_t$, de rayon $R_t$ et de
longueur $\ell_t$, son moment d'inertie autour de l'axe transverse
passant par son centre de gravité vaut
$J_t = m_t\,(3R_t^2 + \ell_t^2)/12$ ; un décalage $d$ entre le centre de
gravité du tuner et la bouche introduit un terme additionnel $m_t d^2$
(théorème d'Huygens) ainsi qu'un couplage masse--rotation $m_t d$ qu'il
convient de prendre en compte pour un tuner long et déporté.

## Modélisation de l'excitation dynamique {#sec:excitation}

L'excitation $\{F(t)\}$ regroupe plusieurs contributions hétérogènes en
amplitude et en bande spectrale :

#### (a) Pression des gaz et recul.

La pression $p(t)$ dans la chambre engendre une force longitudinale
appliquée à la culasse. Si l'arme est tenue de manière non parfaitement
symétrique (ce qui est toujours le cas), cette force induit un moment
transverse $M_0(t)$ au nœud d'encastrement. On modélise typiquement
$p(t)$ par un profil de Pierret ou un ajustement gaussien sur les
courbes pression-temps mesurées, d'amplitude pic $p_{\max} \sim 200$ MPa
pour la .22 LR et de durée caractéristique $\sim 0{,}5$ ms.

Physiquement, ce moment n'est pas un artefact d'asymétrie de tenue : il
existe même pour une arme tenue idéalement, car la *ligne de l'âme est
décalée du centre de gravité* de l'arme. Suivant la modélisation de
Kolbe `\cite{KolbeSim}`{=latex}, l'arme libre de masse $m_r$ recule et
tourne autour de son centre de gravité, situé à une distance
$h_\text{cg}$ *sous* l'axe de l'âme ; la force de recul
$F_\text{rec}(t) = p(t)\,A_\text{bore}$ s'exerce donc avec ce bras de
levier et imprime à la base du canon le moment $$\label{eq:moment_recul}
        M_0(t) \;=\; p(t)\,A_\text{bore}\,h_\text{cg},$$ tandis que
l'amplitude angulaire résultante décroît comme $1/m_r$ (une arme lourde
recule moins). Le modèle encastré ne contient pas le mouvement de corps
rigide du recul ; on relie donc le bras de levier $h_\text{offset}$ de
la section [5.5](#sec:validation_num){reference-type="ref"
reference="sec:validation_num"} au moment de recul physique en
factorisant sa dépendance, $$\label{eq:hoffset_phys}
        h_\text{offset}(m_r, h_\text{cg}) \;=\; h_\text{offset}^\text{réf}\,\frac{h_\text{cg}}{h_\text{cg}^\text{réf}}\,\frac{m_r^\text{réf}}{m_r},$$
où la *seule* constante $h_\text{offset}^\text{réf}$ est calée une fois
pour toutes sur l'enveloppe vibratoire mesurée par Kolbe pour son arme
de référence ($m_r^\text{réf} = 5$ kg,
$h_\text{cg}^\text{réf} = 25{,}4$ mm). Au lieu de recalibrer pour chaque
arme, on *prédit* alors la dépendance à deux grandeurs mesurables sur
l'arme réelle --- le poids total $m_r$ (aisé à peser) et la distance
$h_\text{cg}$ de l'âme au centre de gravité (plus délicate, déterminée
par équilibrage) : monter le même canon sur une arme plus légère, ou
dont l'âme est plus haut placée au-dessus du centre de gravité, amplifie
les vibrations et impose un nouveau réglage du tuner.

#### (a$'$) Conditions aux limites : encastrement *vs* action souple.

Le modèle ci-dessus suppose un encastrement parfait à la culasse
($y(0)=0,\ \theta(0)=0$). C'est une idéalisation :
Kolbe `\cite{KolbeSim}`{=latex} préfère ne *pas* bloquer la base et
représente la souplesse de la boîte de culasse par un tronçon
supplémentaire de poutre ($\sim 100$ mm de long, $\sim 38$ mm de
diamètre, alésage $\sim 25$ mm) ajouté à l'arrière, l'ensemble reculant
*librement* dans l'espace (hypothèse *free recoil*, bien approchée par
un tir au sac mais mise en défaut si une arme légère est fermement
épaulée). Notre encastrement rigide surestime donc légèrement la
fréquence fondamentale et ignore le mouvement de corps rigide du recul ;
il reste néanmoins pertinent pour l'*angle relatif* de bouche
$\theta(L,t)$, seule grandeur qui pilote la compensation
(cf. section [5](#sec:cinematique){reference-type="ref"
reference="sec:cinematique"}), et simplifie l'analyse modale. Le
raffinement « action souple » constitue l'extension naturelle pour un
travail quantitatif sur arme réelle.

#### (b) Charge mobile du projectile (effet de poids).

Le projectile, de masse $m_p$, applique en sa position courante $x_p(t)$
une force transverse due à la gravité (poids) et à l'engagement dans les
rayures. En première approximation, on traite le projectile comme une
force ponctuelle mobile :
$$q_p(x,t) \;=\; -m_p\,g\,\delta(x - x_p(t)),$$ où $g$ est
l'accélération de la pesanteur et $\delta$ la distribution de Dirac. Le
vecteur nodal correspondant s'obtient en évaluant les fonctions de forme
à la position $x_p(t)$.

#### (c) Couple gyroscopique des rayures.

L'accélération angulaire imprimée par les rayures transmet un couple
réactif au canon, modélisable comme un moment réparti proportionnel à
$\dot v(x_p)$. Cet effet, plus faible que (a) et (b), peut être négligé
en première analyse.

#### (d) Position du projectile.

La cinématique interne $x_p(t)$ est obtenue par intégration de
l'équation de l'écoulement interne :
$$m_p\,\ddot x_p \;=\; p(t)\, A_b \;-\; F_\text{frot}(x_p, \dot x_p),$$
où $A_b$ est la section utile du tube et $F_\text{frot}$ la résistance
d'engagement. Le temps de sortie $t_b$ vérifie $x_p(t_b) = L$ ; pour la
.22 LR Match dans un canon de 26\", $t_b \approx 1{,}5$ à $1{,}7$ ms.

## Schéma d'intégration temporelle de Newmark--$\beta$

L'équation semi-discrète de la structure s'écrit $$\label{eq:dynsemi}
        [M]\{\ddot u(t)\} \;+\; [C]\{\dot u(t)\} \;+\; [K]\{u(t)\} \;=\; \{F(t)\},$$
où $[C]$ est usuellement modélisé par un amortissement de Rayleigh
$[C] = \alpha_M [M] + \alpha_K [K]$ avec deux constantes
$(\alpha_M,\alpha_K)$ calibrées sur les taux d'amortissement modaux
mesurés ($\zeta_1 \sim 0{,}5\,\%$ à $1\,\%$ pour un canon en acier).

On utilise le schéma implicite de Newmark--$\beta$ pour avancer du pas
$n$ au pas $n+1$ ($\Delta t$) : $$\begin{aligned}
        \{u_{n+1}\} &= \{u_n\} + \Delta t\, \{\dot u_n\} + \tfrac{\Delta t^2}{2}\!\big((1-2\beta)\{\ddot u_n\} + 2\beta\,\{\ddot u_{n+1}\}\big), \label{eq:newmark1} \\
        \{\dot u_{n+1}\} &= \{\dot u_n\} + \Delta t\,\big((1-\gamma)\{\ddot u_n\} + \gamma\,\{\ddot u_{n+1}\}\big), \label{eq:newmark2}
    
\end{aligned}$$ dans lesquelles $\{\ddot u_{n+1}\}$ est obtenu en
réinjectant [\[eq:newmark1\]](#eq:newmark1){reference-type="eqref"
reference="eq:newmark1"} dans
[\[eq:dynsemi\]](#eq:dynsemi){reference-type="eqref"
reference="eq:dynsemi"} évaluée en $t_{n+1}$, ce qui conduit au système
linéaire :
$$\Big([M] + \gamma\,\Delta t\,[C] + \beta\,\Delta t^2\,[K]\Big)\,\{\ddot u_{n+1}\} \;=\; \{F_{n+1}\} - [C]\{\dot u_n^\star\} - [K]\{u_n^\star\},$$
où $\{u_n^\star\}$ et $\{\dot u_n^\star\}$ sont les prédicteurs
construits à partir des valeurs courantes. Le couple classique
$(\gamma,\beta) = (1/2, 1/4)$ (*average constant acceleration*) est
*inconditionnellement stable* et conserve l'énergie sans dissipation
numérique ; il est privilégié ici pour préserver fidèlement l'amplitude
des modes excités. Le pas de temps doit néanmoins résoudre la fréquence
d'intérêt la plus élevée : on prendra $\Delta t \lesssim T_{\max}/20$,
où $T_{\max}$ est la période du plus haut mode physiquement significatif
(typiquement quelques kHz pour un canon, donc $\Delta t \sim 10\,\mu$s).

# Analyse modale et sensibilité au tuner {#sec:modale}

## Problème aux valeurs propres

En l'absence d'amortissement et d'excitation, on cherche des solutions
harmoniques $\{u(t)\} = \{\phi\}\,e^{\mathrm{i}\omega t}$, ce qui ramène
[\[eq:dynsemi\]](#eq:dynsemi){reference-type="eqref"
reference="eq:dynsemi"} à un problème aux valeurs propres généralisé :
$$\label{eq:eig}
        \big([K] - \omega^2 [M]\big)\,\{\phi\} = \mathbf{0}.$$ Les
solutions $(\omega_n^2,\{\phi_n\})$ donnent les pulsations propres
$\omega_n$ et les modes propres $\{\phi_n\}$. Les fréquences en hertz
sont $f_n = \omega_n/(2\pi)$. Pour le canon de référence du simulateur
(acier, $L=0{,}66$ m, $D_\text{ext}=24$ mm cylindrique,
$D_\text{int}=5{,}6$ mm), une résolution MEF à 20 éléments donne les
premiers modes suivants :

::: center
  Mode                                     Canon nu   Tuner 200 g   Tuner 400 g
  ------------ ------------------------------------ ------------- -------------
  $f_1$ (Hz)                                  39,87         34,16         30,33
  $f_2$ (Hz)                                  230,1        217,66         206,9
  $f_3$ (Hz)                                  633,5        602,39         575,1
  $f_4$ (Hz)     $\approx$`<!-- -->`{=html}1,19 kHz      1,14 kHz      1,09 kHz
  $f_5$ (Hz)     $\approx$`<!-- -->`{=html}1,87 kHz      1,79 kHz      1,72 kHz
:::

La solution analytique pour la poutre cantilever nue,
$f_1 = \tfrac{1{,}875^2}{2\pi}\sqrt{EI/(\rho A L^4)}$, donne
$f_1 \approx 40{,}0$ Hz, en parfait accord avec le calcul MEF
($39{,}87$ Hz). L'ajout d'un tuner de 200 g à la bouche fait chuter
$f_1$ d'environ $14\,\%$ (de 39,9 à 34,2 Hz) ; un tuner de 400 g donne
une réduction de $24\,\%$. Cette plage couvre largement la fenêtre
temporelle utile autour de $t_b \approx 2$--$3$ ms.

![Quatre premiers modes propres $\phi_n(x)$ du canon muni d'un tuner de
200 g (calcul MEF, 20 éléments). On notera la déformée monotone du mode
1 (ventre spatial à la bouche), et le rapprochement progressif des nœuds
spatiaux des modes supérieurs vers la bouche, conséquence de la
concentration de masse à
$x=L$.](plot_modes_propres.png){#fig:modes_propres
width="0.85\\linewidth"}

## Sensibilité par le quotient de Rayleigh {#sec:rayleigh}

Le *quotient de Rayleigh* associé au mode $n$ s'écrit
$$\omega_n^2 \;=\; \frac{\{\phi_n\}^{\!\top}[K]\{\phi_n\}}{\{\phi_n\}^{\!\top}[M]\{\phi_n\}}.$$
Soit $\phi_n^{(L)}$ la composante du mode au nœud de bouche (en
translation). Une perturbation $\delta m_t$ de la masse du tuner modifie
$[M]$ de $\delta m_t\,\phi_n^{(L)\,2}$ sur le dénominateur, sans toucher
à $[K]$. À l'ordre un, on obtient : $$\label{eq:dwdmt}
        \boxed{\;\frac{\mathrm{d}\omega_n^2}{\mathrm{d}m_t} \;\approx\; -\,\omega_n^2\,\frac{\big|\phi_n^{(L)}\big|^2}{\{\phi_n\}^{\!\top}[M]\{\phi_n\}}\;\le 0.\;}$$
Deux conclusions importantes en découlent :

1.  L'ajout de masse à la bouche *abaisse* toujours la fréquence (signe
    négatif).

2.  L'abaissement est d'autant plus marqué que le mode présente une
    grande amplitude à la bouche. Pour le mode fondamental d'une poutre
    encastrée-libre, $\phi_1^{(L)}$ est maximal, donc $f_1$ chute
    fortement ; pour le mode 2, le nœud vibratoire se rapproche de la
    bouche, et l'effet du tuner est plus modeste.

Ces propriétés fondent le mécanisme d'accord : en déplaçant la position
d'une masse mobile (ou en empilant des masselottes), le tireur balaye
continûment $f_1$ sur une plage typique de l'ordre de $\pm 15$ à
$\pm 30\,\%$.

## Décomposition modale et réponse en régime transitoire

En décomposant la réponse sur la base des modes propres normalisés en
masse, $\{u(t)\} = \sum_n q_n(t)\{\phi_n\}$, l'équation
[\[eq:dynsemi\]](#eq:dynsemi){reference-type="eqref"
reference="eq:dynsemi"} se découple en
$$\ddot q_n(t) + 2\zeta_n \omega_n\, \dot q_n(t) + \omega_n^2\, q_n(t) \;=\; f_n(t),\quad f_n(t) = \{\phi_n\}^{\!\top}\{F(t)\}.$$
Pour une excitation impulsionnelle (la signature pression+projectile
s'apparente à une impulsion de quelques centaines de microsecondes),
chaque mode démarre un régime oscillatoire libre. La réponse globale à
la bouche en rotation est alors
$$\theta(L,t) \;=\; \sum_n \phi_n^{\theta(L)}\, q_n(t),$$ où
$\phi_n^{\theta(L)}$ est la composante de rotation du mode $n$ au nœud
de bouche. En pratique, le premier mode domine la cinématique de bouche
aux temps caractéristiques $\sim t_b$ ; c'est donc essentiellement son
réglage qui est ciblé par le tuner.

#### Régime transitoire, et non onde stationnaire établie.

Une mise en garde s'impose, dont Kolbe `\cite{KolbeSim}`{=latex} fait le
cœur de son argumentaire : il est tentant de représenter « la façon dont
vibre un canon » par les *solutions analytiques en ondes stationnaires*
de l'équation des poutres. Or ces ondes ne se forment *pas* pendant la
fenêtre utile. Leur vitesse de phase est trop faible pour qu'un régime
stationnaire ait le temps de s'établir sur la durée
$t_b \approx 1$--$3$ ms du passage de la balle : à cet instant, le canon
est encore dans la *réponse transitoire* au moment impulsionnel de
recul, et non dans un régime périodique installé. C'est précisément
pourquoi le découpage modal ci-dessus est exploité ici en *réponse
forcée transitoire* (intégrée pas à pas par Newmark-$\beta$) plutôt
qu'en superposition de modes d'amplitudes figées, et pourquoi le
raisonnement de la section [5.3](#sec:nodes){reference-type="ref"
reference="sec:nodes"} privilégie le *nœud temporel* (l'état instantané
de l'oscillation à $t_b$) sur le nœud spatial (propriété d'un mode
pleinement établi).

# Cinématique de bouche et compensation positive {#sec:cinematique}

## Temps de parcours et angle de bouche

Soit $v(x,t)$ la vitesse du projectile à l'abscisse $x$ ; le temps de
sortie vaut $$t_b \;=\; \int_0^L \frac{\mathrm{d}x}{v(x)}.$$ L'angle
effectif de lancement, à la sortie, s'écrit
$$\theta_\text{out} \;=\; \theta(L, t_b) \;=\; \left.\frac{\partial y}{\partial x}\right|_{x=L,\,t=t_b}.$$
Une variation $\Delta \theta_\text{out}$ entraîne un déplacement
vertical sur cible (à la distance $D$, sans correction balistique au
premier ordre) : $$\Delta h \;\approx\; D\,\Delta \theta_\text{out}.$$
La géométrie sous-jacente est celle de la
figure [1](#fig:compensation_positive){reference-type="ref"
reference="fig:compensation_positive"} de la vulgarisation : ce que la
section [1](#sec:intro){reference-type="ref" reference="sec:intro"}
décrivait qualitativement est ici reformulé en termes d'angle de bouche
et d'effet sur cible.

## Couplage entre dispersion de vitesse et dispersion verticale {#sec:couplage}

Une variation $\Delta v_0$ de la vitesse initiale induit deux effets
cumulés :

1.  *Chute balistique additionnelle.* À la distance $D$, pour une
    trajectoire tendue, la chute sous la ligne de visée s'écrit
    $h_\text{bal}(v_0) \simeq g D^2/(2 v_0^2)$, d'où $$\label{eq:dhbal}
                \frac{\partial h_\text{bal}}{\partial v_0} \;\simeq\; -\,\frac{g D^2}{v_0^3} \;<\;0.$$
    Numériquement, à $D=50$ m et $v_0 \approx 308$ m/s (1010 ft/s, .22
    LR Match),
    $\partial h_\text{bal}/\partial v_0 \approx -0{,}84$ mm/(m/s), soit
    environ $-0{,}016$ MOA par ft/s ; ce nombre est aussi celui calculé
    par Kolbe `\cite{Kolbe}`{=latex} à partir d'un solveur de
    trajectoire balistique.

2.  *Décalage du temps de sortie.* La cinématique interne couple $t_b$
    et $v_0$ : à la limite, $t_b \propto 1/v_0$ donne
    $\partial t_b/\partial v_0 = -t_b/v_0$, mais la dépendance réelle
    est plus marquée car une charge plus énergique augmente la pression
    *partout* dans le tube et raccourcit $t_b$ davantage qu'un simple
    rapport cinématique. On notera donc
    $$\tau_v \;\equiv\; -\,\frac{\partial t_b}{\partial v_0} \;>\; 0,$$
    quantité directement mesurable à l'aide d'un chronographe couplé à
    un capteur de sortie de bouche. Kolbe
    rapporte `\cite{Kolbe}`{=latex}, pour Eley Tenex en canon de 26\",
    $$\tau_v \;\approx\; \frac{1\ \text{ms}}{375\ \text{ft/s}} \;=\; \frac{1\ \text{ms}}{114\ \text{m/s}} \;\approx\; 8{,}8\ \mu\text{s\,/\,(m/s)},$$
    et avance que cette valeur est, à $\pm 10\,\%$ près, une constante
    pour les canons rimfire de plus de 6 pouces.

La condition de compensation positive impose l'égalité des deux
contributions en cible : $$\label{eq:compensation}
        \boxed{\;\frac{\mathrm{d}\theta_\text{out}}{\mathrm{d}t}\bigg|_{t_b}\cdot \Delta t_b \;+\; \frac{1}{D}\frac{\partial h_\text{bal}}{\partial v_0}\Delta v_0 \;=\; 0.\;}$$
En substituant $\Delta t_b = -\tau_v\,\Delta v_0$ et
[\[eq:dhbal\]](#eq:dhbal){reference-type="eqref" reference="eq:dhbal"},
on obtient la *vitesse angulaire optimale* à la bouche au moment de la
sortie : $$\label{eq:dthetaopt}
        \dot \theta_\text{out}^\star \;=\; \frac{1}{D\,\tau_v}\,\frac{\partial h_\text{bal}}{\partial v_0} \;=\; -\,\frac{g\,D}{v_0^3\,\tau_v}.$$
*Le signe négatif de $\partial h_\text{bal}/\partial v_0$ et de
$\partial t_b/\partial v_0$ se composent pour rendre
$\dot\theta_\text{out}^\star$ positif : la bouche doit donc se relever
($\dot\theta>0$) au moment du départ du coup.* Numériquement, avec
$D=50$ m, $v_0 = 308$ m/s et $\tau_v = 8{,}8\,\mu$s/(m/s) :
$$\dot\theta_\text{out}^\star \;\approx\; 1{,}9\ \text{mrad/ms} \;\approx\; 6{,}6\ \text{MOA/ms},$$
à comparer aux $6{,}0$ MOA/ms mesurés expérimentalement par Kolbe
(l'écart de l'ordre de 10 % s'explique par la non-tendance de la
trajectoire à 50 m et par les approximations sur $\tau_v$).

## Nœud spatial *vs* nœud temporel : précision sémantique {#sec:nodes}

Le vocabulaire de *nœud* et *ventre*, hérité de l'acoustique des ondes
stationnaires, prête souvent à confusion lorsqu'il est appliqué au tir.
Deux significations physiquement distinctes coexistent :

#### Nœud/ventre spatial

--- points particuliers du *mode propre* $\phi_n(x)$ le long du canon :
un *nœud spatial* est un $x^\star$ tel que $\phi_n(x^\star)=0$
(déplacement transverse nul) ; un *ventre spatial* est un $x^\star$ où
$|\phi_n|$ atteint un maximum local. Pour le mode fondamental d'une
poutre encastrée-libre, le profil est monotone croissant de 0
(encastrement) à $\phi_1(L)$ (bouche). **La bouche est donc, par
construction, un ventre spatial du mode 1.** Aucun choix de tuner ne
déplace cette propriété : la bouche reste là où la balle sort, à $x=L$.

Les modes supérieurs ($n\ge 2$) présentent en revanche des nœuds
spatiaux *à l'intérieur* du canon. Placer le tuner précisément à un nœud
spatial du mode $n$ rend ce mode insensible à la masse ajoutée
(cf. extension de [\[eq:dwdmt\]](#eq:dwdmt){reference-type="eqref"
reference="eq:dwdmt"} à une masse positionnée hors de $x=L$). C'est
l'astuce des tuners dits *harmoniquement neutralisés* : on peut ajuster
$f_1$ sans remuer $f_2$ ou $f_3$, ce qui rend le réglage plus monotone
et reproductible.

#### Nœud/ventre temporel

--- points particuliers de l'oscillation *dans le temps* de l'angle de
bouche $\theta(L,t)$ : un *nœud temporel* est un instant $t^\star$ où
$\theta(L,t^\star)=0$ (passage par zéro) ; un *ventre temporel* est un
instant où $|\theta(L,t)|$ atteint un maximum
(i.e. $\dot\theta(L,t^\star)=0$). Pour un mode 1 quasi-sinusoïdal, ces
deux types d'instants alternent à un quart de période.

Le critère [\[eq:dthetaopt\]](#eq:dthetaopt){reference-type="eqref"
reference="eq:dthetaopt"} se traduit alors sans ambiguïté :

::: center
  Instant de sortie $t_b$               $\theta(L,t_b)$          $\dot\theta(L,t_b)$
  ---------------------------------- --------------------- -------------------------------
  **Ventre temporel** (*antinode*)    $\pm\Theta_1$ (max)                $0$
  **Nœud temporel** (*node*)                  $0$           $\pm\Theta_1\,\omega_1$ (max)
:::

**La compensation positive exige une sortie au voisinage d'un *nœud
temporel ascendant***, là où $\theta = 0$ et $\dot\theta > 0$ au
maximum. Sortir au sommet de l'oscillation (ventre temporel) est, au
contraire, le *pire* cas : non seulement $\dot\theta=0$ annule toute
compensation, mais on tire en plus avec un écart angulaire statique
maximal qui décale le point d'impact moyen.

<figure id="fig:noeud_temporel">

<figcaption>Distinction nœud / ventre temporel sur l’angle de bouche
<span class="math inline"><em>θ</em>(<em>L</em>,<em>t</em>)</span>. La
compensation positive impose <span
class="math inline"><em>t</em><sub><em>b</em></sub></span> au voisinage
d’un <em>nœud temporel ascendant</em> (point vert à <span
class="math inline"><em>t</em> = 0</span> et <span
class="math inline"><em>t</em> = <em>T</em><sub>1</sub></span>), où
<span class="math inline"><em>θ</em> = 0</span> mais <span
class="math inline"><em>θ̇</em></span> est maximal et positif. Sortir au
ventre temporel (points rouges, <span
class="math inline"><em>θ̇</em> = 0</span>) annule la
compensation.</figcaption>
</figure>

Ce qu'on appelle upward swing of the vibration at the muzzle dans
Kolbe `\cite{Kolbe,KolbeWeb}`{=latex} correspond précisément à ce *nœud
temporel ascendant*. Dans le jargon des compétiteurs (méthode *ladder*
ou *OCW* `\cite{Audette}`{=latex}), le node d'accord désigne plus
largement une zone de *robustesse* de groupement, c'est-à-dire un
plateau autour duquel un petit décalage de $t_b$ ne dégrade pas le
groupement. Au sens strict de la mécanique vibratoire, ce plateau est
précisément la fenêtre autour d'un nœud temporel ascendant --- où
$\ddot\theta(L,t_b)$ est proche de zéro et $\dot\theta(L,t_b)$ stationne
près de sa valeur extrémale, ce qui rend le tuning peu sensible aux
fluctuations de $t_b$.

## Confirmation expérimentale (Kolbe, 2015)

Kolbe `\cite{KolbeWeb}`{=latex} a vérifié directement le critère
[\[eq:dthetaopt\]](#eq:dthetaopt){reference-type="eqref"
reference="eq:dthetaopt"} sur un banc instrumenté mesurant l'angle de
bouche par voie optique (polariseur croisé) avec une porte
photo-détectrice repérant l'instant exact de sortie. Sur un canon Border
de 26\", calibre .22 LR, munition Eley EPS Tenex, deux configurations
sont comparées :

::: center
  Configuration               $\dot\theta(L,t_b)$ mesuré          Groupements en cible    
  --------------------- -------------------------------------- -------------------------- --
  Canon nu               $-9{,}4$ MOA/ms (bouche descendante)   cordon *vertical* marqué  
  Canon + tuner 200 g    $+6{,}0$ MOA/ms (bouche ascendante)      groupements *ronds*     
:::

La valeur $+6{,}0$ MOA/ms coïncide avec le critère
[\[eq:dthetaopt\]](#eq:dthetaopt){reference-type="eqref"
reference="eq:dthetaopt"} pour 50 m et la munition employée.
L'élimination quasi-complète du cordonnement vertical confirme
expérimentalement que le critère est non seulement nécessaire (signe)
mais également suffisant (magnitude) lorsque les modes d'ordre supérieur
restent faiblement excités.

En particulier, Kolbe note que la *vitesse de translation verticale*
$\dot y(L,t_b)$ de la bouche contribue de manière négligeable à la
dispersion en cible, comparée à $\dot\theta(L,t_b)$ : c'est bien la
*rotation* de la bouche, et non sa translation, qui détermine l'angle de
lancement effectif. Ce constat justifie *a posteriori* le découpage
d.d.l. (translation + rotation) du modèle MEF, et oriente le contrôle
vers l'extraction de la composante rotationnelle au nœud terminal.

## Validation par simulation MEF + Newmark-$\beta$ {#sec:validation_num}

Le simulateur en `Julia` associé à ce document (`simulation.jl`)
implémente l'ensemble du formalisme : assemblage MEF, encastrement,
tuner, balistique interne (modèle exponentiel calé sur $t_b$ et
$v_\text{muzzle}$), excitation par moment de recul à la culasse
$M(t) = p(t)\,A_\text{bore}\,h_\text{offset}$, charge mobile du
projectile, amortissement de Rayleigh et schéma de Newmark-$\beta$. Le
bras de levier $h_\text{offset}$ est calibré automatiquement pour
reproduire l'enveloppe vibratoire mesurée par Kolbe ($\sim 10$ MOA/ms en
pic).

![Réponse transitoire simulée pour la configuration nominale (tuner
200 g, $h_\text{offset} = 16{,}5$ mm). De haut en bas : déflexion
$y(L,t)$, angle de bouche $\theta(L,t)$, et vitesse angulaire
$\dot\theta(L,t)$. Le trait rouge pointillé marque $t_b = 2{,}5$ ms ; le
trait vert pointillé sur le panneau inférieur indique la cible Kolbe (6
MOA/ms). La valeur calculée à $t_b$ est de $+6{,}64$ MOA/ms, en accord à
$\sim 10\,\%$ avec la
mesure.](plot_tir_nominal.png){#fig:reponse_transitoire
width="0.85\\linewidth"}

![Balayage paramétrique sur la masse du tuner $m_t \in [0, 400]$ g
(autres paramètres fixes). De haut en bas : (i) $\dot\theta(L,t_b)$
varie peu (5,9 à 6,7 MOA/ms) car $t_b/T_1 \approx 0{,}085$ reste dans le
quart de cycle ascendant initial ; (ii) $\theta(L,t_b)$ varie en
revanche fortement et monotoniquement ($+166$ à $-255$ µrad), traduisant
un décalage du point d'impact moyen ; (iii) $f_1$ décroît
monotoniquement de 39,9 à 30,3 Hz (relation prévue par le quotient de
Rayleigh, équation [\[eq:dwdmt\]](#eq:dwdmt){reference-type="eqref"
reference="eq:dwdmt"}).](plot_balayage_tuner.png){#fig:balayage_tuner
width="0.85\\linewidth"}

Ces tracés illustrent une distinction essentielle pour le tireur : le
tuner agit *principalement sur le point d'impact moyen* (via
$\theta(L,t_b)$) et *secondairement sur la dispersion verticale* (via
$\dot\theta(L,t_b)$). Lorsque $t_b \ll T_1$, la fenêtre d'optimalité en
$\dot\theta$ est large et plate ; le réglage du tuner sert alors surtout
à compenser la dispersion de vitesse, le réglage du visage prenant en
charge le décalage du POI.

# Conclusions pratiques pour le tireur {#sec:conclusions}

Le formalisme précédent, validé numériquement par la simulation par
éléments finis, conduit à plusieurs conclusions exploitables.

#### 1. Effet d'une masse de bouche.

L'ajout d'une masse à la bouche *abaisse* la fréquence fondamentale
(cf. [\[eq:dwdmt\]](#eq:dwdmt){reference-type="eqref"
reference="eq:dwdmt"}). Pour le canon de référence, la plage d'accord
d'un tuner de 100--400 g couvre une variation de $f_1$ de l'ordre de $8$
à $24\,\%$ (39,9 Hz $\to$ 30,3 Hz). Le décalage temporel équivalent dans
le cycle vibratoire est de l'ordre de quelques fractions de
milliseconde, ce qui suffit à balayer l'intégralité de la fenêtre
d'optimalité dans laquelle $t_b \approx 2$--$3$ ms doit être positionné.

#### 2. Sens du réglage.

Allonger la position effective du tuner (le visser vers l'extérieur)
augmente $J_t$ et le bras de levier, ce qui amplifie l'effet de la masse
et abaisse $f_1$. C'est le geste qui « ralentit » la cinématique de
bouche et retarde le passage par zéro de $\theta(L,t)$ : utile pour des
munitions plus lentes (sortie $t_b$ retardée).

#### 3. Sensibilité à la munition.

Tout changement de lot ou de marque modifie $\bar v$ et donc $t_b$. La
fenêtre d'optimalité ($\dot\theta_\text{out}$ proche de
$\dot\theta^\star_\text{out}$) est étroite : un décalage de 0,1 ms sur
$t_b$ peut représenter une fraction non négligeable de la période
$T_1 = 2\pi/\omega_1$. Tout nouveau lot doit donc faire l'objet d'un
retuning systématique au pas de tuner (méthode *ladder tune* ou
*Audette*).

#### 4. Limites physiques.

Le tuner ne supprime pas la dispersion intrinsèque de la munition : il
convertit une dispersion verticale en un groupement réellement plus
serré uniquement si
[\[eq:compensation\]](#eq:compensation){reference-type="eqref"
reference="eq:compensation"} est vérifiée à $\pm$ quelques MOA/ms près.
La dispersion horizontale, indépendante du couple
($\Delta v_0,\, \theta_\text{out}$), n'est pas compensée. Le gain
attendu en cible reste donc borné par les autres sources d'erreur (vent,
plomb de balle, parallaxe, mire).

#### 5. Robustesse aux modes supérieurs.

Les modes 2 et 3, de fréquences nettement plus élevées (cf. tableau
précédent), produisent une oscillation rapide superposée sur
$\theta(L,t)$. S'ils sont mal amortis et excités significativement, ils
créent une variabilité résiduelle (*flyers*) même après tuning. Une
masse à la bouche distincte d'un harmonique entier du mode 1
(anti-résonance) aide à les filtrer, ce qui motive les conceptions «
accordable » à double anneau.

#### 6. Domaine de validité : centerfire *vs* rimfire.

Une réserve d'honnêteté, formulée par Kolbe
lui-même `\cite{KolbeSim}`{=latex}, mérite d'être rapportée. Le
*principe* de compensation positive
(sections [1](#sec:intro){reference-type="ref" reference="sec:intro"}
et [5](#sec:cinematique){reference-type="ref"
reference="sec:cinematique"}) est général et fut d'abord établi pour le
rimfire. En revanche, le *modèle d'excitation par moment impulsionnel de
recul* retenu ici (et dans le simulateur de Kolbe) est avant tout
pertinent pour les calibres *centerfire* : la courbe pression-temps y
est longue ($\sim 1$ ms, pic $\sim 350$ MPa $\approx 50\,000$ psi,
profil type .308 Win) et le recul élevé, si bien que le moment de recul
domine effectivement la mise en vibration. Pour la *.22 LR*, l'impulsion
de pression est si brève que, selon Kolbe, le moment de recul seul
reproduit mal les vibrations observées : d'*autres* sources d'excitation
(engagement dans les rayures, effet de poids du projectile mobile, jeux
de l'assemblage) y prennent une part comparable. Notre modèle les
intègre partiellement --- c'est tout l'objet des contributions (b) et
(c) de la section [3.4](#sec:excitation){reference-type="ref"
reference="sec:excitation"}, absentes du simulateur originel de Kolbe
--- mais les valeurs numériques rimfire avancées dans ce document
doivent être lues comme des *ordres de grandeur calibrés*, non comme des
prédictions absolues. La validation expérimentale directe
(section [5.5](#sec:validation_num){reference-type="ref"
reference="sec:validation_num"}) reste la référence pour le rimfire.

#### 7. Un modèle planaire : vibration verticale seulement.

Une seconde réserve d'honnêteté porte sur la dimensionnalité du modèle.
Tout le développement précédent est *planaire* : la poutre
d'Euler-Bernoulli [\[eq:EB\]](#eq:EB){reference-type="eqref"
reference="eq:EB"} n'est résolue que dans un seul plan et ne décrit que
la composante *verticale* de l'oscillation de bouche. Rien n'oblige
pourtant le canon à ne vibrer que verticalement : en réalité la bouche
décrit une *orbite bidimensionnelle* (un « whip » elliptique), et sa
composante horizontale --- déjà signalée comme non compensée au point 4
--- n'est corrigée par aucun mécanisme de compensation positive. Elle
s'ajoute en dispersion résiduelle et fixe un plancher à la précision
atteignable. Si le modèle vertical reste néanmoins pertinent au premier
ordre, c'est que *deux effets brisent la symétrie* et privilégient le
plan vertical : la *gravité*, qui impose une flèche statique et oriente
le mode dominant ; et la *géométrie du recul*, dont le couple de
relèvement (âme au-dessus du centre de gravité,
section [1](#sec:intro){reference-type="ref" reference="sec:intro"}) est
essentiellement vertical. Surtout, la compensation positive
[\[eq:compensation\]](#eq:compensation){reference-type="eqref"
reference="eq:compensation"} n'exige pas que la vibration soit
*purement* verticale : il suffit que la composante verticale de l'angle
à la bouche corrèle, dans le bon sens, avec le temps de sortie $t_b$
(donc avec la vitesse initiale). Le mouvement horizontal ajoute du bruit
sans détruire ce bénéfice, ce que confirme le succès empirique des
tuners en benchrest .22 LR. Le modèle planaire est donc un modèle de
*premier ordre* : utile et prédictif, mais qui sous-estime
structurellement la dispersion réelle.

#### 8. Un réglage spécifique à la distance.

La vitesse angulaire optimale
[\[eq:dthetaopt\]](#eq:dthetaopt){reference-type="eqref"
reference="eq:dthetaopt"} est *proportionnelle à la distance de tir* :
$\dot\theta_\text{out}^\star = -\,gD/(v_0^3\,\tau_v) \propto D$. Un
tuner accordé pour $D = 50$ m ne réalise donc l'égalité
[\[eq:compensation\]](#eq:compensation){reference-type="eqref"
reference="eq:compensation"} qu'à cette distance ; à une autre distance,
le terme balistique (en
$\partial h_\text{bal}/\partial v_0 \propto D^2$) et le terme angulaire
(en $\dot\theta_\text{out}\,D$) ne s'annulent plus exactement et la
compensation devient *partielle* --- sous-compensation au-delà de la
distance de réglage, sur-compensation en-deçà ---, la dispersion de
vitesse réapparaissant graduellement en traînée verticale. Il faut en
retenir que la compensation positive ne *réduit pas* la dispersion de
vitesse initiale du lot (l'écart-type des vitesses reste inchangé) :
elle en *masque* l'effet vertical, et seulement au voisinage de la
distance de réglage. Un lot à faible écart-type de vitesse reste donc
préférable, et un tuner gagne à être vérifié à la distance de la
compétition. Le phénomène est néanmoins progressif : un réglage établi à
courte distance conserve une part de son bénéfice à plus longue distance
tant que le groupement intrinsèque reste inférieur à la traînée
balistique qu'il corrige.

#### 9. Un tuner n'est pas toujours bénéfique : l'étude MEF de Harral.

Une dernière mise en garde, empirique celle-là, vient d'une simulation
par éléments finis indépendante conduite par
A. Harral `\cite{Harral}`{=latex} sur une carabine benchrest .22 LR
(l'arme d'« Esten »). En calculant, pour la plage de vitesses
1035--1075 fps, la traînée verticale sur cible avec et sans masse de
bouche, Harral obtient un résultat contre-intuitif : le *meilleur*
groupement vertical est ici obtenu *sans* tuner ($0{,}091$ pouce, soit
$\approx 2{,}3$ mm), l'ajout d'une masse le *dégradant* monotonement ---
$0{,}122$ pouce ($\approx 3{,}1$ mm) avec 4,9 oz ($\approx 139$ g),
$0{,}183$ pouce ($\approx 4{,}6$ mm) avec 16 oz ($\approx 454$ g).
L'interprétation est cohérente avec notre formalisme : ce canon nu place
déjà, par sa géométrie (profil *reverse taper* flexible), la sortie de
balle $t_b$ sur un flanc *ascendant* favorable de $\theta(L,t)$ ;
alourdir la bouche ralentit la cinématique
[\[eq:dwdmt\]](#eq:dwdmt){reference-type="eqref" reference="eq:dwdmt"}
et fait dériver $t_b$ vers un flanc *descendant*, où la combinaison
devient additive (dispersion *aggravée*). On retiendra deux leçons
pratiques : (i) le tuner ne *crée* pas la compensation, il ne fait que
*déplacer* $t_b$ dans le cycle vibratoire --- si le canon est déjà bien
placé, la meilleure position de tuner peut être « aucune masse », et son
rôle se réduit à *réaccorder* après un changement de lot ou de distance
(point 3 et point 8) ; (ii) une masse *légère*, à réglage plus fin, est
préférable à une masse lourde qui risque de « survoler » l'optimum d'un
cran à l'autre. Harral confirme par ailleurs, sur la même étude, le
caractère spécifique à la distance (point 8), un réglage à 50 yd ne
restant pas optimal à 100 yd.

*Réserve sur cette mise en garde.* Tout ce qui précède repose sur une
source *unique et purement numérique*. Kolbe `\cite{KolbeSim}`{=latex}
--- qui, lui, a mesuré --- juge ce travail non confirmé :  His work has
lacked the experimental confirmation needed to verify his computer
modelling, however.  Harral ne publie ni les données de son arme (un
seul chiffre : 10,5 lb), ni ses fréquences propres, ni son
amortissement. Le résultat reste plausible et cohérent avec le
formalisme, mais doit être lu comme une *prédiction de simulation non
validée*, non comme un fait mesuré.

#### 10. À l'étau ou à l'épaule ? La condition d'appui *est* le mécanisme.

La question la plus fréquente en pratique --- faut-il accorder l'arme
serrée dans un étau, pour « éliminer le facteur humain » ? --- appelle
une réponse tranchée : *non*, et pas pour une raison de réalisme, mais
parce qu'un étau rigide *supprime le phénomène même que l'on cherche à
régler*. Le moment excitateur
[\[eq:moment_recul\]](#eq:moment_recul){reference-type="eqref"
reference="eq:moment_recul"} n'existe que parce que l'arme *recule et
pivote autour de son centre de gravité*
(section [3.4](#sec:excitation){reference-type="ref"
reference="sec:excitation"}) ; bloquer rigidement la base annule ce
terme, et il ne subsiste que la flexion élastique du tube, qui ne
produit *aucune* compensation. Le banc de Kolbe lui-même était un étau
--- mais il précise qu'il n'était pas rigide :  The relatively thin base
plate flexed under recoil and allowed the barrel clamp to rotate
backwards, resulting in an upwards vertical muzzle flip.  Un serrage
réellement rigide n'aurait rien montré. Pour le *benchrest*, la question
est close : le tir au sac est, selon Kolbe,  a fair approximation  du
recul libre postulé par le modèle. Pour le *tir épaulé ou couché*, en
revanche, épauler ajoute de la masse effective et contraint la rotation
: puisque l'amplitude varie comme $h_\text{cg}/m_r$
[\[eq:hoffset_phys\]](#eq:hoffset_phys){reference-type="eqref"
reference="eq:hoffset_phys"}, un accord établi en recul libre n'a aucune
raison d'être optimal épaulé. Kolbe reste prudent ( if a small calibre
rifle is gripped tightly or pulled hard into the shoulder then the
recoil dynamics could be affected ) et ne l'a pas mesuré. En l'état des
connaissances, la règle de prudence est donc : *accorder dans la
position de tir*, et re-vérifier l'accord si l'on change de position ou
de tenue.

# Bibliographie {#bibliographie .unnumbered}

::: thebibliography
9

S. S. Rao, *Mechanical Vibrations*, Pearson, 2017.

L. Meirovitch, *Elements of Vibration Analysis*, McGraw-Hill, 1986.

K.-J. Bathe, *Finite Element Procedures*, 2nd ed., Prentice Hall, 2014.

D. Carlucci and S. Jacobson, *Ballistics: Theory and Design of Guns and
Ammunition*, 3rd ed., CRC Press, 2018.

G. Kolbe, *Using barrel vibrations to tune a barrel --- The Vibrations
of a Barrel Tuned for Positive Compensation*, Border Barrels, 2015.

G. Kolbe, *The Vibrations of a Barrel Tuned for Positive Compensation*
(article en ligne, mise à jour du 18 novembre 2015),\
<http://www.geoffrey-kolbe.com/articles/rimfire_accuracy/tuning_a_barrel.htm>.

G. Kolbe, *Barrel Vibrations Simulator* (modèle « lumped parameter » par
éléments finis et notes de modélisation, en ligne, consulté en juin
2026),\
<http://www.geoffrey-kolbe.com/articles/rimfire_accuracy/barrel_vibrations.htm>.

A. Mallock, *Vibrations of Rifle Barrels*, Proceedings of the Royal
Society, Vol. 68, p. 327, 1901.\
<https://www.tireur.org/articles/Mall01.pdf>.

C. Audette, *The Optimum Charge Weight (OCW) Method*, *Precision
Shooting Magazine*, 2005--2010.

A. Harral, *Al's 22LR --- Barrel Tuner Analysis* (étude par éléments
finis d'un tuner de bouche sur carabine benchrest .22 LR, en ligne,
consulté en juillet 2026),\
<https://varmintal.com/a22lr.htm>.

N. M. Newmark, *A Method of Computation for Structural Dynamics*, ASCE
Journal of the Engineering Mechanics Division, 85(3), pp. 67--94, 1959.
:::
