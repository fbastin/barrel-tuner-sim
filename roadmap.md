# Feuille de route (Roadmap) — Dynamique et Optimisation du Tuner de Canon

Ce document consigne l'analyse comparative des modélisations du simulateur de tuner par rapport à la littérature scientifique existante, et trace les perspectives d'amélioration à court, moyen et long terme.

> **Origine et mise à jour.** Rédigé avec Gemini le 20 juillet 2026, puis **révisé le même jour** contre l'état courant du dépôt. Les axes de travail proposés restent pertinents et plusieurs sont **déjà engagés** — ils sont marqués comme tels. En revanche les chiffres de validation du §1 et du §2 dataient d'états antérieurs du modèle et ont été corrigés : ils annonçaient un accord avec Kolbe que le modèle **ne réalise pas**.

---

## 1. État de l'art et positionnement scientifique actuel

Les travaux du dossier [tuners](file:///var/www/tireur.org/tuners) se situent au croisement de la dynamique transitoire des structures et de la balistique. Ils confrontent et valident les travaux historiques et contemporains de la discipline :

*   **Principe de compensation positive (A. Mallock, 1901 ; G. Kolbe, 2015)** :
    *   *Théorie* : Exploiter la vibration de flexion du canon lors du recul pour que les projectiles plus lents sortent lorsque la bouche pointe vers le haut ($\dot{\theta} > 0$), compensant ainsi leur flèche balistique différentielle.
    *   *Validation, **corrigée le 20 juillet 2026*** : il faut distinguer deux choses que la rédaction d'origine confondait.
        *   **Ce qui est validé** : la *chaîne cinématique* de Kolbe. En partant de ses grandeurs mesurées (chute naturelle $0{,}016$, sensibilité $\tau_v = 8{,}8$ $\mu$s/(m/s)), la formule $\dot\theta^\star = gD/(v_0^3\tau_v)$ rend $5{,}96$ MOA/ms contre $6{,}0$ mesurés au banc — accord meilleur que 1 %, par deux chemins indépendants. Ce résultat tient (`kolbe_validation.jl`).
        *   **Ce qui n'est PAS validé** : que *notre* simulateur y parvienne. Il prédit $\tau_v = 7{,}6$ $\mu$s/(m/s), d'où une cible auto-cohérente de **$6{,}88$** MOA/ms — et son $\dot\theta$ **plafonne à $6{,}38$**. La compensation n'est atteinte à aucun réglage. L'accord affiché auparavant venait de ce que la cible était **codée en dur à $6{,}0$**, empruntée à Kolbe : le modèle se comparait à une cible qui n'était pas la sienne.
        *   *Précision utile* : $\tau_v$ **dépend de la longueur du canon** ($7{,}07$ $\mu$s/(m/s) à 610 mm, $8{,}74$ à 762 mm). Ce n'est pas une constante de cartouche, et la longueur sur laquelle Kolbe a mesuré $8{,}8$ n'est documentée dans aucune de nos sources — l'écart de « 13 % » longtemps publié n'est donc pas établi.
*   **Variabilité de l'excitation (H. R. Vaughn, 1998)** :
    *   *Théorie* : Le moment de recul varie de $\pm 30\ \%$ d'un coup à l'autre en raison de la gigue de combustion et des asymétries mécaniques de la culasse.
    *   *Validation* : Notre script [variability.jl](file:///var/www/tireur.org/tuners/variability.jl) simule par Monte-Carlo cet effet aléatoire, prouvant que cette variabilité fixe un plancher de dispersion verticale résiduelle infranchissable par le seul accord du tuner.
*   **Limites des modèles encastrés et simplifiés (A. Harral / Varmint Al)** :
    *   *Analyse critique* : Nous avons montré que les modèles purement encastrés à la culasse ne peuvent pas restituer la compensation positive, car le moment excitateur dépend de la rotation de corps rigide de l'arme sous le recul (modélisé dans [rifle_coupled.jl](file:///var/www/tireur.org/tuners/rifle_coupled.jl)). De plus, l'analyse de Harral confond l'affaissement statique sous gravité du tuner (sag) avec la déflexion dynamique du tir.

---

## 2. Découvertes récentes et limites structurelles du modèle

L'optimisation récente du modèle de balistique interne couplée a mis en évidence une contradiction fondamentale (détaillée dans [ROADMAP.md](file:///var/www/tireur.org/ROADMAP.md)) :

1.  **L'argument d'impossibilité sur le mode de whip fondamental** — *(conclusion confirmée, chiffres corrigés)* :
    *   Les valeurs citées ($-9{,}4$ MOA/ms au canon nu) datent d'avant la correction de `J_tuner` : le canon nu donne aujourd'hui $-2{,}01$ MOA/ms (et $-0{,}32$ avec le tuner de 200 g posé à la bouche, sans porte-à-faux --- les deux avaient été confondus ici). Le raisonnement en déphasage a été refait sur une base plus solide.
    *   **Forme retenue de l'argument** : reproduire le basculement mesuré par Kolbe par un simple décalage de fréquence exigerait $\Delta f/f \approx 72\ \%$, quand un tuner de 200 g ne déplace le fondamental que de $10{,}3\ \%$ — un **facteur 7** d'écart.
    *   *Conclusion inchangée, et c'est le résultat le plus robuste du dossier* : **l'accord ne s'effectue pas sur le mode de flexion fondamental**, mais sur des modes supérieurs (modes 9 et 10, à $2333$ et $3009$ Hz), où le projectile sort après $6$ à $8$ cycles. Deux raisons indépendantes convergent : l'accélération de bouche favorise les hautes fréquences en $\omega^2$, et ces modes portent $28\ \%$ du débattement chacun contre $2\ \%$ pour le fondamental.
2.  **La flexibilité propre du tuner-tube** :
    *   Dans la bande spectrale utile (2-3 kHz), les tuners à tube long (type Starik de 15-30 cm) ne se comportent plus comme des masses ponctuelles rigides car ils présentent leurs propres résonances élastiques internes (vers 600 Hz). Cet effet est désormais **modélisé et quantifié** (poutre élastique dans `frame_v2.jl`/`sensitivity_coupled.jl`), et non plus seulement signalé.
3.  **Le « mur » d'excitation était en grande partie un artefact de l'encastrement** — *(découverte du 21 juillet 2026)* :
    *   Sur le modèle **encastré** (`sensitivity_map.jl`), atteindre la cible de compensation exigeait un bras de levier de **$3{,}4$ à $13{,}8\times$** la cote physique, sur tout le domaine (longueur × amortissement) : un « mur » apparemment structurel.
    *   Mais ce modèle **supprime la rotation d'ensemble de l'arme sous le recul**, mécanisme que Kolbe désigne comme la source de la vibration. Un modèle qui la **restaure** (`sensitivity_coupled.jl` : crosse, appuis, balistique couplée, tube élastique) fait tomber ce facteur à **$1{,}2$ à $2{,}5\times$** selon l'amortissement — un **quasi-manque**, non un mur.
    *   **Conclusion qualitative inchangée** : la compensation n'est atteinte à **aucun réglage physique**, dans aucune des deux versions. Mais sa **magnitude** change du tout au tout, et l'écart de $1{,}2\times$ du meilleur cas tient dans l'incertitude du modèle (l'*ordre* de Kolbe est reproduit, jamais son couple exact $-9{,}4/+6{,}0$). Cette nuance est propagée à tous les supports publiés (wiki, PDF FR/EN, simulateur web).

---

## 3. Feuille de route technique (Prochaines étapes)

Pour lever ces verrous physiques, le simulateur doit évoluer selon les axes suivants :

### Étape 1 : Résolution et maillage haute fréquence ($> 2$ kHz) — **engagée**
*   **Objectif** : Capturer les modes supérieurs ($n \ge 5$) sur lesquels s'opère le réglage fin observé sur le terrain (ladder tuning).
*   **Action** : Affiner le maillage éléments finis du canon et réduire drastiquement le pas de temps d'intégration $\Delta t$ dans le solveur de Newmark-$\beta$.
*   *État* : les modes supérieurs sont **identifiés et quantifiés** (9 et 10, $2333$/$3009$ Hz, $28\ \%$ du débattement chacun). Réserve connue : à 20 éléments il ne reste que $7{,}7$ à $8{,}8$ éléments par longueur d'onde dans cette bande — c'est peu, et le raffinement reste à faire.

### Étape 2 : Modélisation élastique du tuner-tube — **réalisée et composée dans le cadre unifié**
*   **Objectif** : Remplacer l'approximation de masse ponctuelle rigide pour les tuners longs.
*   **Action** : Représenter le tuner comme un segment de poutre flexible doté de ses propres propriétés mécaniques (module d'Young, géométrie creuse, masse linéique).
*   *État (mis à jour le 21 juillet 2026)* : le tube élastique (maillé en poutre, curseur coulissant, aluminium) vit dans `flexible_tube.jl`, dans la branche `:tube` de `build_rifle`, et surtout dans **`frame_v2.jl`, le cadre unifié** — structure fusil × balistique couplée × tube — qui fait désormais tourner les **DEUX idéalisations du tuner côte à côte** (masse ponctuelle rigide *et* tube élastique, à masse totale égale). La de-autorité que le tube impose a été **quantifiée sur balayage de géométrie** (`sensitivity_coupled.jl`) : au meilleur cas, $\dot\theta_\text{max}$ passe de $8{,}7$ (rigide) à $6{,}1$ (tube).
*   ⚠️ **Correction d'un chiffre trop fort** : le « l'idéalisation rigide gonflait l'effet d'un **ordre de grandeur** » mesuré autrefois valait sur le modèle **encastré à géométrie figée**. Dès qu'on libère la géométrie d'arme dans le modèle couplé, la de-autorité tombe à **~×1,2 à 1,5** — la liberté de configuration récupère l'essentiel de l'écart. Corollaire de méthode, matérialisé par un garde-fou explicite dans `frame_v2` : à géométrie **figée**, l'ordre rigide/tube n'est pas robuste et peut même s'inverser ; seule une carte à géométrie balayée quantifie.
*   *Relation fixée le 21 juillet 2026 — pas un remplacement, une division des rôles.* `frame_v2` est promu **modèle unifié de référence** au sens de la physique la plus complète, et sert de **validateur qui borne l'erreur** de `simulation.jl`. Il ne le **remplace pas** : `simulation.jl`, encastré et paramétré, reste la source des grandeurs **publiées**, parce que ses chiffres sont propres, reproductibles et indépendants de toute géométrie d'arme non mesurée — tandis que ceux de `frame_v2` dépendent d'une géométrie inconnue et ne sont robustes qu'en **balayage** (`sensitivity_coupled.jl`), pas à géométrie figée. **Référence de physique ≠ référence de chiffres** : les deux rôles sont distincts et assumés. L'idée initiale d'un remplacement pur et simple est donc **abandonnée** comme mal posée.

### Étape 3 : Calibration expérimentale de l'amortissement modal ($\zeta$) — **protocole publié, mesure à faire (priorité 1)**
*   **Objectif** : Remplacer le modèle d'amortissement de Rayleigh qui sur-amortit artificiellement les hautes fréquences.
*   **Action** : Permettre l'attribution d'un coefficient d'amortissement modal indépendant $\zeta_n$ pour les modes élevés. Ces coefficients devront être mesurés expérimentalement par analyse modale physique (choc au marteau instrumenté et accéléromètre de bouche sur canon réel).
*   *État* : **protocole rédigé et publié** ([protocole_amortissement](https://tireur.org/wiki/doku.php?id=technique:protocole_amortissement)). Reste un verrou majeur — l'amortissement gouverne la réachabilité — mais son statut a été **nuancé le 21 juillet 2026** par la carte couplée : même au $\zeta$ **le plus favorable** ($0{,}3\ \%$), le modèle couplé avec tube élastique reste court d'un facteur $1{,}2\times$. **$\zeta$ seul n'ouvre donc pas la compensation** ; il faudrait qu'un autre facteur (modes supérieurs mal maillés — Étape 1 — ou appui de crosse) joue *aussi* dans le bon sens. Mesurer $\zeta$ dit dans quelle région de la carte se trouve une arme réelle, sans garantir à lui seul l'accès à la compensation. Aucune mesure publiée ne couvre la bande 2–3 kHz — l'obstacle n'est pas l'instrumentation mais la *grandeur choisie* par les auteurs (caméra à $0{,}0177$ mm/pixel pour une amplitude de $3\ \mu$m, analyseur à 800 Hz).

### Étape 4 : Modélisation 3D de la dispersion
*   **Objectif** : Simuler la dérive horizontale et la dispersion radiale en cible.
*   **Action** : Étendre le solveur aux trois dimensions (flexion transversale $y$ et $z$, ainsi qu'au couplage de torsion) pour modéliser les asymétries d'excitation tridimensionnelles.

---

## 4. Littérature scientifique complémentaire à explorer

Pour appuyer ces développements, l'intégration des travaux académiques et industriels suivants est recommandée :

> *État du dépouillement (20 juillet 2026)* : les points 1 et 2 sont **entamés** — le rapport ARCCB-TR-02002 (*Gun Barrel Vibration Absorber*, domaine public) et le brevet US 5 798 473 sont archivés dans `articles/` et exploités ; la thèse d'O'Neil a été lue intégralement et a fourni la raideur de racine élastique ($K_\text{root}$) qui reproduit son rapport $f_2/f_1 = 8{,}16$ mesuré. Les points 3 et 4 restent **entièrement ouverts**, et le point 3 est le plus prometteur : c'est la seule voie qui validerait les amplitudes absolues sans passer par un calage arbitraire — donc la seule qui puisse trancher sur $h_\text{offset}$.

1.  **Dynamique des canons de chars et amortisseurs de bouche (Tuned Vibration Absorbers - TVA)** :
    *   La littérature militaire (notamment du *US Army Research Laboratory*) modélise de façon très rigoureuse les dispositifs amortisseurs de bouche sur les canons de gros calibre, équivalents physiques à grande échelle de nos tuners.
    *   *Utilité* : Fournir des formulations analytiques d'optimisation de placement de masse et de matériaux amortisseurs.
2.  **Étude des raideurs d'interface et d'encastrement** :
    *   Les publications sur le comportement dynamique des assemblages filetés boîtier/canon et de l'amortissement introduit par le couchage (*bedding*) dans la crosse.
    *   *Utilité* : Remplacer les conditions aux limites trop simples (encastrement rigide ou recul 100 % libre) par des ressorts de raideur équivalente à l'interface.
3.  **Vibrométrie Laser Doppler et DIC (Digital Image Correlation)** :
    *   Les études académiques mesurant directement le déplacement et la vitesse angulaire réelle de la bouche lors du tir à l'aide de capteurs optiques ou laser.
    *   *Utilité* : Valider les amplitudes absolues de la simulation sans dépendre d'hypothèses de calage arbitraires.
4.  **Confrontation avec la théorie OBT (Optimal Barrel Time)** :
    *   La théorie OBT (C. Long) postule que la régularité du tir dépend du passage d'ondes de pression *longitudinales* (vitesse du son dans l'acier $\approx 5100$ m/s) à la bouche.
    *   *Utilité* : Évaluer si un couplage physique existe entre la contraction/dilatation longitudinale du canon et son angle transversal de flexion à l'instant de sortie du projectile.
