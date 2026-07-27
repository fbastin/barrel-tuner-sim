# Évaluation Critique — Modélisation et Simulation du Tuner de Canon

Ce document dresse un bilan critique, objectif et honnête du modèle physique de tuner développé dans le dossier [tuners](file:///var/www/tireur.org/tuners). Il distingue les réussites conceptuelles majeures des simplifications ou arrangements méthodologiques qui limitent sa portée prédictive.

> **Origine et mise à jour.** Rédigé avec Gemini le 20 juillet 2026, puis **révisé le même jour** contre l'état courant du dépôt. Les passages marqués *(révisé)* ont été corrigés : le diagnostic de fond était juste, mais plusieurs chiffres dataient d'avant la refonte en balistique intérieure couplée et d'avant le passage à une cible de compensation dérivée. Ces corrections-là **durcissaient** la critique d'origine. Une seconde vague, marquée *(21 juillet 2026)*, va en sens inverse sur un seul point : le **modèle couplé** (rotation d'ensemble de l'arme + tube élastique) montre que le « mur » de $3{,}4$–$13{,}8\times$ était en grande partie un artefact de l'encastrement, et ramène le manque à un quasi-écart de $1{,}2$–$2{,}5\times$. Elle **atténue la magnitude** sans toucher à la conclusion qualitative — le modèle ne reproduit toujours la compensation à aucun réglage physique.

---

## 1. Les Points Forts Majeurs (Ce qui est exceptionnel)

Le simulateur et sa documentation atteignent un niveau de rigueur technique rarement observé dans les publications amateur et de vulgarisation sportive :

*   **Rupture avec l'empirisme dogmatique** : Le modèle substitue des équations différentielles solides (poutre d'Euler-Bernoulli, intégration temporelle transitoire de Newmark-$\beta$) aux théories intuitives mais souvent fausses d'ondes stationnaires acoustiques qui polluent la vulgarisation sur les tuners.
*   **Intégrité méthodologique et déconstruction** : L'analyse critique et la tentative de reproduction des travaux de Harral (Varmint Al) dans [harral_a22lr.jl](file:///var/www/tireur.org/tuners/harral_a22lr.jl) démontrent une véritable démarche scientifique. L'auteur identifie et documente des erreurs conceptuelles persistantes dans la littérature (comme la confusion entre la flèche statique du canon sous gravité et l'amplitude de vibration dynamique).
*   **Élégance de l'implémentation logicielle** : Le choix de Julia permet d'allier simplicité de lecture et performances. L'exploitation de la linéarité du système pour décomposer la réponse en composantes (recul vs projectile) permet de réaliser des analyses de sensibilité et des simulations de Monte-Carlo de 20 000 tirs en quelques secondes.
*   **Intégration de la variabilité (Monte-Carlo)** : L'implémentation de la gigue d'excitation de Vaughn ($\pm 30\ \%$ sur le couple de recul) dans [variability.jl](file:///var/www/tireur.org/tuners/variability.jl) constitue le pont le plus réaliste entre la physique idéale de la simulation et la dispersion statistique observée au stand par les tireurs.

---

## 2. Les Faiblesses Méthodologiques (Ce qui pèche ou reste "bricolé")

Pour appliquer les standards d'une publication scientifique rigoureuse, plusieurs approximations et simplifications doivent être relevées :

### A. Le piège du calibrage *ad-hoc* (Curve Fitting) — *(révisé, et aggravé)*
Le diagnostic était juste, et la suite lui a donné raison au-delà de ce qu'il annonçait.

*   Le modèle encastré utilise un bras de levier d'excitation artificiel ($h_\text{offset} = 62{,}5$ mm dans l'état courant) gonflé d'environ **6×** par rapport à la cote physique ($\sim 9$–$12$ mm, seule valeur adossée à une mesure, celle de Vaughn).
*   Ce facteur sert de variable d'ajustement globale pour compenser ce que le modèle omet : les degrés de liberté de corps rigide et la rotation d'ensemble de l'arme sous le recul.
*   **Ce que la critique n'avait pas vu, et qui est pire.** La « concordance à moins de 1 % avec la cible de Kolbe » n'était pas seulement construite : elle était **circulaire**. La cible de $6{,}0$ MOA/ms était codée en dur, reprise de Kolbe, alors qu'elle se dérive — $\dot\theta^\star = gD/(v_0^3\tau_v)$ — à partir d'un $\tau_v$ que *ce* modèle prédit à $7{,}6$ $\mu$s/(m/s) quand Kolbe en mesure $8{,}8$. Le modèle calculait donc $\dot\theta$ avec sa propre cinématique tout en le comparant à la cible d'une autre.
*   **Cible rendue cohérente le 20 juillet 2026 : elle passe à $6{,}88$ MOA/ms, et $\dot\theta$ plafonne à $6{,}38$.** La compensation n'est donc atteinte **à aucun réglage** — malgré le facteur 6 sur le bras de levier. Le « résultat construit » ne tenait pas seulement à un calage forcé : une fois le calage rendu cohérent, il ne tient plus du tout.
*   **Et la critique se généralise — c'est le point le plus solide du dossier.** Le verdict ayant changé trois fois en deux jours au gré de paramètres non mesurés, ils ont été cartographiés (`sensitivity_map.jl`) sur $4$ longueurs de canon $\times$ $5$ valeurs d'amortissement. Comme $\dot\theta$ est *linéaire* en $h_\text{offset}$, la question s'inverse : quel bras de levier faudrait-il pour atteindre la cible ? **Réponse : de $3{,}4$ à $13{,}8$ fois la cote physique, sur tout le domaine.** Il n'existe donc **aucune** combinaison $(L, \zeta)$ où la compensation soit atteinte avec une excitation défendable. Le calibrage *ad-hoc* n'est pas un raccourci qu'un meilleur réglage lèverait : il **compense un mécanisme absent**.
*   **Précision du 21 juillet 2026 — la magnitude, elle, était surestimée.** Ce « $3{,}4$–$13{,}8\times$ » est le chiffre du modèle **encastré**. Un modèle qui *restaure* le mécanisme absent — la rotation d'ensemble de l'arme — et traite le tube du tuner en poutre élastique (`sensitivity_coupled.jl`, composé dans le cadre unifié `frame_v2.jl`) ramène le manque à un **quasi-écart de $1{,}2$ à $2{,}5\times$** selon l'amortissement. La critique tient donc dans sa *direction* — le calage compense bel et bien un mécanisme absent — mais le « mur » était en grande partie un **artefact de l'encastrement**, et l'écart résiduel du meilleur cas ($1{,}2\times$) tient dans l'incertitude du modèle. C'est la nuance, pas le renversement : la compensation reste inatteinte à tout réglage physique, dans les deux versions.

### B. L'omission des vibrations bidimensionnelles (3D)
Le modèle est strictement planaire (flexion verticale uniquement).
*   Dans la réalité, la dispersion en cible est circulaire ou elliptique. Le couple gyroscopique de la balle dans les rayures et les asymétries des tenons de culasse excitent des modes horizontaux et torsionnels.
*   Optimiser le tuner uniquement sur le plan vertical risque de masquer ou d'aggraver la dispersion sur l'axe horizontal.

### C. Le choix et la sensibilité de l'amortissement
L'amortissement de Rayleigh ($[C] = \alpha_M [M] + \alpha_K [K]$) est appliqué par simplicité numérique :
*   Ce modèle est notoirement connu pour sur-amortir les hautes fréquences.
*   Puisque l'accord réel du tuner s'effectue sur des harmoniques élevés ($> 2$ kHz) et non sur le fondamental (35 Hz), l'amplitude des modes critiques à la sortie du projectile dépend d'un paramètre d'amortissement très mal contraint dans le simulateur.
*   *(révisé)* La sensibilité a depuis été **chiffrée**, et elle est pire qu'une simple imprécision : dans le modèle encastré, à $\zeta = 1\ \%$ la compensation est approchable, à $2\ \%$ elle ne l'est plus, et **deux des quatre estimations publiées sont au-dessus de ce seuil**. Le paramètre ne dégrade pas le résultat, il en décide. Un protocole de mesure au stand a été rédigé et publié ([protocole_amortissement](https://tireur.org/wiki/doku.php?id=technique:protocole_amortissement)) ; **la mesure reste à faire**, et c'est un verrou majeur.
*   *(nuancé le 21 juillet 2026)* Le modèle **couplé** avec tube élastique tempère toutefois le rôle de $\zeta$ : même au $\zeta$ **le plus favorable** ($0{,}3\ \%$), il reste court d'un facteur $1{,}2\times$. **$\zeta$ seul n'ouvre donc pas la compensation** — le mesurer dit dans quelle région de la carte se trouve une arme réelle, sans garantir à lui seul l'accès à la compensation ; il faudrait qu'un autre facteur (modes supérieurs mal maillés, appui de crosse) joue *aussi* dans le bon sens.

### D. L'incohérence géométrique du tuner-tube — *(confirmée, partiellement traitée)*
Le tuner est modélisé dans [simulation.jl](file:///var/www/tireur.org/tuners/simulation.jl) comme une masse ponctuelle rigide avec inertie concentrée à la bouche.
*   Pour des tubes longs de type Starik ou Centra (15 à 30 cm), le tube présente des résonances propres dès 600 Hz.
*   Considérer ce tube comme un solide rigide dans la bande spectrale active de l'accord (2-3 kHz) est une incohérence physique.
*   *(révisé)* Critique **fondée et retenue**. Les documents publiés l'affirmaient encore à l'envers — ils déduisaient la rigidité de ce que 600 Hz est « bien au-dessus du fondamental (35 Hz) », alors que 600 Hz est *en dessous* de la bande qui gouverne l'accord, ce qui inverse la conclusion. Corrigé dans le wiki puis dans les PDF le 20 juillet 2026. Un modèle de tube **élastique** existe (`flexible_tube.jl`), et le cadre unifié `frame_v2.jl` fait désormais tourner les **deux idéalisations du tuner côte à côte** — masse ponctuelle rigide *et* tube élastique, à masse totale égale.
*   *(mis à jour le 21 juillet 2026)* Le « l'idéalisation rigide gonflait l'effet d'un **ordre de grandeur** » valait sur le modèle *encastré à géométrie figée*. La de-autorité du tube a depuis été **quantifiée sur balayage de géométrie** (`sensitivity_coupled.jl` : au meilleur cas $\dot\theta_\text{max}$ passe de $8{,}7$ rigide à $6{,}1$ tube) : dès qu'on libère la géométrie d'arme, elle tombe à **~$\times 1{,}2$–$1{,}5$**, la liberté de configuration récupérant l'essentiel de l'écart. Garde-fou de méthode : à géométrie **figée** l'ordre rigide/tube n'est pas robuste et peut s'inverser (matérialisé dans `frame_v2`). *Relation fixée le 21 juillet 2026* : `frame_v2` est le **modèle unifié de référence** (physique la plus complète) et le **validateur qui borne l'erreur** de `simulation.jl` — mais il ne le **remplace pas** comme source des chiffres publiés. Ses nombres dépendent d'une géométrie d'arme inconnue et ne sont robustes qu'en balayage ; ceux du modèle encastré, paramétrés et reproductibles, restent la référence chiffrée. **Référence de physique ≠ référence de chiffres.**

---

## 3. Conclusion générale

Le travail actuel est un **excellent modèle didactique de démonstration de concept**, d'une qualité largement supérieure à tout ce qui est publié en ligne à destination du public des tireurs. 

Cependant, il ne s'agit pas encore d'un modèle prédictif autonome. Sa réussite repose sur des variables d'ajustement d'excitation ($h_\text{offset}$) et d'amortissement calées sur des résultats empiriques connus. Pour passer du statut de démonstrateur à celui d'outil de conception, le simulateur doit s'affranchir de ces raccourcis en adoptant des modèles d'amortissement et de liaisons élastiques validés expérimentalement.

**Mise à jour du 20 juillet 2026 — la conclusion doit être durcie.** Tant que la cible restait empruntée à Kolbe, on pouvait dire que le modèle « reproduisait » la compensation moyennant un calage. Ce n'est plus exact : rendue auto-cohérente, la cible n'est atteinte à aucun réglage. La formulation juste est donc : **ce modèle ne sait pas encore reproduire la compensation positive sans être forcé.** Ce qui ne remet pas en cause le phénomène — Kolbe le mesure au banc, et deux méthodes indépendantes se recoupent chez lui — mais situe exactement ce qui manque, et c'est utile : un mécanisme d'excitation qui n'a pas besoin d'un bras de levier six fois trop grand.

**Complément du même jour — de « pas encore » à « pas ainsi ».** La carte de sensibilité ci-dessus permet de resserrer une dernière fois. « Pas encore » suggérait qu'un meilleur réglage finirait par y arriver ; le balayage montre que non, sur tout le domaine plausible. La conclusion n'est donc plus suspendue à des paramètres à mesurer : **le mécanisme d'excitation retenu — moment de recul appliqué à une poutre encastrée — ne peut pas produire la compensation, quel que soit son réglage.** Ce qui manque est identifié : la rotation de corps rigide de l'arme, hors de portée d'un modèle encastré.

**Mise à jour du 21 juillet 2026 — le mécanisme manquant a été restauré, et il confirme le « pas ainsi » en corrigeant sa magnitude.** Ce que `rifle_coupled.jl` puis `sensitivity_coupled.jl` (cadre unifié `frame_v2.jl`) explorent n'est plus une piste mais un résultat : rendre au modèle la rotation d'ensemble de l'arme — et l'élasticité du tube du tuner — fait tomber le manque du « mur » encastré ($3{,}4$–$13{,}8\times$) à un **quasi-écart de $1{,}2$ à $2{,}5\times$** selon l'amortissement. Le « mur » était donc en grande partie un artefact de l'encastrement. **Mais le « pas ainsi » survit au modèle le plus riche** : la cible n'est atteinte à aucun réglage physique dans aucune des deux versions. La conclusion juste n'est donc ni « mur structurel » ni « compensation atteinte », mais : le modèle est **robustement, modestement court** — à $1{,}2\times$ près au meilleur cas, écart qui tient dans son incertitude propre (l'*ordre* de Kolbe est reproduit, jamais son couple exact $-9{,}4/+6{,}0$).

Il faut noter la conséquence méthodologique, moins confortable. Deux questions circulaient ici sous le même nom :

*   *« Le modèle atteint-il la cible à $h_\text{offset} = 62{,}5$ mm ? »* — indécidable (le verdict bascule sur $L$ et sur $\zeta$) **et circulaire**, $62{,}5$ étant le nombre calé pour que ce soit le cas.
*   *« Le modèle reproduit-il la compensation avec une excitation physique ? »* — tranchée **négativement sur le modèle encastré** ; le modèle couplé en abaisse fortement l'enjeu (quasi-écart $1{,}2$–$2{,}5\times$) sans l'annuler.

La première a longtemps tenu lieu de résultat. Ce n'en était pas un, et c'est elle qui explique les trois retournements de conclusion en deux jours : une question circulaire donne la réponse qu'on lui met dedans.
