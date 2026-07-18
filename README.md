# barrel-tuner-sim

Simulation par éléments finis des **vibrations transversales d'un canon de
carabine**, et analyse de la **compensation positive** et de l'accord par
*tuner* (masse de bouche).

Modèle : poutre d'Euler–Bernoulli (éléments d'Hermite cubiques, 2 d.d.l. par
nœud), intégration temporelle de Newmark-β, excitation par le moment de recul
et par la charge mobile du projectile. Le mécanisme physique suit
G. Kolbe : les vibrations naissent du **recul, qui fait pivoter l'arme autour
de son centre de gravité** et imprime un moment à la culasse.

## Contenu

| Fichier | Rôle |
| --- | --- |
| `simulation.jl` | Simulation principale : modes propres, tir nominal, balayage de masse de tuner. Produit les figures `plot_*.png`. |
| `harral_a22lr.jl` | Recréation de l'analyse FEA de A. Harral (canon *reverse taper* .22 LR d'« Esten »), à section variable, culasse encastrée. |
| `harral_rifle_sweep.jl` | Modèle du **fusil entier libre** sur deux sacs (canon + boîte + crosse + avant-bras), appuis unilatéraux, balayage géométrique. |
| `kolbe_validation.jl` | Validation de la chaîne cinématique de Kolbe (chute naturelle, sensibilité du temps de sortie, taux d'accord de 6,0 MOA/ms à 50 m). |
| `kolbe_amplitude.jl` | Amplitude du taux angulaire de bouche : le fusil libre atteint-il l'ordre de grandeur mesuré par Kolbe ? |
| `variability.jl` | Dispersion prédite quand l'excitation ne se répète pas d'un coup à l'autre (±30 % mesurés par Vaughn) : Monte-Carlo sur vitesse **et** amplitude, plancher non compensable, critère d'accord révisé. |
| `tuner_fr.tex` / `tuner_en.tex` | Document explicatif (physique, résultats, bibliographie) — versions française et anglaise. |
| `v1.tex`, `v2.tex`, `fig1.tex`, `code.jl.txt` | Versions antérieures, conservées pour l'historique. |
| `Makefile` | Génère les PDF (pdflatex, 2 passes) et le `.md` dérivé (pandoc). |

## Exécution

```sh
julia simulation.jl        # figures plot_*.png + tir nominal + balayage
julia harral_a22lr.jl      # table de dispersion façon Harral (culasse encastrée)
julia kolbe_validation.jl  # chaîne cinématique de Kolbe
julia kolbe_amplitude.jl   # amplitude du taux angulaire (fusil libre)
julia variability.jl     # dispersion prédite avec excitation variable (Monte-Carlo)
julia harral_rifle_sweep.jl  # balayage du fusil libre (long)
```

Dépendances Julia : `LinearAlgebra` et `Printf` (bibliothèque standard).
`simulation.jl` installe `Plots.jl` à la première exécution ; les autres
scripts s'en passent si le paquet est absent (les tableaux sont affichés, les
tracés ignorés).

Documents :

```sh
make pdf    # tuner_fr.pdf + tuner_en.pdf   (pdflatex)
make md     # tuner_fr.md                   (pandoc 3.1.3 de référence)
make clean  # auxiliaires LaTeX
```

## Ce que ce code reproduit — et ce qu'il ne reproduit pas

Ce simulateur vise les **ordres de grandeur** et la **structure physique**, pas
la reproduction au centième d'études dont les données ne sont pas publiées.

**Reproduit :**
- La **chaîne cinématique de Kolbe** : chute naturelle 0,016 MOA/ft/s à 50 m,
  sensibilité du temps de sortie ≈ 375 ft/s/ms, taux d'accord 6,0 MOA/ms
  (`kolbe_validation.jl`), à ~1 %.
- L'**ordre de grandeur de l'amplitude** du taux angulaire de bouche : le
  fusil libre, avec un moment de recul physiquement ancré, atteint un pic de
  quelques MOA/ms, celui de Kolbe (`kolbe_amplitude.jl`).
- Le rôle **décisif de la condition d'appui** : culasse encastrée → aucune
  compensation ; arme libre → dynamique de bouche radicalement différente.

**Ne reproduit pas :**
- La **table de dispersion de Harral** ni les **amplitudes absolues de Kolbe**
  au centième : les deux dépendent de données non publiées (géométrie exacte
  de l'arme, souplesse du banc d'essai, amortissement). La dispersion prédite
  est notamment très sensible à l'amortissement du canon.
- La compensation est par ailleurs **spécifique à la distance** — elle masque
  l'écart de vitesse autour d'une distance, elle ne l'élimine pas.

Le modèle est **planaire** (vibration verticale seule) : la composante
horizontale du mouvement de bouche fixe un plancher de dispersion résiduelle
qu'il ne décrit pas.

## Sources

- A. Harral, « 22 Long Rifle Barrel Tuner Analysis — FEA Dynamic Analysis »,
  <https://varmintal.com/a22lr.htm>
- G. Kolbe, « Using barrel vibrations to tune a barrel »,
  <http://www.geoffrey-kolbe.com/articles/rimfire_accuracy/tuning_a_barrel.htm>
- A. Mallock, « The Recoil of Guns », *Proc. Roy. Soc.*, 1901.

Les bibliographies complètes figurent dans `tuner_fr.tex` / `tuner_en.tex`.

## Licence

- **Code** (`*.jl`, `Makefile`) : MIT — voir [`LICENSE`](LICENSE).
- **Documentation et figures** (`*.tex`, `*.md`, `*.pdf`, `*.png`) :
  CC BY-SA 4.0 — voir [`LICENSE-docs`](LICENSE-docs).

© 2026 Fabian Bastin. Extrait du site [tireur.org](https://tireur.org).
