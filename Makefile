# Chaîne de génération du dossier « tuners » (vibrations de canon).
#
#   Sources   : tuner_fr.tex, tuner_en.tex
#               (+ plot_*.png produits par simulation.jl — voir ce script)
#               harral_a22lr.jl — script de recréation FEA (Julia), servi tel quel
#   Produits  : *.pdf          — dans ce dossier ET dans ../files/tuners (web)
#               tuner_fr.md     — dérivé pandoc (artefact interne, non servi)
#               ../files/tuners/harral_a22lr.jl — copie servie du script
#
#   Dépendances : pdflatex (TeX Live), pandoc.
#
#   NB pandoc : la MISE EN FORME du .md dépend de la version de pandoc.
#   Génération de référence = pandoc 3.1.3. Une autre version reformatera
#   (puces, sauts de ligne) sans changer le fond ; l'option
#   « -t markdown-citations » est indispensable pour préserver le rendu
#   \cite brut, cohérent avec la bibliographie manuelle (thebibliography).
#
#   Cibles usuelles :
#     make            → PDF (2 passes) + tuner_fr.md + copie web
#     make pdf        → seulement les PDF
#     make md         → seulement tuner_fr.md
#     make publish    → recopie les PDF ET le script vers ../files/tuners
#     make clean      → supprime les auxiliaires LaTeX

PDFLATEX = pdflatex -interaction=nonstopmode -halt-on-error
PANDOC   = pandoc -s -t markdown-citations
WEBDIR   = ../files/tuners

DOCS   = tuner_fr tuner_en
PDFS   = $(addsuffix .pdf,$(DOCS))
AUX    = $(foreach d,$(DOCS),$(d).aux $(d).log $(d).out $(d).toc)

# Figures incluses par les deux documents (générées par simulation.jl).
IMAGES = plot_modes_propres.png plot_tir_nominal.png plot_balayage_tuner.png

# Script servi tel quel (lié depuis le wiki) : recréation FEA de Harral.
SCRIPT = harral_a22lr.jl

.PHONY: all pdf md publish clean
all: pdf md publish

pdf: $(PDFS)

# Deux passes : références croisées (\ref/\eqref) + bibliographie manuelle.
%.pdf: %.tex $(IMAGES)
	$(PDFLATEX) $<
	$(PDFLATEX) $<
	$(RM) $*.aux $*.log $*.out $*.toc

md: tuner_fr.md
tuner_fr.md: tuner_fr.tex
	$(PANDOC) $< -o $@

# Copie des PDF et du script vers le dossier servi par le site.
publish: $(PDFS) $(SCRIPT)
	cp $(PDFS) $(SCRIPT) $(WEBDIR)/

clean:
	$(RM) $(AUX)
