 
all: $(patsubst %,Figures/figure_%.png,1 2 3 4 5) $(patsubst %,Figures/si_figure_%.png,1 2)

.PHONY: update_figures
update_figures:
	cp Output/Figures/* Figures/

Figures/%.png: Figures/%.svg update_figures
	flatpak run org.inkscape.Inkscape $< --export-background-opacity=1.0 --export-dpi=300 -o $@
