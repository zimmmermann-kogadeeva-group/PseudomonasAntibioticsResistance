PNGS=$(patsubst %,figure_%.png,1 2 3 4 5) $(patsubst %,si_figure_%.png,1 2)

all: ${PNGS}

%.png: %.svg
	flatpak run org.inkscape.Inkscape $< --export-background-opacity=1.0 --export-dpi=300 -o $@

clean:
	rm -f ${PNGS}
