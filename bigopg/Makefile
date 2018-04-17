CFLAGS = -Wall -std=gnu1x -lm $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)


tex: plot BP.tex
	pdflatex BP.tex
	$$(echo xdg-open BP.pdf)

plot: plot.gpi data.txt
	gnuplot $<


data.txt: main
	./$< > data.txt


main: main.c


.PHONEY: clean
clean:
	rm main data.txt *.pdf *.log *.out
