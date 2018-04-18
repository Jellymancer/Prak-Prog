'CFLAGS = -Wall -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)

all: erf.pdf
	xdg-open erf.pdf

erf.pdf: final.pdf erf.tex final.pdf
	pdflatex erf.tex


final.pdf: runinteg plot.gpi
	gnuplot plot.gpi


runinteg: integ
	echo -8 8 0.01 | ./integ


integ: integ.o


.PHONEY: clean

clean:
	-rm *.txt final.pdf *.o *.log *.out integ erf.pdf
