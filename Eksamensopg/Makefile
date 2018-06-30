CFLAGS = -Wall -std=gnu1x -O2
LDLIBS = -lm


plot: splines.dat plots.svg final.tex Makefile
	pdflatex final.tex
	pdflatex final.tex
	xdg-open final.pdf

plots.svg: splines.dat Makefile
	gnuplot plots.gpi

splines.dat: main
	./$< > out.txt

main: main.o aspline.o cspline.o binsearch.o quadspline.o ospline.o


.PHONEY:clean

clean:
	find . -type f -executable -delete
	$(RM) *.o *.out *.log *.dat* *.svg *.pdf *.aux *.txt *~
