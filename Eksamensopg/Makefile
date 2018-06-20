''CFLAGS = -Wall -std=gnu1x -O2
LDLIBS = -lm


plot: splines.dat splines.svg dsplines.svg d2splines.svg final.tex Makefile
	pdflatex final.tex
	xdg-open final.pdf

splines.svg: splines.dat Makefile
	echo "\
        set term pdf;\
        set out 'splines.pdf';\
        set xlabel 'x';\
        set ylabel 'y';\
        set key right bottom;\
        set tics out;\
        set title 'Fig. 1: Discrete data points with akima spline';\
        plot 'out.txt' using 1:2 with points title 'discrete points'\
        ,'aspline.dat' using 1:2 with lines title 'aspline'\
        ,'aspline.dat' using 1:3 with lines title 'cspline'\
        ,'aspline.dat' using 1:4 with lines title 'ospline'\
        " | tee log.txt | gnuplot
	sed -i 's/;/;\n/g;s/,/,\n/g' log.txt

dsplines.svg: splines.dat Makefile
	echo "\
        set term pdf;\
        set out 'dsplines.pdf';\
        set xlabel 'x';\
        set ylabel 'y';\
        set key right bottom;\
        set tics out;\
        set title 'Fig. 1: Derivatives of aspline and cspline';\
        plot 'sinspline.dat' index 0 using 1:2 with lines title 'cspline derivative',\
	'sinspline.dat' index 0 using 1:3 with lines title 'aspline derivative',\
	'sinspline.dat' index 0 using 1:4 with lines title 'ospline derivative',\
	'sinspline.dat' index 1 using 1:2 with points title 'Analytic derivative: $df(x)=0.5 cos(0.5*x)';\
        " | tee log.txt | gnuplot
	sed -i 's/;/;\n/g;s/,/,\n/g' log.txt

d2splines.svg: splines.dat Makefile
	echo "\
        set term pdf;\
        set out 'dsplines2.pdf';\
        set xlabel 'x';\
        set ylabel 'y';\
        set key right bottom;\
        set tics out;\
        set title 'Fig. 1: Second derivatives of aspline and cspline';\
        plot 'sinspline.dat' index 0 using 1:5 with lines title 'cspline second derivative',\
	'sinspline.dat' index 0 using 1:6 with lines title 'aspline second derivative',\
	'sinspline.dat' index 2 using 1:2 with points title 'analytic second derivative: d²f(x) = -0.25*sin(0.5*x)';\
        " | tee log.txt | gnuplot
	sed -i 's/;/;\n/g;s/,/,\n/g' log.txt


splines.dat: main
	./$< > out.txt

main: main.o aspline.o cspline.o binsearch.o quadspline.o ospline.o


.PHONEY:clean

clean:
	find . -type f -executable -delete
	$(RM) *.o out* log* *.dat* *.svg *.pdf *~
