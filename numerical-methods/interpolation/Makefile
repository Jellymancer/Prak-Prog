'CFLAGS = -Wall -std=gnu1x -O2
LDLIBS = -lm


plot: splines.svg splinederiv.svg splines.dat
	xdg-open splines.svg
	xdg-open splinederiv.svg
	cat out.txt

splines.svg: splines.dat Makefile sine.dat
	echo "\
	set term svg size 800,600 font 'Times,18' fsize 18;\
	set terminal svg enhanced background rgb 'white';\
	set out '$@';\
	set xlabel 'x';\
	set ylabel 'y';\
	set key right bottom;\
	set tics out;\
	set title 'Fig. 1: Discrete data points of a sine curve with linear and quadratic splines';\
	plot 'sine.dat' using 1:2 with points title 'f(x)=sin(0.5x)'\
	,'splines.dat' with lines title 'linear spline'\
	,'qspline.dat' using 1:2 with lines title 'quadratic spline'\
	" | tee log.txt | gnuplot
	sed -i 's/;/;\n/g;s/,/,\n/g' log.txt

splinederiv.svg: splines.dat Makefile sine.dat
	echo "\
	set term svg size 800,600 font 'Times,18' fsize 18;\
	set terminal svg enhanced background rgb 'white';\
	set out '$@';\
	set xlabel 'x';\
	set ylabel 'y';\
	set key right bottom;\
	set tics out;\
	set title 'Fig. 2: Discrete data points of the derivative of sine (cosine) and the quadratic spline derivative';\
	plot 'sine.dat' using 1:3 with points title 'df(x)/dx = 0.5*cos(0.5*x)'\
	,'qspline.dat' using 1:3 with lines title 'quadratic spline derivative'\
	" | tee log.txt | gnuplot
	sed -i 's/;/;\n/g;s/,/,\n/g' log.txt




splines.dat: main
	./$< > out.txt

main: main.o lin.o quadspline.o binsearch.o


.PHONEY:clean

clean:
	find . -type f -executable -delete
	$(RM) *.o out* log* *.dat* *.svg *~
