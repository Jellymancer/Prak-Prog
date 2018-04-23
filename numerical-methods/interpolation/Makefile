'CFLAGS = -Wall -std=gnu1x -O2
LDLIBS = -lm


plot: splines.svg
	xdg-open splines.svg
	cat out.txt

splines.svg: splines.dat Makefile random.dat
	echo "\
	set term svg size 800,600 font 'Times,18' fsize 18;\
	set terminal svg enhanced background rgb 'white';\
	set out '$@';\
	set xlabel 'x';\
	set ylabel 'y';\
	set key right bottom;\
	set tics out;\
	set title 'Random discrete data points with linear, cubic and quadratic splines';\
	plot 'random.dat' with points title 'data'\
	,'splines.dat' with lines title 'linear spline'\
	,'qspline.dat' using 1:2 with lines title 'quadratic spline'\
	,'qspline.dat' using 1:3 with lines title 'quadratic spline derivative'\
	#,'$<' index 0 smooth csplines title 'gnuplot cubic spline'\
	" | tee log.txt | gnuplot
	sed -i 's/;/;\n/g;s/,/,\n/g' log.txt




splines.dat: main
	./$< > out.txt

main: main.o lin.o quadspline.o


.PHONEY:clean

clean:
	find . -type f -executable -delete
	$(RM) *.o out* log* *.dat* *.svg *~
