CFLAGS += -Wall -std=c11 `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

all: out.txt as1.svg plotairy.svg
	cat out.txt
	xdg-open as1.svg
	xdg-open plotairy.svg

as1.svg: out.txt Makefile
	echo '\
        set term svg size 800,600 background rgb "white";\
        set out "as1.svg";\
        set xlabel "x";\
        set ylabel "y";\
        set tics out;\
        set key out;\
        set title "Fig1: Solution to dydx=3x+4" ;\
	plot "as1ana.txt" using 1:2 title "Analytical y(x)" with points,\
	"as1.txt" using 1:2 with lines title "Numerical y(x)";\
	'|gnuplot

plotairy.svg: out.txt Makefile
	echo '\
        set term svg size 1200,800 background rgb "white";\
        set out "plotairy.svg";\
        set xlabel "x";\
        set ylabel "y";\
        set tics out;\
        set key out;\
        set title "Fig 2: Airy equation solution" ;\
        plot "airy.txt" using 1:2 with lines title "Numeric-Airy_B",\
	"GSLairy.txt" using 1:3 with points title "GSL-airy_B";\
        '|gnuplot

out.txt: main
	./main > out.txt

main: main.o print.o rkstepX.o driver.o


clean:
	rm -f main *.o out* *.svg *.txt
