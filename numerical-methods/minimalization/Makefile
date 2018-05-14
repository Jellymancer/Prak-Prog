CFLAGS += -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS += `gsl-config --libs` -lm

all: out.txt rad.svg
	cat out.txt
	xdg-open rad.svg

rad.svg: out.txt Makefile
	echo '\
        set term svg size 1200,800 background rgb "white";\
        set out "rad.svg";\
        set xlabel "x";\
        set ylabel "y";\
        set tics out;\
        set key out;\
        set title "Fig1: Radioactive decay data with best fit" ;\
	plot "exp.dat" index 0 with errorbars title "Experimental data points",\
	"exp.dat" index 1 using 1:2 with lines title "Best fit";\
	'|gnuplot


out.txt: main
	./$< > $@

main:main.o newton.o print.o func.o

clean:
	rm -f main *.o out*
