CFLAGS = -std=gnu99 $$(gsl-config --cflags) 
LDLIBS = $$(gsl-config --libs) 
CC = cc

plot.svg: result.txt Makefile
	echo '\
	set term svg size 640,480 background rgb "white";\
	set out "$@";\
	set key out;\
	set tics out;\
	set title "neural network curve fitting";\
	set xlabel "x";\
	set ylabel "y";\
	plot\
	 "$<" index 0 with points pointtype 6 title "data"\
	,"$<" index 1 with lines title "fit"\
	'|gnuplot


result.txt log.txt: main
	./$< 1> $@ 2> log.txt

main: main.o neural_network.o

.PHONEY: clean
clean:
	rm *.o main 
