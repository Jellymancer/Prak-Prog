'CFLAGS += -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs) -lm


all: main plot.svg plot2.svg
	xdg-open plot.svg
	xdg-open plot2.svg

plot.svg: out.txt Makefile
	echo '\
	set term svg size 1200,800 background rgb "white";\
	set out "plot.svg";\
	set xlabel "x";\
	set ylabel "y";\
	set tics out;\
	set key out;\
	set title "Fig1: Interpolating wavelet (cos5*x-1)*exp(-x²) using activation function x*exp(-x²)" ;\
	plot "nn_out.txt" index 1 with lines title "ANN interpolation",\
	"nn_out.txt" index 0 with points title "Initial known values of f(x)=cos(5*x-1)*exp(x²)";\
	'|gnuplot

plot2.svg: out.txt Makefile
	echo '\
	set term svg size 1200,800 background rgb "white";\
	set out "plot2.svg";\
	set xlabel "x";\
	set ylabel "y";\
	set tics out;\
	set key out;\
	set title "Fig2: Interpolating x³ using activation function cos(x)*exp(-x²)" ;\
	plot "nn_out2.txt" index 1 with lines title "ANN interpolation",\
	"nn_out2.txt" index 0 with points title "Initial known function values of f(x)=x³";\
	'|gnuplot


out.txt: main
	./main

main: main.o neurons.o

.PHONEY: clean
clean:
	rm main *.o *.txt
