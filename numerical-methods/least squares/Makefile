CFLAGS += -O -march=native -Wall -std=gnu11 `gsl-config --cflags`
QRDIR = ../qr.gramschmidt.gsl
CFLAGS += -I$(QRDIR)
LDLIBS += `gsl-config --libs`

all: out.txt plot.svg
	xdg-open plot.svg
	cat opg2.txt
plot.svg: out.txt ogdata.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set tics out;\
	set key out;\
	set title "least-squares fit of F_c(x)=c_0log(x)+c_1+c_2*x to some data" ;\
	plot "ogdata.txt" index 0 with errorbars title "data"\
	,"out.txt" using 1:2 index 1 with lines title "F_c(x)"\
	,"out.txt" using 1:3 index 1 with lines title "F_c(x)_{c0/c1/c2 - dc0/dc1/dc2}"\
	,"out.txt" using 1:4 index 1  with lines title "F_c(x)_{c0/c1/c2 + dc0/dc1/dc2}"\
	'|gnuplot

out.txt: main
	./$< >out.txt

main: main.o matrixprep.o fit.o

clean:; $(RM) *.o main out* *.svg
