CFLAGS = -Wall -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`
PLOTTER = gnuplot

all: main
	xdg-open opg1.pdf
	xdg-open opg2i.pdf
	xdg-open opg2ii.pdf
	xdg-open opg2iii.pdf

main: opg1.dat opg1.gpi opg2i.dat opg2ii.dat opg2iii.dat
	gnuplot opg1.gpi
	gnuplot opg2i.gpi
	gnuplot opg2ii.gpi
	gnuplot opg2iii.gpi

opg2i.dat: opg2i
	./opg2i > opg2i.dat

opg2ii.dat: opg2ii
	./opg2ii >opg2ii.dat

opg2iii.dat: opg2iii
	./opg2iii>opg2iii.dat



opg2i: opg2.c
	gcc opg2.c -o opg2i $(CFLAGS) $(LDLIBS)


opg2ii: opg2ii.c
	gcc opg2ii.c -o opg2ii $(CFLAGS) $(LDLIBS)

opg2iii: opg2iii.c
	gcc opg2iii.c -o opg2iii $(CFLAGS) $(LDLIBS)


opg1.dat: opg1
	./opg1 >opg1.dat

opg1: opg1.c
	gcc opg1.c -o opg1 $(CFLAGS) $(LDLIBS)



.PHONEY:clean

clean:
	rm opg1 opg2 opg2i.dat opg2ii.dat opg2iii.dat *.pdf *.dat
