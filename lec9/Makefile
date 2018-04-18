CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs` -lm


all: main qm.dat
	xdg-open qm.svg

main: integ1.dat qm.dat qm.gpi
	cat integ1.dat
	gnuplot qm.gpi

integ1.dat: integ1
	./integ1 > integ1.dat

qm.dat: QM
	./QM > qm.dat

integ1: integ1.o


QM: QM.o

.PHONEY : clean

clean:
	rm *.o *.dat *.svg QM integ1
