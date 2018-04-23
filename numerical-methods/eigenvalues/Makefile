CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)



all: noofsweeps.txt noofsweeps.svg
	cat opgA.txt
	cat opgB.txt
	xdg-open noofsweeps.svg

noofsweeps.svg: noofsweeps.gpi noofsweeps.txt
	gnuplot noofsweeps.gpi

noofsweeps.txt: main
	./main > noofsweeps.txt

main: main.o myprint.o jac.o jac2.o

.PHONEY: clean

clean:
	rm -f main *.o out* *.png *.svg log* *log
