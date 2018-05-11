CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs) -lm



all: noofsweeps.svg out1.txt out2.txt
	cat out1.txt
	cat out2.txt
	xdg-open noofsweeps.svg
	xdg-open time.svg

noofsweeps.svg: noofsweeps.gpi out2.txt time.txt
	gnuplot noofsweeps.gpi
	gnuplot time.gpi


time.txt: main1
	for n in $$(seq 50 10 300); do time --format "$$n %U" --append --output $@ ./main1 $$n 2> /dev/null; done

out1.txt: main1
	./main1 3 > out1.txt

out2.txt: main2
	./main2 > out2.txt

main1: main1.o myprint.o jac.o jac2.o

main2: main2.o myprint.o jac.o jac2.o

.PHONEY: clean

clean:
	rm -f main *.o *.txt out* *.png *.svg log* *log
