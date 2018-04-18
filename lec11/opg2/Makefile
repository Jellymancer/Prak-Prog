CFLAGS = -Wall -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)
PLOTTER = gnuplot



show: main final.svg
	xdg-open final.svg
	cat data.txt


final.svg:  plot.gpi
	gnuplot plot.gpi

main: ls
	./ls > data.txt


ls: ls.o

.PHONEY: clean

clean:
	-rm *.txt *.o ls fit.dat *.svg
