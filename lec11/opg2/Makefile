CFLAGS = -Wall -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)
PLOTTER = gnuplot



show: main final.png
	xdg-open final.png
	cat data.txt


final.png:  plot.gpi
	gnuplot plot.gpi

main: ls
	./ls > data.txt


ls: ls.o

.PHONEY: clean

clean:
	-rm *.txt *.o ls fit.dat
