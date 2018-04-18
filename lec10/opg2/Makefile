CFLAGS = -Wall -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS = `gsl-config --libs`
ECHO = /bin/echo


all: plot data.txt
	xdg-open final.svg

plot: data.txt plot.gpi
	gnuplot plot.gpi


data.txt: hyd.o
	gcc hyd.o -o program $(CFLAGS) $(LDLIBS)
	./program > data.txt

main: hyd.o


.PHONEY: clean

clean:
	rm program *.o *.svg data.txt
