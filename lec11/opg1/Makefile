CFLAGS = -Wall -std=gnu99 -O
CFLAGS += `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)
PLOTTER = gnuplot


show: main
	cat data.txt

main: ros
	./ros > data.txt

ros: ros.o

.PHONEY: clean

clean:
	rm *.txt *.o
