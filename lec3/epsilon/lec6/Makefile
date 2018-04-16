CFLAGS = -Wall -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

all: gnuplot
	xdg-open plot.svg

gnuplot: program
	gnuplot plot.gpi

program: main-stdio main-cmdline input.data.txt
	echo `seq 0 0.5 4`| ./main-stdio > test.io.out.txt
	./main-cmdline 1 2 3 4 5 > test.cmd.out.txt
	cat input.data.txt | ./main-stdio > out.io.txt
	./main-cmdline `cat input.data.txt` >out.cmd.txt

main-stdio: main-stdio.o

main-cmdline: main-cmdline.o


input.data.txt:
	echo `seq 1 0.2 5` > input.data.txt

.PHONEY:clean

clean:
	rm *.o *.txt
