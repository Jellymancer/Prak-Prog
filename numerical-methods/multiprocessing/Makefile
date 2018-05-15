CFLAGS += -Wall -std=gnu1x `gsl-config --cflags` -fopenmp
LDLIBS += -lm `gsl-config --libs` -lgomp

all: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

main: main.o integ.o

clean:
	rm -f main *.o out*

