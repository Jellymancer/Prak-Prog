CFLAGS += -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS += -lm `gsl-config --libs`

all: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

main: main.o integ.o

clean:
	rm -f main *.o out*

