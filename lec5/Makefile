CFLAGS = -Wall -std=c99 -O
LDLIBS = -lm


all: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

main: nvector.o main.o


.PHONEY:clean

clean:
	rm *.o *.txt main
