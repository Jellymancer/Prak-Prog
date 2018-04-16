CFLAGS = -Wall -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

all: opg1-2.txt
	cat opg1-2.txt

opg1-2.txt: main
	./main -> opg1-2.txt

main: main.o opg3i.o


.PHONEY: clean

clean:
	rm *.txt main.o
